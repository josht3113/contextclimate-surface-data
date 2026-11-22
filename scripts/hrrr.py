"""HRRR surface fields → compact regular-grid frames for the ContextClimate
Surface Analysis Builder.

Source: NOAA HRRR CONUS `wrfsfc` GRIB2 on the NOAA Open Data AWS bucket. Only the
needed messages are downloaded (HTTP byte ranges from the .idx files).

Frame file layout (one file per valid time):
  int16 little-endian planes, order = VARS, each plane [lat][lon] on GRID
  (lat ascending from GRID south, lon ascending from GRID west, edges included).
  Bytes are shuffled (all low bytes, then all high bytes), then gzip.
  value = stored * scale + offset; −32768 = no data (outside the HRRR domain).
"""
from __future__ import annotations

import gzip
import time
from datetime import datetime, timedelta, timezone

import numpy as np
import requests

BUCKET = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
UA = {"User-Agent": "contextclimate-surface-data (github.com/josht3113)"}

GRID = {"west": -126.0, "east": -65.0, "south": 22.0, "north": 51.0, "res": 0.125}
GRID["nx"] = int(round((GRID["east"] - GRID["west"]) / GRID["res"])) + 1
GRID["ny"] = int(round((GRID["north"] - GRID["south"]) / GRID["res"])) + 1

# name, (GRIB short name, level string), scale, offset, unit
VARS = [
    ("t",    ("TMP",   "2 m above ground"),   0.1,  0.0,    "°F"),
    ("td",   ("DPT",   "2 m above ground"),   0.1,  0.0,    "°F"),
    ("u",    ("UGRD",  "10 m above ground"),  0.01, 0.0,    "m/s"),
    ("v",    ("VGRD",  "10 m above ground"),  0.01, 0.0,    "m/s"),
    ("p",    ("MSLMA", "mean sea level"),     0.1,  1000.0, "mb"),
    ("refc", ("REFC",  "entire atmosphere"),  0.5,  0.0,    "dBZ"),
]
MISSING = -32768

_s = requests.Session()
_s.headers.update(UA)


def utcnow() -> datetime:
    return datetime.now(timezone.utc).replace(tzinfo=None)


def tkey(t: datetime) -> str:
    return t.strftime("%Y%m%d%H")


def parse_tkey(s: str) -> datetime:
    return datetime.strptime(s, "%Y%m%d%H")


def file_url(cycle: datetime, fh: int) -> str:
    return f"{BUCKET}/hrrr.{cycle:%Y%m%d}/conus/hrrr.t{cycle:%H}z.wrfsfcf{fh:02d}.grib2"


def get(url: str, headers: dict | None = None, tries: int = 4, timeout: int = 90):
    """GET with retries; None when the file isn't posted (403/404)."""
    for k in range(tries):
        try:
            r = _s.get(url, headers=headers or {}, timeout=timeout)
            if r.status_code in (403, 404):
                return None
            if r.status_code in (200, 206):
                return r
        except requests.RequestException:
            pass
        time.sleep(2 * (k + 1))
    raise RuntimeError(f"GET failed: {url}")


def idx(cycle: datetime, fh: int):
    r = get(file_url(cycle, fh) + ".idx")
    if r is None:
        return None
    rows = [ln.split(":") for ln in r.text.strip().splitlines()]
    out = []
    for i, p in enumerate(rows):
        start = int(p[1])
        end = int(rows[i + 1][1]) - 1 if i + 1 < len(rows) else -1
        out.append((start, end, p[3], p[4], p[5] if len(p) > 5 else ""))
    return out


def posted(cycle: datetime, fh: int) -> bool:
    return get(file_url(cycle, fh) + ".idx", tries=2) is not None


def pick_messages(rows):
    """(start, end) for each wanted field: the instantaneous message (no ave/max/acc)."""
    want = {}
    for start, end, var, lev, fc in rows:
        for name, (gv, gl), *_ in VARS:
            if var == gv and lev == gl and name not in want and not any(w in fc for w in ("ave", "max", "min", "acc")):
                want[name] = (start, end)
    missing = [v[0] for v in VARS if v[0] not in want]
    if missing:
        raise RuntimeError(f"idx missing fields: {missing}")
    return want


def download(cycle: datetime, fh: int):
    rows = idx(cycle, fh)
    if rows is None:
        return None
    msgs = {}
    for name, (a, b) in pick_messages(rows).items():
        rng = f"bytes={a}-" if b < 0 else f"bytes={a}-{b}"
        r = get(file_url(cycle, fh), headers={"Range": rng})
        if r is None:
            return None
        msgs[name] = r.content
    return msgs


# ── decoding and regridding ────────────────────────────────
class Regridder:
    """Mean of the source points within ~1 target cell of each target point
    (16 nearest, distance-capped), on the unit sphere."""

    def __init__(self, lat_src: np.ndarray, lon_src: np.ndarray):
        from scipy.spatial import cKDTree
        lon_t = GRID["west"] + GRID["res"] * np.arange(GRID["nx"])
        lat_t = GRID["south"] + GRID["res"] * np.arange(GRID["ny"])
        LON, LAT = np.meshgrid(lon_t, lat_t)
        xyz = lambda la, lo: np.stack([np.cos(np.radians(la)) * np.cos(np.radians(lo)),
                                       np.cos(np.radians(la)) * np.sin(np.radians(lo)),
                                       np.sin(np.radians(la))], axis=-1)
        tree = cKDTree(xyz(lat_src.ravel(), lon_src.ravel()))
        d, i = tree.query(xyz(LAT.ravel(), LON.ravel()), k=16)
        cell = np.radians(GRID["res"]) * 0.75          # radius ≈ 0.75 target cells
        self.w = (d <= cell).astype(np.float32)
        self.i = i
        self.n = self.w.sum(axis=1)
        self.shape = LAT.shape

    def __call__(self, field: np.ndarray) -> np.ndarray:
        v = field.ravel()[self.i]
        ok = np.isfinite(v) & (self.w > 0)
        s = np.where(ok, v, 0).sum(axis=1)
        n = ok.sum(axis=1)
        out = np.where(n > 0, s / np.maximum(n, 1), np.nan)
        return out.reshape(self.shape).astype(np.float32)


def decode(buf: bytes):
    import eccodes
    gid = eccodes.codes_new_from_message(buf)
    try:
        nx = eccodes.codes_get(gid, "Nx")
        ny = eccodes.codes_get(gid, "Ny")
        vals = eccodes.codes_get_values(gid).astype(np.float32)
        miss = eccodes.codes_get(gid, "missingValue")
        meta = {
            "lov": eccodes.codes_get(gid, "LoVInDegrees"),
            "latin1": eccodes.codes_get(gid, "Latin1InDegrees"),
            "uvRelativeToGrid": eccodes.codes_get(gid, "uvRelativeToGrid") if eccodes.codes_is_defined(gid, "uvRelativeToGrid") else 0,
        }
        lat = lon = None
        if not hasattr(decode, "_ll"):
            lat = eccodes.codes_get_array(gid, "latitudes").astype(np.float64)
            lon = eccodes.codes_get_array(gid, "longitudes").astype(np.float64)
    finally:
        eccodes.codes_release(gid)
    vals[vals == miss] = np.nan
    a = vals.reshape(ny, nx)
    if lat is not None:
        lon = np.where(lon > 180, lon - 360, lon)
        decode._ll = (lat.reshape(ny, nx), lon.reshape(ny, nx))
    return a, meta


def earth_winds(u, v, lon, meta):
    """Rotate Lambert-conformal grid-relative winds to earth-relative (NCEP convention)."""
    if not meta["uvRelativeToGrid"]:
        return u, v
    lov = meta["lov"] - 360 if meta["lov"] > 180 else meta["lov"]
    a = np.radians(np.sin(np.radians(meta["latin1"])) * (lon - lov))
    return np.cos(a) * u + np.sin(a) * v, -np.sin(a) * u + np.cos(a) * v


def build_frame(msgs, regrid_cache: dict):
    fields, metas = {}, {}
    for name, buf in msgs.items():
        fields[name], metas[name] = decode(buf)
    lat, lon = decode._ll
    if "rg" not in regrid_cache:
        regrid_cache["rg"] = Regridder(lat, lon)
    rg = regrid_cache["rg"]
    u, v = earth_winds(fields["u"], fields["v"], lon, metas["u"])
    k2f = lambda k: (k - 273.15) * 9 / 5 + 32
    out = {
        "t": rg(k2f(fields["t"])), "td": rg(k2f(fields["td"])),
        "u": rg(u), "v": rg(v),
        "p": rg(fields["p"] / 100.0),
        "refc": rg(np.maximum(fields["refc"], -10.0)),
    }
    return out


def encode(frame: dict) -> bytes:
    planes = []
    for name, _, scale, offset, _ in VARS:
        a = frame[name]
        q = np.where(np.isfinite(a), np.round((a - offset) / scale), MISSING)
        planes.append(np.clip(q, MISSING, 32767).astype("<i2"))
    raw = np.stack(planes).tobytes()
    b = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 2)
    return gzip.compress(np.concatenate([b[:, 0], b[:, 1]]).tobytes(), compresslevel=9, mtime=0)


def decode_frame(blob: bytes) -> dict:
    raw = np.frombuffer(gzip.decompress(blob), dtype=np.uint8)
    n = raw.size // 2
    b = np.empty((n, 2), np.uint8); b[:, 0] = raw[:n]; b[:, 1] = raw[n:]
    a = b.view("<i2").reshape(len(VARS), GRID["ny"], GRID["nx"]).astype(np.float32)
    out = {}
    for k, (name, _, scale, offset, _) in enumerate(VARS):
        p = a[k]
        out[name] = np.where(p == MISSING, np.nan, p * scale + offset)
    return out


def sample(frame: dict, name: str, lat: float, lon: float) -> float:
    x = (lon - GRID["west"]) / GRID["res"]; y = (lat - GRID["south"]) / GRID["res"]
    i, j = int(np.floor(x)), int(np.floor(y))
    if i < 0 or j < 0 or i >= GRID["nx"] - 1 or j >= GRID["ny"] - 1:
        return float("nan")
    tx, ty = x - i, y - j
    f = frame[name]
    return float((f[j, i] * (1 - tx) + f[j, i + 1] * tx) * (1 - ty) + (f[j + 1, i] * (1 - tx) + f[j + 1, i + 1] * tx) * ty)
