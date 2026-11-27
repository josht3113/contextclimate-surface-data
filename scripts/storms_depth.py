"""How deep each WPC low is relative to its surroundings.

depth = mean sea-level pressure on a 750 km ring around the low − WPC's central pressure.

Fields:
  - NCEP/NCAR Reanalysis 1, 6-hourly, 2.5° (NOAA PSL), 2003 → 17 Mar 2026, when it stopped.
    Frozen, so the depths are precomputed once: static/storm-depth-r1.csv.gz
    (t, lat×10, lon×10, depth×10) — regenerate with `python scripts/storms_depth.py --r1 OUT`.
  - GFS analyses, 1°, PRMSL only (byte ranges from the NOAA GFS archive on AWS), after that.
    Each analysis is kept, cropped to the northern hemisphere, in the Actions cache.
    GFS reads 0.6 mb deeper than R1 over the Jan–Mar 2026 overlap (r = 0.99); that is removed.
Analyses at 03/09/15/21 UTC use the mean of the fields either side.
"""
from __future__ import annotations

import gzip
import io
import math
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import requests

RING_KM = 750.0
GFS_BIAS = 0.6
R1_END = "2026031718"
HERE = Path(os.path.dirname(os.path.abspath(__file__)))
R1_TABLE = HERE.parent / "static" / "storm-depth-r1.csv.gz"
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
UA = {"User-Agent": "contextclimate-surface-data (github.com/josht3113)"}
AZ = np.radians(np.arange(0, 360, 15))


def ring(lat, lon, r_km=RING_KM):
    d = r_km / 6371.0
    la, lo = math.radians(lat), math.radians(lon)
    la2 = np.arcsin(np.sin(la) * np.cos(d) + np.cos(la) * np.sin(d) * np.cos(AZ))
    lo2 = lo + np.arctan2(np.sin(AZ) * np.sin(d) * np.cos(la), np.cos(d) - np.sin(la) * np.sin(la2))
    return np.degrees(la2), (np.degrees(lo2) + 360) % 360


def sample(field, lat0, dlat, dlon, la, lo):
    """Bilinear sample; field rows run from lat0 southward by dlat, columns from 0°E by dlon."""
    y = (lat0 - la) / dlat
    x = lo / dlon
    ny, nx = field.shape
    y0 = np.clip(np.floor(y).astype(int), 0, ny - 2)
    x0 = np.floor(x).astype(int) % nx
    x1 = (x0 + 1) % nx
    fy = np.clip(y - y0, 0, 1)
    fx = x - np.floor(x)
    return (field[y0, x0] * (1 - fx) * (1 - fy) + field[y0, x1] * fx * (1 - fy)
            + field[y0 + 1, x0] * (1 - fx) * fy + field[y0 + 1, x1] * fx * fy)


def ring_mean(field, lat0, dlat, dlon, lat, lon):
    a, b = ring(lat, lon)
    return float(np.mean(sample(field, lat0, dlat, dlon, a, b)))


# ── R1 (frozen table) ────────────────────────────────────────────────
def load_r1():
    out = {}
    if not R1_TABLE.exists():
        return out
    with gzip.open(R1_TABLE, "rt") as f:
        next(f)
        for line in f:
            t, la, lo, d = line.split(",")
            out[(t, int(la), int(lo))] = int(d) / 10.0
    return out


# ── GFS ──────────────────────────────────────────────────────────────
def gfs_path(cache: Path, key: str):
    return cache / f"{key}.npy"


def fetch_gfs(key: str, cache: Path):
    """PRMSL (mb, float16, 0–90°N × 360) for analysis key YYYYMMDDHH, cached. None if not out yet."""
    p = gfs_path(cache, key)
    if p.exists():
        return np.load(p).astype(np.float32)
    ymd, hh = key[:8], key[8:10]
    base = f"{GFS}/gfs.{ymd}/{hh}/atmos/gfs.t{hh}z.pgrb2.1p00.f000"
    for k in range(3):
        if k:
            time.sleep(3 * k)
        try:
            r = requests.get(base + ".idx", headers=UA, timeout=60)
            if r.status_code in (403, 404):
                return None
            r.raise_for_status()
            lines = r.text.strip().split("\n")
            for i, l in enumerate(lines):
                if ":PRMSL:mean sea level:" in l:
                    a = int(l.split(":")[1])
                    b = int(lines[i + 1].split(":")[1]) - 1 if i + 1 < len(lines) else ""
                    g = requests.get(base, headers={**UA, "Range": f"bytes={a}-{b}"}, timeout=60)
                    g.raise_for_status()
                    import pygrib                       # noqa: PLC0415
                    tmp = cache / f".{key}.grb2"
                    tmp.write_bytes(g.content)
                    v = pygrib.open(str(tmp)).read(1)[0].values / 100.0
                    tmp.unlink()
                    nh = np.asarray(v[:91], dtype=np.float16)     # rows 90°N → 0°
                    cache.mkdir(parents=True, exist_ok=True)
                    np.save(p, nh)
                    return nh.astype(np.float32)
            return None
        except Exception:                               # noqa: BLE001
            continue
    return None


def six_hourly(key):
    d = datetime.strptime(key, "%Y%m%d%H")
    if d.hour % 6 == 0:
        return [key]
    return [(d - timedelta(hours=3)).strftime("%Y%m%d%H"), (d + timedelta(hours=3)).strftime("%Y%m%d%H")]


class Depth:
    """depth(t_key, lat, lon, p) → mb or None."""

    def __init__(self, cache: Path, keys_needed, log=print, fetch=True):
        self.r1 = load_r1()
        self.cache = cache
        self.fields = {}
        new = [k for k in sorted({k6 for k in keys_needed if k > R1_END for k6 in six_hourly(k)})]
        got = miss = 0
        for k in new:
            if fetch or gfs_path(cache, k).exists():
                f = fetch_gfs(k, cache) if fetch else np.load(gfs_path(cache, k)).astype(np.float32)
            else:
                f = None
            if f is None:
                miss += 1
            else:
                self.fields[k] = f; got += 1
        log(f"  depth: R1 table {len(self.r1)} lows; GFS analyses {got} ready, {miss} not available")

    def field(self, key):
        ks = six_hourly(key)
        fs = [self.fields.get(k) for k in ks]
        if any(f is None for f in fs):
            return None
        return fs[0] if len(fs) == 1 else (fs[0] + fs[1]) / 2

    extra = None

    def __call__(self, key, lat, lon, p):
        if self.extra is not None and key > R1_END:
            v = self.extra.get((key, int(round(lat * 10)), int(round(lon * 10))))
            if v is not None:
                return v
        if key <= R1_END:
            return self.r1.get((key, int(round(lat * 10)), int(round(lon * 10))))
        f = self.field(key)
        if f is None:
            return None
        return ring_mean(f, 90.0, 1.0, 1.0, lat, lon) - p - GFS_BIAS


# ── one-off: rebuild the R1 table ────────────────────────────────────
def build_r1_table(bulletin_root: Path, out: Path):
    """Needs netCDF4. Downloads each year of NCEP R1 slp.YYYY.nc (≈19 MB) from NOAA PSL."""
    import glob
    import json
    import netCDF4                                       # noqa: PLC0415
    frames = []
    for f in sorted(glob.glob(str(bulletin_root / "*" / "*.json.gz"))):
        for b in json.load(gzip.open(f, "rt"))["b"]:
            if "2003010100" <= b["t"] <= R1_END:
                frames.append((b["t"], b["lows"]))
    w = gzip.open(out, "wt", compresslevel=9)
    w.write("t,la,lo,d\n")
    cur, F, T0 = None, None, None
    for key, lows in frames:
        d = datetime.strptime(key, "%Y%m%d%H")
        if d.year != cur:
            cur = d.year
            url = f"https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/slp.{cur}.nc"
            r = requests.get(url, headers=UA, timeout=600); r.raise_for_status()
            fn = f"/tmp/r1_{cur}.nc"; Path(fn).write_bytes(r.content)
            ds = netCDF4.Dataset(fn)
            F = np.array(ds.variables["slp"][:], dtype=np.float32) / 100.0
            tt = ds.variables["time"]; t0 = netCDF4.num2date(tt[0], tt.units)
            T0 = datetime(t0.year, t0.month, t0.day, t0.hour)
            ds.close(); os.remove(fn)
        i = (d - T0).total_seconds() / 3600 / 6
        if i > F.shape[0] - 1:
            continue
        i0 = int(math.floor(i)); f = i - i0
        fld = F[i0] if f == 0 else F[i0] * (1 - f) + F[min(i0 + 1, F.shape[0] - 1)] * f
        for la, lo, p in lows:
            dep = ring_mean(fld, 90.0, 2.5, 2.5, la, (lo + 360) % 360) - p
            w.write(f"{key},{int(round(la * 10))},{int(round(lo * 10))},{int(round(dep * 10))}\n")
    w.close()


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--r1":
        build_r1_table(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        print("usage: storms_depth.py --r1 BULLETIN_ROOT OUT.csv.gz")
