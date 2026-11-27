"""Storm tracks for the ContextClimate Storm Tracks map.

  python scripts/storms_build.py --data data            # hourly run (skips when nothing new)
  python scripts/storms_build.py --data data --force    # rebuild every year

Inputs (kept in the Actions cache with the rest of data/, refetched when missing):
  data/storms/bulletins/YYYY/MM.json.gz   WPC coded surface bulletins (storms_csb.py)
  data/storms/side/                       IBTrACS (NA + EP), CPC ONI / NAO / PNA

Published (data/site/storms/):
  index.json            years, latest analysis, climate indices, method constants
  tracks/YYYY.json.gz   storms that formed in YYYY
  b/YYYY/MM.json.gz     the month's bulletins (fronts + lows) at 00 and 12 UTC; the two most
                        recent months keep every 3-hourly analysis

Storm record (all times are whole hours since 1970-01-01 UTC):
  id   "YYYYMMDDHH-k" (genesis analysis, k-th storm that formed then)
  t    genesis time
  ty   P Pacific · A Alberta clipper · C Colorado/lee · G Gulf · E East Coast · T tropical · O other
  f    flags: n nor'easter (Oct–Apr, up or out along the coast, within 600 km of 40°N 70°W) · b bomb · j redevelopment jump · q quasi-stationary or heat low ·
       m merged into another low · a still active · y active and younger than 24 h ·
       s significant (see SIG_*) · u passed within 300 km of the lower 48
  mi   Miller type for nor'easters ("A"/"B"), when it can be told
  tc   tropical cyclone name (IBTrACS), when the storm was or became one
  pm   lowest pressure (mb), tm its hour offset
  bg   largest 24-h deepening in Bergerons (1 = 24 mb · sin φ / sin 60°), dp the 24-h fall in mb
  bm   nor'easters: [signed km from 40°N 70°W (+ inside / northwest, − outside), hour offset]
  L    path length (km)
  h    hour offsets of the main-path points; la, lo: tenths of a degree; p: mb
  j    point indices that follow a jump
  bw   bomb windows [[i0, i1], ...] (point indices, 24 h apart)
  sd   faint pieces [{h, la, lo, p}]: the primary after a jump, the secondary before it
  x    point-index ranges [[i0, i1], ...] taken from IBTrACS (tropical / post-tropical positions)
  d    depth (mb) at each point: 750 km ring-mean pressure minus central pressure (-99 = unknown); dm its maximum
  Positions are smoothed (storms_smooth.py); pressures are WPC's.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import storms_csb as csb                      # noqa: E402
import storms_track as trk                    # noqa: E402
import storms_smooth as smo                   # noqa: E402
import storms_depth as dep                    # noqa: E402

VERSION = 1
FIRST_YEAR = csb.FIRST_YEAR
EPOCH = datetime(1970, 1, 1)
BENCH = (40.0, -70.0)
OVER_KM = 80.0               # within 50 mi = "over the benchmark"
UA = csb.UA

IBTRACS = ("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/"
           "v04r01/access/csv/")
IB_FULL = ("ibtracs.NA.list.v04r01.csv", "ibtracs.EP.list.v04r01.csv")
IB_RECENT = "ibtracs.last3years.list.v04r01.csv"
CPC = {
    "oni": "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
    "nao": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii",
    "pna": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii",
}
ONI_SEASONS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]

# significant storm: the deepest low within PRIMARY_KM, at least SIG_DEPTH mb below the ring mean,
# for SIG_HOURS in all (12 h for a storm still going); bombs, nor'easters and tropical cyclones always
PRIMARY_KM = 800.0
SIG_DEPTH = 8.0
SIG_HOURS = 24
SIG_HOURS_ACTIVE = 12
US_MASK = None


def us48(lat, lon):
    """Within 300 km of the lower 48 (precomputed 0.25° mask)."""
    global US_MASK
    if US_MASK is None:
        import base64
        m = json.loads(Path(os.path.dirname(os.path.abspath(__file__)), "us48_mask.json").read_text())
        bits = np.unpackbits(np.frombuffer(base64.b64decode(m["bits"]), np.uint8))[: m["ny"] * m["nx"]]
        US_MASK = (m, bits.reshape(m["ny"], m["nx"]).astype(bool))
    m, g = US_MASK
    i, j = int((lat - m["lat0"]) / m["res"]), int((lon - m["lon0"]) / m["res"])
    return 0 <= i < m["ny"] and 0 <= j < m["nx"] and bool(g[i, j])


TYPES = [
    ["P", "Pacific"], ["A", "Alberta clipper"], ["C", "Colorado / lee"], ["G", "Gulf"],
    ["E", "East Coast"], ["T", "Tropical"], ["O", "Other"],
]


def log(*a):
    print(*a, flush=True)


def hkey(t: int) -> str:
    return (EPOCH + timedelta(hours=t)).strftime("%Y%m%d%H")


def thours(key: str) -> int:
    d = datetime(int(key[:4]), int(key[4:6]), int(key[6:8]), int(key[8:10]))
    return int((d - EPOCH).total_seconds() // 3600)


def month_of(t: int) -> int:
    return (EPOCH + timedelta(hours=t)).month


def year_of(t: int) -> int:
    return (EPOCH + timedelta(hours=t)).year


# ── side data ────────────────────────────────────────────────────────
def fetch(url, timeout=300):
    err = None
    for k in range(3):
        if k:
            time.sleep(5 * k)
        try:
            r = requests.get(url, headers=UA, timeout=timeout)
            if r.status_code == 200:
                return r.text
            err = RuntimeError(f"HTTP {r.status_code} {url}")
        except requests.RequestException as x:
            err = x
    raise err


OFFLINE = False
DEPTH_EXTRA = None          # dev only: {(key, lat10, lon10): depth} for analyses without a cached GFS field


def refresh(path: Path, url: str, max_age_h: float):
    """Cached download; a failed refresh keeps the old copy."""
    if path.exists() and (OFFLINE or time.time() - path.stat().st_mtime < max_age_h * 3600):
        return path.read_text()
    try:
        txt = fetch(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(txt)
        return txt
    except Exception as e:                    # noqa: BLE001
        if path.exists():
            log(f"  keep cached {path.name}: {e}")
            return path.read_text()
        log(f"  no {path.name}: {e}")
        return None


def ibtracs_rows(side: Path):
    """{hour: [(lat, lon, mb, name, sid)]} for tropical, subtropical and post-tropical points
    since FIRST_YEAR (3-hourly, North Atlantic and eastern North Pacific)."""
    out = {}
    files = [(f, 24 * 30) for f in IB_FULL] + [(IB_RECENT, 24)]
    for name, age in files:
        p = side / (name + ".slim")
        if not (p.exists() and (OFFLINE or time.time() - p.stat().st_mtime < age * 3600)):
            try:
                raw = fetch(IBTRACS + name, timeout=600)
                rd = csv.reader(io.StringIO(raw))
                hdr = next(rd); next(rd)
                ix = {h: i for i, h in enumerate(hdr)}
                rows = []
                for r in rd:
                    if r[ix["SEASON"]] < str(FIRST_YEAR - 1) or r[ix["BASIN"]] not in ("NA", "EP"):
                        continue
                    rows.append([r[ix["SID"]], r[ix["NAME"]], r[ix["ISO_TIME"]], r[ix["NATURE"]],
                                 r[ix["LAT"]], r[ix["LON"]], r[ix["TRACK_TYPE"]], r[ix["USA_PRES"]], r[ix["WMO_PRES"]]])
                buf = io.StringIO(); csv.writer(buf).writerows(rows)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(buf.getvalue())
            except Exception as e:            # noqa: BLE001
                log(f"  IBTrACS {name}: {e}")
        if not p.exists():
            continue
        for row in csv.reader(io.StringIO(p.read_text())):
            sid, nm, iso, nature, lat, lon, ttype, up, wp = row
            add_tc(out, sid, nm, iso, nature, lat, lon, ttype, up or wp)
    return tc_by_time(out)


def add_tc(out, sid, nm, iso, nature, lat, lon, ttype, pres):
    if nature not in ("TS", "SS", "ET", "MX") and not (nature == "NR" and ttype.startswith("PROV")):
        return
    try:
        d = datetime.strptime(iso, "%Y-%m-%d %H:%M:%S")
        la, lo, pr = float(lat), float(lon), float(pres)
    except ValueError:
        return
    if d.minute or d.hour % 3 or not (850 < pr < 1030):
        return
    t = int((d - EPOCH).total_seconds() // 3600)
    out[(sid, t)] = (t, la, lo, pr, nm.title() if nm and nm not in ("NOT_NAMED", "UNNAMED") else "", sid)


def tc_by_time(out):
    by_t = defaultdict(list)
    for t, la, lo, pr, nm, sid in out.values():
        by_t[t].append((la, lo, pr, nm, sid))
    return by_t


def inject_tc(frames, tc_by_t):
    """Tropical cyclones are not in WPC's list of lows: add IBTrACS positions as lows (src = storm
    number), and tag a WPC low within 300 km as that storm. → frames with a 4th column, names."""
    sids, names = {}, [""]
    out = []
    for t, lows in frames:
        lows = np.hstack([lows, np.zeros((len(lows), 1))]) if len(lows) else np.zeros((0, 4))
        extra = []
        for la, lo, pr, nm, sid in tc_by_t.get(t, ()):
            k = sids.get(sid)
            if k is None:
                k = sids[sid] = len(names); names.append(nm)
            elif nm and not names[k]:
                names[k] = nm
            if len(lows):
                d = trk.hav(la, lo, lows[:, 0], lows[:, 1])
                i = int(d.argmin())
                if d[i] <= 300 and lows[i, 3] == 0:
                    lows[i, 3] = k
                    continue
            extra.append([la, lo, pr, k])
        if extra:
            lows = np.vstack([lows, np.array(extra)])
        out.append((t, lows))
    return out, names


def parse_oni(txt):
    out = []
    for line in (txt or "").splitlines():
        p = line.split()
        if len(p) == 4 and p[0] in ONI_SEASONS:
            try:
                out.append([int(p[1]), ONI_SEASONS.index(p[0]) + 1, round(float(p[3]), 2)])   # centre month
            except ValueError:
                pass
    return out


def parse_monthly(txt):
    out = []
    for line in (txt or "").splitlines():
        p = line.split()
        if len(p) == 3:
            try:
                y, m, v = int(p[0]), int(p[1]), float(p[2])
            except ValueError:
                continue
            if y >= FIRST_YEAR - 1 and -9 < v < 9:
                out.append([y, m, round(v, 2)])
    return out


# ── geography ────────────────────────────────────────────────────────
def east_coast(lat, lon):
    if not (25 <= lat <= 46 and lon <= -60):
        return False
    if lat < 34:
        return lon > -82
    if lat < 38:
        return lon > -79.5
    if lat < 42:
        return lon > -76
    return lon > -72


def genesis_type(lat, lon):
    if lon < -124 or lon > 0:
        return "P"
    if 46 <= lat <= 60 and -124 <= lon <= -100:
        return "A"
    if 30 <= lat < 46 and -114 <= lon <= -98:
        return "C"
    if 18 <= lat <= 31 and -98 <= lon <= -80:
        return "G"
    if east_coast(lat, lon):
        return "E"
    return "O"


def in_ne_box(lat, lon):
    return 34 <= lat <= 46 and -78 <= lon <= -60


def benchmark(points):
    """Closest approach of the main path to 40°N 70°W: (signed km, hour offset) — + inside (NW)."""
    c = math.cos(math.radians(BENCH[0])) * 111.2
    xy = [((lo - BENCH[1]) * c, (la - BENCH[0]) * 111.2, q[0]) for q in points for la, lo in [(q[1], q[2])]]
    best = None
    for k in range(len(xy) - 1):
        if points[k + 1][4]:                       # never across a jump
            continue
        (x0, y0, t0), (x1, y1, t1) = xy[k], xy[k + 1]
        dx, dy = x1 - x0, y1 - y0
        L2 = dx * dx + dy * dy
        f = 0.0 if L2 == 0 else max(0.0, min(1.0, -(x0 * dx + y0 * dy) / L2))
        x, y = x0 + f * dx, y0 + f * dy
        d = math.hypot(x, y)
        if best is None or d < best[0]:
            best = (d, x, y, t0 + f * (t1 - t0))
    if best is None:
        return None
    d, x, y, t = best
    side = (-x + y) / math.sqrt(2)                 # + toward the coast (northwest)
    return round(d if side >= 0 else -d), t


def bergeron(points):
    """[(i, j, N, dp)] for every pair 24 h apart within a jump-free stretch."""
    out = []
    times = [p[0] for p in points]
    seg = [0] * len(points)
    for k in range(1, len(points)):
        seg[k] = seg[k - 1] + (1 if points[k][4] else 0)
    for i, (t, la, lo, p, *_) in enumerate(points):
        j = None
        for jj in range(i + 1, len(points)):
            if times[jj] - t >= 24:
                if times[jj] - t <= 27:
                    j = jj
                break
        if j is None or seg[j] != seg[i]:
            continue
        dt = times[j] - t
        dp = (p - points[j][3]) * 24 / dt
        phi = (la + points[j][1]) / 2
        if phi < 15:
            continue
        n = dp * math.sin(math.radians(60)) / (24 * math.sin(math.radians(phi)))
        out.append((i, j, n, dp))
    return out


# ── build ────────────────────────────────────────────────────────────
def load_frames(broot: Path):
    frames, months = [], []
    for ydir in sorted(broot.glob("[0-9]" * 4)):
        for f in sorted(ydir.glob("*.json.gz")):
            with gzip.open(f, "rt") as fh:
                doc = json.load(fh)
            months.append((doc["month"], doc))
            for b in doc["b"]:
                lows = np.array(b["lows"], float).reshape(-1, 3)
                frames.append((thours(b["t"]), lows))
    frames.sort(key=lambda f: f[0])
    return frames, months


def primary(lows, lat, lon, p):
    """No deeper low within PRIMARY_KM at the same analysis."""
    if lows is None or not len(lows):
        return True
    deeper = lows[:, 2] < p - 0.01
    if not deeper.any():
        return True
    d = trk.hav(lat, lon, lows[deeper, 0], lows[deeper, 1])
    return bool((d > PRIMARY_KM).all())


def analyse(storms, latest_t, tc_names, depth=None, lows_at=None):
    out = []
    end_index = defaultdict(list)                  # hour → [(storm k, lat, lon)] for merge checks
    pts_all, deps, prims = [], [], []
    for s in storms:
        pts = s.points()
        dd, pp = [], []
        for t, la, lo, p, _j, _src in pts:
            key = hkey(t)
            dd.append(depth(key, la, lo, p) if depth else None)
            pp.append(primary(lows_at.get(t) if lows_at else None, la, lo, p))
        deps.append(dd); prims.append(pp)
        coarse = all(float(p[1]).is_integer() and float(p[2]).is_integer() for p in pts if not p[5])
        pts_all.append(smo.smooth_points(pts, coarse))
    for k, pts in enumerate(pts_all):
        for t, la, lo, p, *_ in pts:
            end_index[t].append((k, la, lo))
    for k, (s, pts) in enumerate(zip(storms, pts_all)):
        if len(pts) < 2:
            continue
        dd, pp = deps[k], prims[k]
        t0, t1 = pts[0][0], pts[-1][0]
        life = t1 - t0
        active = t1 >= latest_t - trk.MAX_COAST_H + 3
        if life < trk.MIN_LIFE_H and not (active and life >= trk.YOUNG_LIFE_H):
            continue
        la = np.array([p[1] for p in pts]); lo = np.array([p[2] for p in pts]); pr = np.array([p[3] for p in pts])
        flags = set()
        if active:
            flags.add("a")
            if life < trk.MIN_LIFE_H:
                flags.add("y")
        # quasi-stationary / heat lows
        cla, clo = float(la.mean()), float(lo.mean())
        rmax = float(trk.hav(cla, clo, la, lo).max())
        warm = month_of(t0) in (4, 5, 6, 7, 8, 9, 10)
        if rmax < 400 or (warm and 20 <= cla <= 42 and -125 <= clo <= -100 and rmax < 700):
            flags.add("q")
        L = float(trk.hav(la[:-1], lo[:-1], la[1:], lo[1:]).sum()) if len(pts) > 1 else 0.0
        jumps = [i for i, p in enumerate(pts) if p[4]]
        # legs across a jump are not travelled distance
        for i in jumps:
            L -= float(trk.hav(la[i - 1], lo[i - 1], la[i], lo[i]))
        if jumps:
            flags.add("j")
        # tropical: positions that came from IBTrACS (or a WPC low tagged as that cyclone)
        tc_pts = Counter(p[5] for p in pts if p[5])
        tc = None
        if tc_pts:
            k, c = tc_pts.most_common(1)[0]
            if c >= 2:
                tc = tc_names[k] or "Unnamed"
        ty = "T" if tc else genesis_type(pts[0][1], pts[0][2])
        # nor'easter: a cold-season, non-tropical low moving up or out along the Northeast coast
        mi = ""
        box = [i for i, p in enumerate(pts) if in_ne_box(p[1], p[2])]
        if ty != "T" and len(box) >= 2 and month_of(pts[box[0]][0]) in (10, 11, 12, 1, 2, 3, 4):
            vx = vy = 0.0
            for i in box:
                if i + 1 < len(pts) and (i + 1) in box and not pts[i + 1][4]:
                    vx += (pts[i + 1][2] - pts[i][2]) * math.cos(math.radians(pts[i][1]))
                    vy += pts[i + 1][1] - pts[i][1]
            v = math.hypot(vx, vy)
            bmk = benchmark(pts)
            if v > 0 and vy >= 0.1 * v and vx >= -0.3 * v and bmk and abs(bmk[0]) <= 600:
                flags.add("n")
                # Miller B: the primary (Ohio Valley / Great Lakes) hands over to a coastal low;
                # Miller A: one low from the Gulf or the Southeast coast
                pre = [i for i in jumps if i <= box[0] + 1]
                if pre and pts[pre[-1] - 1][1] >= 36 and pts[pre[-1] - 1][2] < -78 and (
                        east_coast(pts[pre[-1]][1], pts[pre[-1]][2]) or in_ne_box(pts[pre[-1]][1], pts[pre[-1]][2])):
                    mi = "B"
                elif pts[0][1] < 36 and pts[0][2] > -100 and all(pts[i - 1][1] < 36 for i in pre):
                    mi = "A"
        # bombs
        bw, bg, dp = [], None, None
        for i, j, n, d in bergeron(pts):
            if bg is None or n > bg:
                bg, dp = n, d
            if n >= 1:
                if bw and i <= bw[-1][1]:
                    bw[-1][1] = max(bw[-1][1], j)
                else:
                    bw.append([i, j])
        if bw:
            flags.add("b")
        # merged: ends close to another low that carries on
        if not active:
            for kk, a, b in end_index.get(t1 + 3, []) + end_index.get(t1 + 6, []):
                if kk != k and pts_all[kk][-1][0] > t1 + 6 and trk.hav(pts[-1][1], pts[-1][2], a, b) <= 400:
                    flags.add("m")
                    break
        # depth, significance, lower-48 scope
        known = [d for d in dd if d is not None]
        deep_h = 3 * sum(1 for d, q in zip(dd, pp) if d is not None and q and d >= SIG_DEPTH)
        if deep_h >= (SIG_HOURS_ACTIVE if active else SIG_HOURS) or "b" in flags or "n" in flags or ty == "T":
            flags.add("s")
        if any(us48(p[1], p[2]) for p in pts):
            flags.add("u")
        im = int(pr.argmin())
        rec = {
            "t": t0, "ty": ty, "f": "".join(sorted(flags)),
            "pm": int(pr[im]), "tm": pts[im][0] - t0,
            "L": int(round(L)),
            "h": [p[0] - t0 for p in pts],
            "la": [int(round(p[1] * 10)) for p in pts],
            "lo": [int(round(p[2] * 10)) for p in pts],
            "p": [int(p[3]) for p in pts],
        }
        if known:
            rec["d"] = [-99 if d is None else int(round(d)) for d in dd]
            rec["dm"] = int(round(max(known)))
        if mi:
            rec["mi"] = mi
        if tc:
            rec["tc"] = tc
        if bg is not None:
            rec["bg"] = round(bg, 2); rec["dp"] = int(round(dp))
        if "n" in flags:
            b = benchmark(pts)
            if b:
                rec["bm"] = [b[0], int(round(b[1] - t0))]
        if jumps:
            rec["j"] = jumps
        xs = [i for i, p in enumerate(pts) if p[5]]
        if xs:
            rng = [[xs[0], xs[0]]]
            for i in xs[1:]:
                if i == rng[-1][1] + 1:
                    rng[-1][1] = i
                else:
                    rng.append([i, i])
            rec["x"] = rng
        if bw:
            rec["bw"] = bw
        sd = []
        for seg in s.side_points():
            sd.append({"h": [q[0] - t0 for q in seg], "la": [int(round(q[1] * 10)) for q in seg],
                       "lo": [int(round(q[2] * 10)) for q in seg], "p": [int(q[3]) for q in seg]})
        if sd:
            rec["sd"] = sd
        out.append(rec)
    out.sort(key=lambda r: (r["t"], r["la"][0], r["lo"][0]))
    seen = Counter()
    for r in out:
        key = hkey(r["t"])
        r["id"] = f"{key}-{seen[key]}"
        seen[key] += 1
    return out


def write_gz(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wt", compresslevel=9) as f:
        json.dump(obj, f, separators=(",", ":"))
    os.replace(tmp, path)


def code_hash():
    h = hashlib.sha1()
    for f in ("storms_csb.py", "storms_track.py", "storms_smooth.py", "storms_depth.py", "storms_build.py"):
        h.update(Path(os.path.dirname(os.path.abspath(__file__)), f).read_bytes())
    return h.hexdigest()[:12]


def run(data: Path, force=False, now=None, offline=False):
    t_start = time.time()
    now = now or datetime.now(timezone.utc).replace(tzinfo=None)
    broot, side, site = data / "storms" / "bulletins", data / "storms" / "side", data / "site" / "storms"
    global OFFLINE
    OFFLINE = offline
    if not offline:
        log("storms: syncing bulletins")
        n = csb.sync(str(broot), now=now, log=log)
        log(f"  {n} month files refreshed ({time.time() - t_start:.0f} s)")
    frames, months = load_frames(broot)
    if not frames:
        raise SystemExit("storms: no bulletins")
    latest_t = frames[-1][0]
    depth = dep.Depth(side / "gfs", [hkey(t) for t, _ in frames if year_of(t) >= FIRST_YEAR], log=log, fetch=not offline)
    if DEPTH_EXTRA:
        depth.extra = DEPTH_EXTRA
    stamp = {"latest": hkey(latest_t), "code": code_hash(), "n": len(frames), "gfs": len(depth.fields)}
    prev = site / "index.json"
    if not force and prev.exists():
        try:
            old = json.loads(prev.read_text())
            if old.get("stamp") == stamp:
                log("storms: nothing new")
                return
        except Exception:                     # noqa: BLE001
            pass
    tc_by_t = ibtracs_rows(side)
    oni = parse_oni(refresh(side / "oni.txt", CPC["oni"], 24))
    nao = parse_monthly(refresh(side / "nao.txt", CPC["nao"], 24))
    pna = parse_monthly(refresh(side / "pna.txt", CPC["pna"], 24))

    log(f"storms: linking {len(frames)} analyses")
    frames = [f for f in frames if year_of(f[0]) >= FIRST_YEAR]
    frames, tc_names = inject_tc(frames, tc_by_t)
    tracks = trk.link(frames)
    storms = trk.build_storms(tracks)
    lows_at = {t: lows for t, lows in frames}
    recs = analyse(storms, latest_t, tc_names, depth, lows_at)
    log(f"  {len(tracks)} tracklets → {len(recs)} storms ({time.time() - t_start:.0f} s)")

    by_year = defaultdict(list)
    for r in recs:
        by_year[year_of(r["t"])].append(r)
    years = {}
    bul_year = Counter(); hr_year = Counter()
    for _m, doc in months:
        y = int(doc["month"][:4])
        bul_year[y] += len(doc["b"]); hr_year[y] += sum(b["hr"] for b in doc["b"])
    tmp = site.with_name("storms.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    for y in sorted(set(bul_year) | set(by_year)):
        if y < FIRST_YEAR:
            continue
        rs = by_year.get(y, [])
        write_gz(tmp / "tracks" / f"{y}.json.gz", {"v": VERSION, "year": y, "storms": rs})
        years[str(y)] = {"n": len(rs), "bulletins": bul_year[y], "hires": round(hr_year[y] / max(1, bul_year[y]), 3)}
    # bulletin months: 00/12 UTC, every analysis in the two newest months
    newest = [m for m, _ in months][-2:]
    for m, doc in months:
        y = int(m[:4])
        if y < FIRST_YEAR:
            continue
        keep = doc["b"] if m in newest else [b for b in doc["b"] if b["t"][8:10] in ("00", "12")]
        write_gz(tmp / "b" / m[:4] / f"{m[5:]}.json.gz", {"v": 1, "month": m, "b": keep})
    index = {
        "v": VERSION,
        "generated": now.strftime("%Y-%m-%dT%H:%MZ"),
        "latest": hkey(latest_t),
        "first": hkey(frames[0][0]) if year_of(frames[0][0]) >= FIRST_YEAR else f"{FIRST_YEAR}010100",
        "stamp": stamp,
        "years": years,
        "types": TYPES,
        "method": {
            "minLifeH": trk.MIN_LIFE_H, "gateKmh": trk.MAX_SPEED, "jumpKm": list(trk.JUMP_KM),
            "bench": list(BENCH), "overKm": OVER_KM,
            "ringKm": dep.RING_KM, "sigDepth": SIG_DEPTH, "sigHours": SIG_HOURS, "primaryKm": PRIMARY_KM,
        },
        "indices": {"oni": oni, "nao": nao, "pna": pna},
        "bulletinMonths": [m for m, _ in months if int(m[:4]) >= FIRST_YEAR],
        "fullMonths": newest,
    }
    (tmp / "index.json").write_text(json.dumps(index, separators=(",", ":")))
    shutil.rmtree(site, ignore_errors=True)
    os.replace(tmp, site)
    log(f"storms: wrote {len(recs)} storms, {len(years)} years ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--offline", action="store_true", help="use the bulletins and side files already in --data")
    a = ap.parse_args()
    run(Path(a.data), force=a.force, offline=a.offline)
