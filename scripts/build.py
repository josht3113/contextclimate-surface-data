"""Hourly build: HRRR analyses for the past 24 h, the latest complete HRRR run's
18-hour forecast, WPC coded surface bulletins, and a check against METARs.

  python scripts/build.py --data data [--max-new 30]

`data/site` is the published folder and is kept between runs in the Actions cache,
so each run only downloads what is new.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
import time
from datetime import timedelta
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).parent))
import fronts  # noqa: E402
import hrrr  # noqa: E402

WINDOW_H = 24
FCST_H = 18
ROSTER = "https://raw.githubusercontent.com/josht3113/context-climate/main/public/station-roster.js"
IEM_ASOS = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"


def log(*a):
    print(*a, flush=True)


def latest_cycle(now):
    c = now.replace(minute=0, second=0, microsecond=0)
    for _ in range(5):
        if hrrr.posted(c, 0):
            return c
        c -= timedelta(hours=1)
    return None


def write_frame(path: Path, cycle, fh, cache):
    t0 = time.time()
    msgs = hrrr.download(cycle, fh)
    if msgs is None:
        return False
    frame = hrrr.build_frame(msgs, cache)
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = hrrr.encode(frame)
    path.write_bytes(blob)
    log(f"  {path.relative_to(path.parents[2])}: {len(blob) / 1024:.0f} KB in {time.time() - t0:.1f} s")
    return True


# ── validation against routine METARs ──────────────────────
def roster():
    txt = requests.get(ROSTER, timeout=60).text
    rows = re.findall(r'\["([A-Z0-9]{4})"\s*,\s*"[^"]*"\s*,\s*"[A-Z]{2}"\s*,\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)', txt)
    return [(r[0], float(r[1]), float(r[2])) for r in rows]


def validate(site: Path, valid):
    try:
        st = [s for s in roster() if hrrr.GRID["south"] < s[1] < hrrr.GRID["north"] and hrrr.GRID["west"] < s[2] < hrrr.GRID["east"]]
    except Exception as e:  # noqa: BLE001
        log(f"validation skipped (roster): {e}")
        return None
    s0, s1 = valid - timedelta(minutes=20), valid + timedelta(minutes=5)
    params = [("station", s[0]) for s in st] + [("data", d) for d in ("tmpf", "dwpf", "drct", "sknt", "mslp")] + [
        ("year1", s0.year), ("month1", s0.month), ("day1", s0.day), ("hour1", s0.hour), ("minute1", s0.minute),
        ("year2", s1.year), ("month2", s1.month), ("day2", s1.day), ("hour2", s1.hour), ("minute2", s1.minute),
        ("tz", "Etc/UTC"), ("format", "onlycomma"), ("missing", "empty"), ("report_type", "3")]
    try:
        txt = requests.get(IEM_ASOS, params=params, timeout=120).text
    except Exception as e:  # noqa: BLE001
        log(f"validation skipped (IEM): {e}")
        return None
    lines = txt.strip().splitlines()
    head = lines[0].split(",")
    obs = {}
    for ln in lines[1:]:
        p = dict(zip(head, ln.split(",")))
        obs[p["station"]] = p                     # last report in the window
    frame = hrrr.decode_frame((site / "a" / f"{hrrr.tkey(valid)}.bin.gz").read_bytes())
    pos = {s[0][1:] if s[0][0] == "K" else s[0]: s for s in st}
    dT, dD, dP, dDir, dDirAlt = [], [], [], [], []
    lov, cone = -97.5, math.sin(math.radians(38.5))
    for k, p in obs.items():
        if k not in pos:
            continue
        _, la, lo = pos[k]
        f = lambda n: hrrr.sample(frame, n, la, lo)
        g = lambda c: float(p[c]) if p.get(c) not in (None, "", "M") else None
        if g("tmpf") is not None and math.isfinite(f("t")):
            dT.append(f("t") - g("tmpf"))
        if g("dwpf") is not None and math.isfinite(f("td")):
            dD.append(f("td") - g("dwpf"))
        if g("mslp") is not None and math.isfinite(f("p")):
            dP.append(f("p") - g("mslp"))
        u, v = f("u"), f("v")
        if g("sknt") and g("sknt") >= 8 and g("drct") is not None and math.isfinite(u) and math.hypot(u, v) > 3.5:
            md = (math.degrees(math.atan2(-u, -v)) + 360) % 360
            dDir.append(abs((md - g("drct") + 540) % 360 - 180))
            a = -2 * math.radians(cone * (lo - lov))          # winds rotated the other way
            ua, va = math.cos(a) * u + math.sin(a) * v, -math.sin(a) * u + math.cos(a) * v
            mda = (math.degrees(math.atan2(-ua, -va)) + 360) % 360
            dDirAlt.append(abs((mda - g("drct") + 540) % 360 - 180))
    st_ = lambda d: {"n": len(d), "bias": round(float(np.mean(d)), 2), "rmse": round(float(np.sqrt(np.mean(np.square(d)))), 2)} if d else None
    res = {"valid": valid.strftime("%Y-%m-%dT%H:00Z"), "temp_F": st_(dT), "dew_F": st_(dD), "mslp_mb": st_(dP),
           "wind_dir_mae_deg": round(float(np.mean(dDir)), 1) if dDir else None,
           "wind_dir_mae_deg_if_rotated_other_way": round(float(np.mean(dDirAlt)), 1) if dDirAlt else None, "wind_n": len(dDir)}
    log("validation vs METARs:", json.dumps(res))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--max-new", type=int, default=30)
    a = ap.parse_args()
    site = Path(a.data) / "site"
    (site / "a").mkdir(parents=True, exist_ok=True)
    (site / "f").mkdir(parents=True, exist_ok=True)
    now = hrrr.utcnow()
    C = latest_cycle(now)
    if C is None:
        sys.exit("no HRRR analysis posted in the last 5 hours")
    log(f"latest HRRR analysis: {C:%Y-%m-%d %HZ}")
    cache = {}

    # analyses: one f00 per hour, newest first so a short run still publishes the present
    new = 0
    for k in range(0, WINDOW_H + 1):
        H = C - timedelta(hours=k)
        path = site / "a" / f"{hrrr.tkey(H)}.bin.gz"
        if path.exists() or new >= a.max_new:
            continue
        if write_frame(path, H, 0, cache):
            new += 1
    keep = {hrrr.tkey(C - timedelta(hours=k)) for k in range(WINDOW_H + 1)}
    for p in (site / "a").glob("*.bin.gz"):
        if p.name.split(".")[0] not in keep:
            p.unlink()
    analyses = sorted(p.name.split(".")[0] for p in (site / "a").glob("*.bin.gz"))

    # forecast: the newest run with every hour posted (HRRR posts hours out of order)
    fc = None
    for k in range(0, 4):
        c = C - timedelta(hours=k)
        if all(hrrr.posted(c, fh) for fh in range(FCST_H, 0, -1)):
            fc = c
            break
    if fc:
        fdir = site / "f" / hrrr.tkey(fc)
        todo = [fh for fh in range(1, FCST_H + 1) if not (fdir / f"{fh:02d}.bin.gz").exists()]
        if todo:
            log(f"forecast run {fc:%Y-%m-%d %HZ}")
        for fh in todo:
            if not write_frame(fdir / f"{fh:02d}.bin.gz", fc, fh, cache):
                log(f"  f{fh:02d} not available yet; will resume next run")
                break
    complete = lambda d: all((d / f"{fh:02d}.bin.gz").exists() for fh in range(1, FCST_H + 1))
    fdirs = sorted(d for d in (site / "f").iterdir() if complete(d))
    newest = fdirs[-1].name if fdirs else ""
    for d in list((site / "f").iterdir()):
        if d.name < newest or (d.name > newest and fc and d.name != hrrr.tkey(fc)):
            shutil.rmtree(d, ignore_errors=True)
    fdirs = fdirs[-1:]
    forecast = {"cycle": fdirs[-1].name, "hours": list(range(1, FCST_H + 1))} if fdirs else None

    # WPC bulletins in the window
    try:
        bl = fronts.fetch(now, limit=14)
        t0 = (C - timedelta(hours=WINDOW_H + 3)).strftime("%Y-%m-%dT%H:00Z")
        bl = [b for b in bl if b["valid"] >= t0]
        (site / "fronts.json").write_text(json.dumps({"source": "NWS Weather Prediction Center coded surface bulletins (CODSUS)", "bulletins": bl}, separators=(",", ":")))
        log(f"WPC bulletins: {len(bl)} ({', '.join(b['valid'] for b in bl)})")
    except Exception as e:  # noqa: BLE001
        log(f"WPC bulletins unavailable: {e}")

    val = validate(site, C) if (site / "a" / f"{hrrr.tkey(C)}.bin.gz").exists() else None
    old = {}
    if (site / "manifest.json").exists():
        try:
            old = json.loads((site / "manifest.json").read_text())
        except ValueError:
            pass
    man = {
        "updated": now.strftime("%Y-%m-%dT%H:%MZ"),
        "source": "NOAA High-Resolution Rapid Refresh (HRRR) CONUS, NOAA Open Data Dissemination (AWS)",
        "grid": {k: hrrr.GRID[k] for k in ("west", "south", "res", "nx", "ny")},
        "vars": [{"name": n, "grib": g[0], "level": g[1], "scale": s, "offset": o, "unit": u} for n, g, s, o, u in hrrr.VARS],
        "missing": hrrr.MISSING,
        "analyses": analyses,
        "forecast": forecast,
        "fronts": "fronts.json" if (site / "fronts.json").exists() else None,
        "validation": val or old.get("validation"),
    }
    (site / "manifest.json").write_text(json.dumps(man, indent=1))
    log(f"manifest: {len(analyses)} analyses, forecast {forecast['cycle'] if forecast else 'none'}")
    if len(analyses) < 2:
        sys.exit("too few analyses built")


if __name__ == "__main__":
    main()
