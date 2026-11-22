"""WPC coded surface bulletins (fronts, highs, lows) → JSON.

Source: NWS Weather Prediction Center "Coded Surface Frontal Positions" (AWIPS id
CODSUS), retrieved from the Iowa Environmental Mesonet text-product archive.

Coordinates: 4 or 5 digits (low-res: 2-digit latitude, 2- or 3-digit west longitude,
whole degrees) or 7 digits (high-res: tenths of a degree). Split at len // 2, as in
MetPy's parse_wpc_surface_bulletin.
"""
from __future__ import annotations

import re
from datetime import datetime

import requests

IEM = "https://mesonet.agron.iastate.edu/cgi-bin/afos/retrieve.py"
KINDS = {"COLD": "cold", "WARM": "warm", "STNRY": "stationary", "OCFNT": "occluded", "TROF": "trough"}
STRENGTH = {"WK", "MDT", "STG"}


def decode_pos(s: str):
    if not s.isdigit() or len(s) not in (4, 5, 7):
        return None
    k = len(s) // 2
    la, lo = s[:k], s[k:]
    lat = float(la[:2] + "." + la[2:]) if len(la) > 2 else float(la)
    lon = float(lo[:3] + "." + lo[3:]) if len(lo) > 3 else float(lo)
    return [round(lat, 2), round(-lon, 2)]


def parse(text: str, ref: datetime):
    """One bulletin → {valid, highs, lows, fronts} or None."""
    m = re.search(r"VALID\s+(\d{2})(\d{2})(\d{2})Z", text)
    if not m:
        return None
    mon, day, hr = int(m.group(1)), int(m.group(2)), int(m.group(3))
    year = ref.year - (1 if mon > ref.month + 1 else 0)
    try:
        valid = datetime(year, mon, day, hr)
    except ValueError:
        return None
    body = text[m.end():].split("$$")[0]
    toks = body.split()
    out = {"valid": valid.strftime("%Y-%m-%dT%H:00Z"), "highs": [], "lows": [], "fronts": []}
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in ("HIGHS", "LOWS"):
            key = "highs" if t == "HIGHS" else "lows"
            i += 1
            while i + 1 < len(toks) and toks[i].isdigit() and len(toks[i]) in (3, 4):
                p, pos = toks[i], decode_pos(toks[i + 1])
                if pos is None:
                    break
                pr = int(p) + (1000 if int(p) < 100 else 0)
                out[key].append(pos + [pr])
                i += 2
            continue
        if t in KINDS:
            kind = KINDS[t]; i += 1
            strength = None
            if i < len(toks) and toks[i] in STRENGTH:
                strength = toks[i]; i += 1
            pts = []
            while i < len(toks) and toks[i].isdigit():
                p = decode_pos(toks[i])
                if p:
                    pts.append(p)
                i += 1
            if len(pts) >= 2:
                out["fronts"].append({"type": kind, "strength": strength, "pts": pts})
            continue
        i += 1
    return out


def fetch(ref: datetime, limit: int = 12):
    r = requests.get(IEM, params={"pil": "CODSUS", "fmt": "text", "limit": limit}, timeout=60,
                     headers={"User-Agent": "contextclimate-surface-data (github.com/josht3113)"})
    r.raise_for_status()
    parts = re.split(r"(?=CODED SURFACE FRONTAL POSITIONS)", r.text)
    seen, out = set(), []
    for p in parts:
        b = parse(p, ref)
        if b and b["valid"] not in seen and (b["fronts"] or b["highs"] or b["lows"]):
            seen.add(b["valid"]); out.append(b)
    out.sort(key=lambda b: b["valid"])
    return out


if __name__ == "__main__":
    sample = """CODED SURFACE FRONTAL POSITIONS
NWS WEATHER PREDICTION CENTER COLLEGE PARK MD
1023 AM EST SUN NOV 22 2026

VALID 112215Z
HIGHS 1032 41108 1028 35112 1030 4975
LOWS 1004 4579 1010 28102
COLD WK 4579 4478 4376 4274 4072
WARM 4579 4475 4371
STNRY 3290 3195 30100 29105
OCFNT 4880 4579
TROF 2710 2612
$$
"""
    b = parse(sample, datetime(2026, 11, 22, 16))
    import json
    print(json.dumps(b, indent=1))
    assert b["valid"] == "2026-11-22T15:00Z"
    assert b["highs"][0] == [41.0, -108.0, 1032]
    assert b["highs"][2] == [49.0, -75.0, 1030]
    assert b["lows"][1] == [28.0, -102.0, 1010]
    assert b["fronts"][0]["strength"] == "WK" and b["fronts"][0]["pts"][0] == [45.0, -79.0]
    assert b["fronts"][2]["pts"][2] == [30.0, -100.0]
    hr = parse("VALID 112218Z\nCOLD 4101080 4051072\n$$", datetime(2026, 11, 22, 20))
    assert hr["fronts"][0]["pts"] == [[41.0, -108.0], [40.5, -107.2]], hr
    print("ok")
