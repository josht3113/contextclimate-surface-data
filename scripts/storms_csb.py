"""WPC coded surface bulletins → monthly JSON files for the storm-track builder.

Source: NWS Weather Prediction Center "Coded Surface Frontal Positions" (AWIPS id
CODSUS), every 3 hours, from the Iowa Environmental Mesonet text-product archive.
Two versions go out for each analysis: whole degrees (WMO ASUS01 / ASUS1) and, from
2011 in this archive, tenths of a degree (ASUS02). The tenths version wins; among
copies of the same version, the latest issuance (a correction) wins.

Month file  bulletins/YYYY/MM.json.gz:
  {"v": 1, "month": "YYYY-MM", "b": [
     {"t": "YYYYMMDDHH", "hr": 0|1,
      "lows": [[lat, lon, mb], ...],
      "fronts": [[kind, [lat, lon, lat, lon, ...]], ...]}, ...]}
  kind: c cold, w warm, s stationary, o occluded, t trough. Longitudes are negative (°W).
"""
from __future__ import annotations

import gzip
import json
import os
import re
import time
from datetime import datetime, timedelta, timezone

import requests

IEM = "https://mesonet.agron.iastate.edu/cgi-bin/afos/retrieve.py"
UA = {"User-Agent": "contextclimate-surface-data (github.com/josht3113)"}
KINDS = {"COLD": "c", "WARM": "w", "STNRY": "s", "OCFNT": "o", "TROF": "t"}
STRENGTH = {"WK", "MDT", "STG"}
KEYWORDS = set(KINDS) | {"HIGHS", "LOWS"}
FIRST_YEAR = 2003


def decode_pos(s: str):
    """4/5 digits = whole degrees (2-digit lat + 2/3-digit °W lon); 7 digits = tenths."""
    if not s.isdigit() or len(s) not in (4, 5, 7):
        return None
    if len(s) == 7:
        lat, lon = int(s[:3]) / 10, int(s[3:]) / 10
    else:
        lat, lon = int(s[:2]), int(s[2:])
    if not (0 <= lat <= 89.9 and 0 <= lon <= 360):
        return None
    if lon > 180:                       # east longitudes come through as >180 °W
        lon -= 360
    return [round(lat, 1), round(-lon, 1)]


def valid_mb(s: str):
    if not s.isdigit():
        return None
    v = int(s)
    if len(s) == 3 and 870 <= v <= 999:
        return v
    if len(s) == 4 and 1000 <= v <= 1065:
        return v
    return None


def parse_one(text: str, ref: datetime):
    """One bulletin's text → dict or None. ref: any time within a month of the valid time."""
    m = re.search(r"VALID\s+(\d{2})(\d{2})(\d{2})Z", text)
    if not m:
        return None
    mon, day, hr = int(m.group(1)), int(m.group(2)), int(m.group(3))
    year = ref.year
    if mon - ref.month > 6:
        year -= 1
    elif ref.month - mon > 6:
        year += 1
    try:
        valid = datetime(year, mon, day, hr)
    except ValueError:
        return None
    if hr % 3:
        return None
    body = text[m.end():].split("$$")[0]
    toks = body.split()
    lows, fronts, n7 = [], [], 0
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in ("HIGHS", "LOWS"):
            i += 1
            j = i
            while j < len(toks) and toks[j].isdigit():
                j += 1
            seq = toks[i:j]
            k = 0
            while k + 1 < len(seq):
                p, pos = valid_mb(seq[k]), decode_pos(seq[k + 1])
                if p is not None and pos is not None:
                    if t == "LOWS":
                        lows.append(pos + [p])
                    if len(seq[k + 1]) == 7:
                        n7 += 1
                    k += 2
                else:
                    k += 1                 # analyst typo (missing or garbled field): resync
            i = j
            continue
        if t in KINDS:
            kind = KINDS[t]
            i += 1
            if i < len(toks) and toks[i] in STRENGTH:
                i += 1
            pts = []
            while i < len(toks) and toks[i].isdigit():
                p = decode_pos(toks[i])
                if p:
                    pts.extend(p)
                    if len(toks[i]) == 7:
                        n7 += 1
                i += 1
            if len(pts) >= 4:
                fronts.append([kind, pts])
            continue
        i += 1
    if not lows and not fronts:
        return None
    return {"t": valid.strftime("%Y%m%d%H"), "hr": 1 if n7 else 0, "lows": lows, "fronts": fronts}


HEADER = re.compile(r"\b(ASUS\d{1,2})\s+KWBC\s+(\d{6})")


def split_products(text: str):
    """IEM text dump → [(issue DDHHMM or '', product text)]."""
    out = []
    pieces = re.split(r"(?=\bASUS\d{1,2}\s+KWBC\s+\d{6})", text)
    for p in pieces:
        h = HEADER.search(p[:80])
        if h:
            out.append((h.group(2), p))
        elif "VALID" in p:
            out.append(("", p))
    return out


def dedupe(products, ref: datetime):
    """Keep one bulletin per valid time: tenths version first, then latest issuance."""
    best = {}
    for issue, text in products:
        b = parse_one(text, ref)
        if not b:
            continue
        key = (b["hr"], issue, len(b["lows"]) + len(b["fronts"]))
        cur = best.get(b["t"])
        if cur is None or key > cur[0]:
            best[b["t"]] = (key, b)
    return [best[t][1] for t in sorted(best)]


def fetch_range(s: datetime, e: datetime, tries: int = 4):
    """Every CODSUS product issued in [s, e) (UTC) as text."""
    params = {"pil": "CODSUS", "fmt": "text", "limit": 9999,
              "sdate": s.strftime("%Y-%m-%dT%H:%MZ"), "edate": e.strftime("%Y-%m-%dT%H:%MZ")}
    err = None
    for k in range(tries):
        if k:
            time.sleep(5 * k)
        try:
            r = requests.get(IEM, params=params, headers=UA, timeout=120)
            if r.status_code == 200:
                return r.text
            err = RuntimeError(f"IEM HTTP {r.status_code}")
        except requests.RequestException as x:
            err = x
    raise err


def month_bounds(y: int, m: int):
    s = datetime(y, m, 1)
    e = datetime(y + (m == 12), m % 12 + 1, 1)
    return s, e


def build_month(y: int, m: int):
    """Fetch and parse one month (products issued up to 6 h into the next month are included,
    so the month's last analyses are never lost to a late issuance)."""
    s, e = month_bounds(y, m)
    text = fetch_range(s, e + timedelta(hours=6))
    bs = dedupe(split_products(text), s + timedelta(days=14))
    lo, hi = s.strftime("%Y%m%d%H"), e.strftime("%Y%m%d%H")
    return [b for b in bs if lo <= b["t"] < hi]


def month_path(root: str, y: int, m: int):
    return os.path.join(root, f"{y:04d}", f"{m:02d}.json.gz")


def write_month(root: str, y: int, m: int, bs, complete: bool):
    p = month_path(root, y, m)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    doc = {"v": 1, "month": f"{y:04d}-{m:02d}", "complete": complete, "b": bs}
    with gzip.open(p + ".tmp", "wt") as f:
        json.dump(doc, f, separators=(",", ":"))
    os.replace(p + ".tmp", p)


def read_month(root: str, y: int, m: int):
    p = month_path(root, y, m)
    if not os.path.exists(p):
        return None
    with gzip.open(p, "rt") as f:
        return json.load(f)


def sync(root: str, now: datetime | None = None, first_year: int = FIRST_YEAR, log=print):
    """Make the month files complete through now. Past months are fetched once (and again only
    if marked incomplete); the current month and the one before are refreshed every call."""
    now = now or datetime.now(timezone.utc).replace(tzinfo=None)
    y, m = first_year, 1
    fetched = 0
    while (y, m) <= (now.year, now.month):
        s, e = month_bounds(y, m)
        recent = e > now - timedelta(days=3)
        doc = read_month(root, y, m)
        if doc is None or not doc.get("complete") or recent:
            bs = build_month(y, m)
            if not bs and doc is None and not recent:
                log(f"  {y}-{m:02d}: no bulletins")
            write_month(root, y, m, bs, complete=not recent)
            fetched += 1
            time.sleep(0.5)
        m += 1
        if m == 13:
            y, m = y + 1, 1
    return fetched


if __name__ == "__main__":
    sample = """
ASUS02 KWBC 152100
CODSUS
CODED SURFACE FRONTAL POSITIONS
VALID 011521Z
HIGHS 1029 4191318 1038 4081137
LOWS 993 5930859 996 6310726 2797 1008 3221115
COLD 3070699 3030716
$$
ASUS01 KWBC 152221
CODSUS
VALID 011521Z
LOWS 993 5986 996 6373 1008 32111
COLD 3170 3072
$$
ASUS1 KWBC 152327
CODSUS
VALID 011518Z
LOWS 1009 36100 1006 5379 664 7652 996 53115 1005 3581 52135 1009 40105
STNRY WK 37102 37103 38104
$$
"""
    bs = dedupe(split_products(sample), datetime(2025, 1, 15))
    assert [b["t"] for b in bs] == ["2025011518", "2025011521"], bs
    b18, b21 = bs
    assert b21["hr"] == 1 and b21["lows"] == [[59.3, -85.9, 993], [63.1, -72.6, 996], [32.2, -111.5, 1008]], b21
    assert b21["fronts"] == [["c", [30.7, -69.9, 30.3, -71.6]]]
    assert b18["lows"] == [[36, -100, 1009], [53, -79, 1006], [53, -115, 996], [35, -81, 1005], [40, -105, 1009]], b18
    assert b18["fronts"][0][0] == "s"
    print("ok")
