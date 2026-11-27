"""Low-pressure tracks from WPC coded surface bulletins.

Every low WPC analysed is linked to the next analysis (3 h later) by solving an
assignment problem; a track may miss one analysis. Short-lived tracks are dropped.
A secondary low that takes over from a dying primary east of the Rockies (e.g. a
Miller B coastal transfer) continues the same storm, with a jump in the path.

Storm fields written by storms_build.py; see that file for the published format.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment

R_EARTH = 6371.0
STEP_H = 3
MAX_COAST_H = 6            # a track may skip one analysis
MAX_SPEED = 100.0          # km/h for the gate (fast clippers run 80–100)
GATE_MIN = 250.0           # km: allows for analysis jitter and whole-degree positions
P_WEIGHT = 12.0            # km per mb of pressure change in the link cost
MIN_LIFE_H = 24            # a storm must last this long
YOUNG_LIFE_H = 12          # an ongoing low is shown once it has lasted this long

# redevelopment (jump): the secondary forms while the primary still exists (or just after),
# downstream (bearing NE → S), within this distance, and becomes as deep as the primary
JUMP_KM = (200.0, 900.0)
JUMP_BEARING = (45.0, 170.0)
JUMP_LEAD_H = 18           # the primary dies within this long after the secondary forms
JUMP_GAP_H = 9             # or the secondary forms this soon after the primary ends
JUMP_TX_H = 18             # the hand-over happens within this long of the secondary forming
JUMP_WEST_LIMIT = -105.0   # no jumps across the Rockies
JUMP_DEEPEN = 4            # mb the secondary deepens in its first 24 h
JUMP_PRIMARY_MB = 1012     # the primary must be at least this deep
CONT_KM = 300.0            # a new low this close to a just-ended track is the same low


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.minimum(1, np.sqrt(a)))


def bearing(lat1, lon1, lat2, lon2):
    p1, p2, dl = math.radians(lat1), math.radians(lat2), math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


@dataclass
class Tracklet:
    id: int
    t: list = field(default_factory=list)      # hours since 1970
    lat: list = field(default_factory=list)
    lon: list = field(default_factory=list)
    p: list = field(default_factory=list)
    src: list = field(default_factory=list)    # 0 = WPC low, k > 0 = tropical cyclone k (IBTrACS)

    @property
    def t0(self): return self.t[0]

    @property
    def t1(self): return self.t[-1]

    def at(self, t):
        """Position/pressure at time t (linear between analyses), or None outside the track."""
        if t < self.t[0] or t > self.t[-1]:
            return None
        i = int(np.searchsorted(self.t, t))
        if self.t[i] == t:
            return self.lat[i], self.lon[i], self.p[i]
        f = (t - self.t[i - 1]) / (self.t[i] - self.t[i - 1])
        return (self.lat[i - 1] + f * (self.lat[i] - self.lat[i - 1]),
                self.lon[i - 1] + f * (self.lon[i] - self.lon[i - 1]),
                self.p[i - 1] + f * (self.p[i] - self.p[i - 1]))

    def index(self, t):
        i = int(np.searchsorted(self.t, t))
        return i if i < len(self.t) and self.t[i] == t else None


def _add(tr, t, row):
    tr.t.append(t); tr.lat.append(float(row[0])); tr.lon.append(float(row[1])); tr.p.append(float(row[2]))
    tr.src.append(int(row[3]) if len(row) > 3 else 0)


def link(frames):
    """frames: [(t_hours, ndarray[n,3|4] lat/lon/p[/src])] sorted by time → list[Tracklet]."""
    done, active, nid = [], [], 0
    for t, lows in frames:
        # retire tracks that have missed too many analyses
        keep = []
        for tr in active:
            (keep if t - tr.t1 <= MAX_COAST_H else done).append(tr)
        active = keep
        n = len(lows)
        used = np.zeros(n, bool)
        if active and n:
            m = len(active)
            cost = np.full((m, n), 1e9)
            for a, tr in enumerate(active):
                dt = t - tr.t1
                la, lo = tr.lat[-1], tr.lon[-1]
                if len(tr.t) >= 2:
                    k = len(tr.t) - 2
                    span = tr.t1 - tr.t[k]
                    vla = (tr.lat[-1] - tr.lat[k]) / span
                    vlo = (tr.lon[-1] - tr.lon[k]) / span
                    pla, plo = la + 0.8 * vla * dt, lo + 0.8 * vlo * dt
                else:
                    pla, plo = la, lo
                d_raw = hav(la, lo, lows[:, 0], lows[:, 1])
                d_pred = hav(pla, plo, lows[:, 0], lows[:, 1])
                gate = max(GATE_MIN, MAX_SPEED * dt)
                c = np.minimum(d_pred, d_raw * 1.15) + P_WEIGHT * np.abs(lows[:, 2] - tr.p[-1]) / max(1, dt / STEP_H)
                if lows.shape[1] > 3:            # a tropical cyclone keeps its own identity
                    src = lows[:, 3]
                    c = c + np.where((src > 0) & (tr.src[-1] > 0) & (src != tr.src[-1]), 1e9, 0)
                ok = (d_raw <= gate) & (c <= gate * 1.3)
                cost[a, ok] = c[ok]
            ri, ci = linear_sum_assignment(cost)
            for a, b in zip(ri, ci):
                if cost[a, b] >= 1e9:
                    continue
                _add(active[a], t, lows[b])
                used[b] = True
        for b in np.where(~used)[0]:
            tr = Tracklet(nid); nid += 1
            _add(tr, t, lows[b])
            active.append(tr)
    return done + active


@dataclass
class Storm:
    parts: list                    # [(Tracklet, t_from, t_to, join)] main path; join None | "cont" | "jump"
    side: list = field(default_factory=list)   # [(Tracklet, t_from, t_to)] faint pieces (primary remnant, early secondary)

    def points(self):
        """[[t, lat, lon, p, jump_before]] along the main path, strictly increasing in time."""
        pts = []
        for tr, a, b, join in self.parts:
            first = True
            for i, tt in enumerate(tr.t):
                if a <= tt <= b and (not pts or tt > pts[-1][0]):
                    pts.append([tt, tr.lat[i], tr.lon[i], tr.p[i], join == "jump" and first, tr.src[i]])
                    first = False
        return pts

    def side_points(self):
        out = []
        for tr, a, b in self.side:
            seg = [[tt, tr.lat[i], tr.lon[i], tr.p[i]] for i, tt in enumerate(tr.t) if a <= tt <= b]
            if len(seg) >= 2:
                out.append(seg)
        return out


def build_storms(tracks):
    """Tracklets → storms: continuations and redevelopment jumps joined, then filtered."""
    tracks = [t for t in tracks if len(t.t) >= 2]
    by_start = sorted(tracks, key=lambda t: t.t0)
    starts = np.array([t.t0 for t in by_start])

    # 1. continuations: a low that appears close to where a track ended, within the coast window
    #    plus one more analysis (the gate missed it) — join outright
    nxt = {}
    taken = set()
    for A in sorted(tracks, key=lambda t: t.t1):
        lo, hi = np.searchsorted(starts, A.t1 + STEP_H), np.searchsorted(starts, A.t1 + 9, side="right")
        best, bd = None, CONT_KM
        for B in by_start[lo:hi]:
            if B.id in taken:
                continue
            d = float(hav(A.lat[-1], A.lon[-1], B.lat[0], B.lon[0]))
            if d < bd and abs(B.p[0] - A.p[-1]) <= 6 and not (A.src[-1] and B.src[0] and A.src[-1] != B.src[0]):
                best, bd = B, d
        if best is not None:
            nxt[A.id] = ("cont", best, A.t1)
            taken.add(best.id)

    # 2. redevelopment jumps (candidates found through an index on end times)
    by_end = sorted(tracks, key=lambda t: t.t1)
    ends = np.array([t.t1 for t in by_end])
    cands = []
    for B in tracks:
        if B.id in taken or B.t1 - B.t0 < 12 or B.lat[0] < 24 or B.src[0]:
            continue
        k36 = int(np.searchsorted(B.t, B.t0 + 36, side="right"))
        pB = min(B.p[:max(1, k36)])
        k24 = int(np.searchsorted(B.t, B.t0 + 24, side="right"))
        if B.p[0] - min(B.p[:max(1, k24)]) < JUMP_DEEPEN:
            continue                         # a secondary that never deepens is a frontal wave
        lo, hi = np.searchsorted(ends, B.t0 - JUMP_GAP_H), np.searchsorted(ends, B.t0 + JUMP_LEAD_H, side="right")
        for A in by_end[lo:hi]:
            if A is B or A.id in nxt or A.t0 > B.t0 - 12 or A.src[-1]:
                continue
            pa = A.at(B.t0) or (A.lat[-1], A.lon[-1], A.p[-1])
            if pa[1] < JUMP_WEST_LIMIT or pB > pa[2] or pa[2] > JUMP_PRIMARY_MB:
                continue
            if abs(pa[0] - B.lat[0]) > 9 or abs(pa[1] - B.lon[0]) > 14:
                continue
            d = float(hav(pa[0], pa[1], B.lat[0], B.lon[0]))
            if not (JUMP_KM[0] <= d <= JUMP_KM[1]):
                continue
            br = bearing(pa[0], pa[1], B.lat[0], B.lon[0])
            if not (JUMP_BEARING[0] <= br <= JUMP_BEARING[1]):
                continue
            # hand-over: the first analysis at which the secondary is as deep as the primary,
            # else when the primary ends — either way within JUMP_TX_H of the secondary forming
            tx = None
            for i, tt in enumerate(B.t):
                if tt > B.t0 + JUMP_TX_H:
                    break
                q = A.at(tt)
                if q is None:
                    tx = tt
                    break
                if B.p[i] <= q[2]:
                    tx = tt
                    break
            if tx is None:
                continue
            qa = A.at(tx) or (A.lat[-1], A.lon[-1], A.p[-1])
            qb = B.at(tx)
            if float(hav(qa[0], qa[1], qb[0], qb[1])) > 1000:
                continue
            cands.append((d + 40 * (pB - pa[2]), A, B, tx))
    cands.sort(key=lambda c: c[0])
    for _score, A, B, tx in cands:
        if A.id in nxt or B.id in taken:
            continue
        nxt[A.id] = ("jump", B, tx)
        taken.add(B.id)

    # 3. chains → storms
    storms = []
    for A in tracks:
        if A.id in taken:
            continue
        parts, side = [], []
        cur, t_from, join = A, A.t0, None
        while True:
            nx = nxt.get(cur.id)
            if nx is None:
                parts.append((cur, t_from, cur.t1, join))
                break
            kind, B, tx = nx
            if kind == "cont":
                parts.append((cur, t_from, cur.t1, join))
                cur, t_from, join = B, B.t0, "cont"
                continue
            # jump at tx: the primary until its last analysis at or before tx, then the secondary
            ia = max(0, int(np.searchsorted(cur.t, tx, side="right")) - 1)
            ta = cur.t[ia]
            parts.append((cur, t_from, ta, join))
            if ta < cur.t1:
                side.append((cur, ta, cur.t1))
            if tx > B.t0:
                side.append((B, B.t0, tx))
            cur, t_from, join = B, tx, "jump"
        storms.append(Storm(parts=parts, side=side))
    return storms
