"""Track smoothing.

WPC places each low by hand every 3 hours, so a raw track kinks at nearly every analysis
(median turn ~67° per step). Each jump-free stretch is smoothed with a Gaussian-weighted
local linear fit in time on a local flat plane; the fit is linear, so stretch ends are not
pulled inward. Whole-degree positions (before 2011) get a wider window.
"""
from __future__ import annotations

import math

import numpy as np

SIGMA_H = 5.0          # 0.1° positions
SIGMA_H_COARSE = 7.5   # whole-degree positions
HALF_WIN = 2.5         # window half-width in sigmas


def _fit(t, v, w):
    """Weighted linear fit value at t=0 (t already centred)."""
    sw = w.sum()
    if sw <= 0:
        return float(v[np.argmin(np.abs(t))])
    mt = (w * t).sum() / sw
    mv = (w * v).sum() / sw
    var = (w * (t - mt) ** 2).sum()
    if var < 1e-9:
        return mv
    slope = (w * (t - mt) * (v - mv)).sum() / var
    return mv + slope * (0 - mt)


def smooth_stretch(t, lat, lon, sigma):
    """t (hours), lat, lon arrays for one jump-free stretch → smoothed lat, lon."""
    n = len(t)
    if n < 3:
        return np.array(lat, float), np.array(lon, float)
    t = np.asarray(t, float); lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    lat0 = float(lat.mean())
    k = math.cos(math.radians(lat0))
    x = lon * k * 111.2
    y = lat * 111.2
    xs, ys = np.empty(n), np.empty(n)
    for i in range(n):
        dt = t - t[i]
        m = np.abs(dt) <= HALF_WIN * sigma
        w = np.exp(-0.5 * (dt[m] / sigma) ** 2)
        xs[i] = _fit(dt[m], x[m], w)
        ys[i] = _fit(dt[m], y[m], w)
    return ys / 111.2, xs / (k * 111.2)


def smooth_points(pts, coarse: bool):
    """pts: [[t, lat, lon, p, jump, src], ...] → same list with lat/lon smoothed in place
    (positions only; pressures are WPC's). Stretches break at jumps."""
    sigma = SIGMA_H_COARSE if coarse else SIGMA_H
    starts = [0] + [i for i, p in enumerate(pts) if p[4]] + [len(pts)]
    for a, b in zip(starts[:-1], starts[1:]):
        if b - a < 3:
            continue
        seg = pts[a:b]
        la, lo = smooth_stretch([p[0] for p in seg], [p[1] for p in seg], [p[2] for p in seg], sigma)
        for j, p in enumerate(seg):
            p[1], p[2] = float(la[j]), float(lo[j])
    return pts


def turning(lat, lon):
    """Mean absolute change of heading between steps (degrees)."""
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    x = lon * np.cos(np.radians(lat)); v = np.diff(np.c_[x, lat], axis=0)
    a = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
    d = np.abs((np.diff(a) + 180) % 360 - 180)
    return float(d.mean()) if len(d) else float("nan")
