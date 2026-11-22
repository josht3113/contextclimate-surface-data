# contextclimate-surface-data

Hourly surface weather grids and NWS frontal analyses for the ContextClimate
Surface Analysis Builder, published to GitHub Pages.

## What's published

- `manifest.json` — grid, variables, the analysis hours available, the forecast run, a check against METARs.
- `a/YYYYMMDDHH.bin.gz` — HRRR analysis (forecast hour 0) for each of the past 24 hours (UTC).
- `f/YYYYMMDDHH/FF.bin.gz` — the newest complete HRRR run, forecast hours 1–18.
- `fronts.json` — NWS Weather Prediction Center coded surface bulletins (fronts, highs, lows) in the same window.
- `terrain.bin.gz` — HRRR model surface height (m), one int16 plane, written once.

Frame files: int16 little-endian planes `t, td, u, v, p, refc` on a 0.125° grid
(126°W–65°W, 22°N–51°N; rows south→north, columns west→east). Bytes shuffled
(all low bytes, then all high bytes), then gzip. value = stored × scale + offset
(see `manifest.json`); −32768 = outside the HRRR domain.

- `t`, `td`: 2 m temperature and dewpoint, °F
- `u`, `v`: 10 m wind, m/s, earth-relative (rotated from the Lambert grid)
- `p`: MAPS mean sea-level pressure (MSLMA), mb
- `refc`: composite reflectivity, dBZ

Regridding: mean of the 3 km HRRR points within ~0.09° of each grid point.

## Sources

- NOAA High-Resolution Rapid Refresh (HRRR), NOAA Open Data Dissemination program (AWS); only the needed GRIB2 messages are downloaded (byte ranges from the `.idx` files).
- NWS Weather Prediction Center coded surface bulletins (AWIPS `CODSUS`), via the Iowa Environmental Mesonet text archive.

## Workflow

**Update surface data** runs hourly (and on demand): restores the window from the Actions cache, downloads only new hours, deploys Pages.

## One-time setup

Settings → Pages → Build and deployment → Source: **GitHub Actions**.
