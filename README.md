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

## Storm tracks (`storms/`)

Low-pressure tracks for the ContextClimate Storm Tracks map, rebuilt from every WPC
coded surface bulletin since 2003 whenever a new analysis is in.

- `storms/index.json` — years and storm counts, the latest analysis, CPC ONI / NAO / PNA.
- `storms/tracks/YYYY.json.gz` — storms that formed in that year (record format in `scripts/storms_build.py`).
- `storms/b/YYYY/MM.json.gz` — the month's bulletins (lows + fronts), every 3-hourly analysis.

Positions are smoothed (`scripts/storms_smooth.py`: a Gaussian-weighted local linear fit over
about ±12 h, wider for the whole-degree positions before 2011); pressures are WPC's.

Depth (`scripts/storms_depth.py`): mean sea-level pressure on a 750 km ring around each low minus
its central pressure — NCEP/NCAR Reanalysis 1 (6-hourly, 2.5°) through 17 Mar 2026, when it stopped
(precomputed: `static/storm-depth-r1.csv.gz`), GFS 1° analyses from the NOAA archive on AWS after
that (cached; 0.6 mb bias to R1 removed). A storm is *significant* when it is the deepest low within
800 km and at least 8 mb deep for 24 h in all, or a bomb, nor'easter or tropical cyclone; flag `u`
marks storms that passed within 300 km of the lower 48 (`scripts/us48_mask.json`).

Method (`scripts/storms_track.py`): each analysis's lows are matched to the next analysis by
minimum-cost assignment (distance from the predicted position plus pressure change, gated at
100 km/h); a low may miss one analysis. Storms last at least 24 h. A secondary low east of the
Rockies that forms 200–900 km downstream of a primary of 1012 mb or lower, deepens at least
4 mb in 24 h and becomes as deep as the primary while it dies takes over the storm (drawn as a
jump). Tropical cyclones are not in WPC's list of lows, so IBTrACS positions are added where
WPC has none within 300 km.

Bulletin months are kept in the Actions cache and fetched again from the IEM archive when missing.

## Sources

- NOAA High-Resolution Rapid Refresh (HRRR), NOAA Open Data Dissemination program (AWS); only the needed GRIB2 messages are downloaded (byte ranges from the `.idx` files).
- NWS Weather Prediction Center coded surface bulletins (AWIPS `CODSUS`), via the Iowa Environmental Mesonet text archive.
- IBTrACS v04r01 (NOAA NCEI), North Atlantic and eastern North Pacific — tropical cyclone positions for the storm tracks.
- NOAA CPC Oceanic Niño Index and monthly NAO / PNA indices.
- NCEP/NCAR Reanalysis 1 sea-level pressure (NOAA PSL) and GFS analyses (NOAA Open Data on AWS) — storm depth.

## Workflow

**Update surface data** runs hourly (and on demand): restores the window from the Actions cache, downloads only new hours, deploys Pages.

## One-time setup

Settings → Pages → Build and deployment → Source: **GitHub Actions**.
