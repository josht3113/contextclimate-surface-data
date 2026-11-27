import requests, os, time
H={"User-Agent":"contextclimate-surface-data (github.com/josht3113)"}
os.makedirs('fx',exist_ok=True)
for u in ['https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/slp.2025.nc',
          'https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/slp.2026.nc',
          'https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis/surface/slp.2026.nc',
          'https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface/slp.mon.ltm.1991-2020.nc']:
    t=time.time()
    try:
        r=requests.get(u,headers=H,timeout=300); print(u, r.status_code, len(r.content), r.headers.get('last-modified'), round(time.time()-t,1),'s', flush=True)
        if r.ok and u.endswith('surface/slp.2026.nc') and 'downloads' in u: open('fx/slp.2026.nc','wb').write(r.content)
    except Exception as e: print(u, e)
import netCDF4, numpy as np
d=netCDF4.Dataset('fx/slp.2026.nc'); print(d.file_format, list(d.variables))
v=d.variables['slp']; print(v.shape, v.dtype, getattr(v,'units',None), getattr(v,'scale_factor',None))
t=d.variables['time']; print(t.units, netCDF4.num2date(t[-1],t.units))
print(d.variables['lat'][:3], d.variables['lon'][:3])
