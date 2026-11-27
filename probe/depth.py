"""Dev: ring depth (field ring-mean minus WPC central pressure) for every WPC low, from NCEP R1
(2003 – Mar 2026) and GFS 1° analyses on AWS (overlap Jan–Mar 2026 + after)."""
import requests, os, io, time, gzip, json, glob, math, sys
import numpy as np, netCDF4
H={"User-Agent":"contextclimate-surface-data (github.com/josht3113)"}
OUT='depth_out'; os.makedirs(OUT, exist_ok=True)
log=open(f'{OUT}/log.txt','w')
def P(*a):
    s=' '.join(map(str,a)); print(s, flush=True); log.write(s+'\n'); log.flush()
# --- probe R2 and GFS
for u in ['https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/gaussian_grid/mslp.2026.nc',
          'https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/surface/mslp.2026.nc']:
    try:
        r=requests.head(u,headers=H,timeout=60,allow_redirects=True); P('R2', u, r.status_code, r.headers.get('content-length'), r.headers.get('last-modified'))
    except Exception as e: P('R2', u, e)
GB='https://noaa-gfs-bdp-pds.s3.amazonaws.com'
def gfs_prmsl(ymd, hh):
    base=f'{GB}/gfs.{ymd}/{hh}/atmos/gfs.t{hh}z.pgrb2.1p00.f000'
    r=requests.get(base+'.idx',headers=H,timeout=60)
    if not r.ok: return None, f'idx {r.status_code}'
    lines=r.text.strip().split('\n')
    for i,l in enumerate(lines):
        if ':PRMSL:mean sea level:' in l:
            a=int(l.split(':')[1]); b=int(lines[i+1].split(':')[1])-1 if i+1<len(lines) else ''
            g=requests.get(base,headers={**H,'Range':f'bytes={a}-{b}'},timeout=60)
            import pygrib
            f=open('/tmp/x.grb2','wb'); f.write(g.content); f.close()
            m=pygrib.open('/tmp/x.grb2').read(1)[0]
            return m.values/100.0, len(g.content)   # mb, (181,360) lat 90..-90, lon 0..359
    return None, 'no PRMSL'
for ymd,hh in [('20210301','00'),('20260115','12'),('20261126','18')]:
    t=time.time(); v,n=gfs_prmsl(ymd,hh); P('GFS',ymd,hh, None if v is None else v.shape, n, round(time.time()-t,1),'s', None if v is None else (round(float(v.min()),1), round(float(v.max()),1)))
# --- bulletins
frames=[]
for f in sorted(glob.glob('dev_out/bulletins/*/*.json.gz')):
    for b in json.load(gzip.open(f,'rt'))['b']:
        if b['t']>='2003010100': frames.append((b['t'], b['lows']))
P('bulletins', len(frames))
AZ=np.radians(np.arange(0,360,15)); RADII=(500,750,1000)
def ring_pts(lat, lon, r):
    d=r/6371.0; la=math.radians(lat); lo=math.radians(lon)
    la2=np.arcsin(np.sin(la)*np.cos(d)+np.cos(la)*np.sin(d)*np.cos(AZ))
    lo2=lo+np.arctan2(np.sin(AZ)*np.sin(d)*np.cos(la), np.cos(d)-np.sin(la)*np.sin(la2))
    return np.degrees(la2), (np.degrees(lo2)+360)%360
def sample(field, lat0, dlat, dlon, la, lo):   # field[lat_i, lon_j], lat from lat0 descending by dlat
    y=(lat0-la)/dlat; x=lo/dlon; ny,nx=field.shape
    y0=np.clip(np.floor(y).astype(int),0,ny-2); x0=np.floor(x).astype(int)%nx; x1=(x0+1)%nx
    fy=y-y0; fx=x-np.floor(x)
    return (field[y0,x0]*(1-fx)*(1-fy)+field[y0,x1]*fx*(1-fy)+field[y0+1,x0]*(1-fx)*fy+field[y0+1,x1]*fx*fy)
def depths(field, lat0, dlat, dlon, lows):
    rows=[]
    for la,lo,p in lows:
        r=[]
        for R in RADII:
            a,b=ring_pts(la,lo,R); r.append(round(float(np.mean(sample(field,lat0,dlat,dlon,a,b)))-p,1))
        c=round(float(sample(field,lat0,dlat,dlon,np.array([la]),np.array([(lo+360)%360]))[0]),1)
        rows.append([la,lo,p]+r+[c])
    return rows
from datetime import datetime, timedelta
def th(k): return datetime(int(k[:4]),int(k[4:6]),int(k[6:8]),int(k[8:10]))
# --- R1 years
out=gzip.open(f'{OUT}/depth_r1.csv.gz','wt'); out.write('t,lat,lon,p,d500,d750,d1000,fc\n')
cur=None; F=None; T0=None
for k,lows in frames:
    d=th(k); y=d.year
    if y!=cur:
        u=f'https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/slp.{y}.nc'
        r=requests.get(u,headers=H,timeout=600)
        if not r.ok: P('R1',y,r.status_code); cur=y; F=None; continue
        fn=f'/tmp/r1_{y}.nc'; open(fn,'wb').write(r.content); ds=netCDF4.Dataset(fn)
        F=np.array(ds.variables['slp'][:],dtype=np.float32)/100.0; tt=ds.variables['time']
        T0=netCDF4.num2date(tt[0],tt.units); nT=F.shape[0]; cur=y; P('R1',y,F.shape, str(netCDF4.num2date(tt[-1],tt.units)))
        ds.close(); os.remove(fn)
        nextF=None
    if F is None: continue
    hrs=(d-datetime(T0.year,T0.month,T0.day,T0.hour)).total_seconds()/3600; i=hrs/6
    if i>F.shape[0]-1: continue
    i0=int(math.floor(i)); f=i-i0
    fld=F[i0] if f==0 else (F[i0]*(1-f)+F[min(i0+1,F.shape[0]-1)]*f)
    for row in depths(fld, 90.0, 2.5, 2.5, lows): out.write(k+','+','.join(map(str,row))+'\n')
out.close()
# --- GFS overlap Jan 1 – Mar 17 2026 (00/06/12/18 only) + a later sample month
out=gzip.open(f'{OUT}/depth_gfs.csv.gz','wt'); out.write('t,lat,lon,p,d500,d750,d1000,fc\n')
n=0
for k,lows in frames:
    if not (('2026010100'<=k<'2026031800') or ('2026100100'<=k)): continue
    if k[8:10] not in ('00','06','12','18'): continue
    try: v,_=gfs_prmsl(k[:8],k[8:10])
    except Exception as e: P('gfs err',k,e); continue
    if v is None: continue
    # coarsen to 2.5 deg like R1 (box mean) for a fair comparison? keep 1 deg; ring means are smooth anyway
    for row in depths(v.astype(np.float32), 90.0, 1.0, 1.0, lows): out.write(k+','+','.join(map(str,row))+'\n')
    n+=1
out.close(); P('gfs analyses', n)
