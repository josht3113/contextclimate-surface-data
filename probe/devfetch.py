import sys, os, time, json, csv, io, gzip
sys.path.insert(0, 'scripts')
import requests
from storms_csb import sync, UA
from datetime import datetime
root = 'dev_out/bulletins'
t = time.time()
n = sync(root, first_year=2002, log=print)
print('months fetched', n, round(time.time() - t), 's')
os.makedirs('dev_out/side', exist_ok=True)
for name, u in {"oni": "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
                "nao": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii",
                "pna": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii"}.items():
    open(f'dev_out/side/{name}.txt', 'w').write(requests.get(u, headers=UA, timeout=60).text)
base = "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/"
rows = []
for f in ("ibtracs.NA.list.v04r01.csv", "ibtracs.EP.list.v04r01.csv", "ibtracs.last3years.list.v04r01.csv"):
    txt = requests.get(base + f, headers=UA, timeout=300).text
    rd = csv.reader(io.StringIO(txt)); hdr = next(rd); next(rd)
    ix = {h: i for i, h in enumerate(hdr)}
    for r in rd:
        if r[ix['SEASON']] < '2002': continue
        rows.append([r[ix['SID']], r[ix['SEASON']], r[ix['BASIN']], r[ix['NAME']], r[ix['ISO_TIME']], r[ix['NATURE']], r[ix['LAT']], r[ix['LON']], r[ix['USA_WIND']], r[ix['USA_PRES']], r[ix['TRACK_TYPE']], r[ix['USA_SSHS']] if 'USA_SSHS' in ix else ''])
with gzip.open('dev_out/side/ibtracs.csv.gz', 'wt') as f:
    w = csv.writer(f); w.writerow(['sid','season','basin','name','time','nature','lat','lon','wind','pres','track','sshs']); w.writerows(rows)
print('ibtracs rows', len(rows))
