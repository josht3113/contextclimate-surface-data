import requests, os
H={"User-Agent":"contextclimate-surface-data (github.com/josht3113)"}
B='https://josht3113.github.io/contextclimate-hourly-data/'
os.makedirs('fx/hourly/hourly/KISP',exist_ok=True)
open('fx/hourly/index.json','wb').write(requests.get(B+'index.json',headers=H,timeout=60).content)
for D in (2000,2010,2020):
    r=requests.get(f'{B}hourly/KISP/{D}.bin.gz',headers=H,timeout=120); print(D,r.status_code,len(r.content))
    if r.ok: open(f'fx/hourly/hourly/KISP/{D}.bin.gz','wb').write(r.content)
u='https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=KISP&data=wxcodes&data=p01i&year1=2025&month1=12&day1=31&year2=2026&month2=11&day2=29&tz=UTC&format=onlycomma&latlon=no&missing=empty&trace=T&direct=no&report_type=3&report_type=4'
r=requests.get(u,headers=H,timeout=180); print('iem',r.status_code,len(r.content))
open('fx/iem_kisp_2026.csv','wb').write(r.content)
