import requests, json, time, os, re, sys
OUT = 'probe_out'; os.makedirs(OUT, exist_ok=True)
H = {"User-Agent": "contextclimate-surface-data (github.com/josht3113)"}
IEM = "https://mesonet.agron.iastate.edu/cgi-bin/afos/retrieve.py"
log = {}
def get(url, params=None, name=None, keep=200000):
    t = time.time()
    try:
        r = requests.get(url, params=params, headers=H, timeout=180)
        info = {"status": r.status_code, "bytes": len(r.content), "sec": round(time.time()-t, 1), "url": r.url,
                "ctype": r.headers.get("content-type"), "acao": r.headers.get("access-control-allow-origin")}
        if name:
            open(os.path.join(OUT, name), 'wb').write(r.content[:keep])
        return r, info
    except Exception as e:
        return None, {"error": str(e)}
# coverage: one day per year (Jan 15) 2003..2026
cov = {}
for y in range(2002, 2027):
    r, info = get(IEM, {"pil": "CODSUS", "sdate": f"{y}-01-15", "edate": f"{y}-01-16", "fmt": "text", "limit": 9999}, f"day_{y}.txt", 400000)
    n = r.text.count("CODED SURFACE") if r is not None else -1
    nv = len(re.findall(r"VALID\s+\d{6}Z", r.text)) if r is not None else -1
    hires = len(re.findall(r"\b\d{7}\b", r.text)) if r is not None else -1
    info.update(n=n, nvalid=nv, hires_tokens=hires); cov[y] = info
log["coverage"] = cov
# July sample too
for y in (2005, 2015, 2025):
    r, info = get(IEM, {"pil": "CODSUS", "sdate": f"{y}-07-15", "edate": f"{y}-07-16", "fmt": "text", "limit": 9999}, f"jul_{y}.txt", 400000)
    log[f"jul_{y}"] = info
# full month timing
r, info = get(IEM, {"pil": "CODSUS", "sdate": "2010-01-01", "edate": "2010-02-01", "fmt": "text", "limit": 9999}, "month_2010_01.txt", 20000000)
info["n"] = r.text.count("CODED SURFACE") if r is not None else -1
log["month_2010_01"] = info
r, info = get(IEM, {"pil": "CODSUS", "sdate": "2025-01-01", "edate": "2025-02-01", "fmt": "text", "limit": 9999}, "month_2025_01.txt", 20000000)
info["n"] = r.text.count("CODED SURFACE") if r is not None else -1
log["month_2025_01"] = info
# other possible PILs
for pil in ("CODSUSHR", "CODSUS2", "CODSAS", "CODSFP"):
    r, info = get(IEM, {"pil": pil, "fmt": "text", "limit": 3}, f"pil_{pil}.txt", 50000)
    log["pil_"+pil] = info
# metadata search for list of products
r, info = get("https://mesonet.agron.iastate.edu/api/1/nws/afos/list.json", {"pil": "CODSUS", "date": "2026-11-20"}, "afos_list.json", 100000)
log["afos_list"] = info
# zenodo
for rec in ("2642801", "2646544"):
    r, info = get(f"https://zenodo.org/api/records/{rec}", None, f"zenodo_{rec}.json", 2000000)
    log["zenodo_"+rec] = info
# WPC live
for u in ("https://www.wpc.ncep.noaa.gov/discussions/codsus", "https://www.wpc.ncep.noaa.gov/discussions/codsus_hr",
          "https://tgftp.nws.noaa.gov/data/raw/as/asus02.kwbc.cod.sus.txt", "https://tgftp.nws.noaa.gov/data/raw/as/asus01.kwbc.cod.sus.txt"):
    r, info = get(u, None, "wpc_" + re.sub(r'\W', '_', u[-30:]) + ".txt", 50000)
    log["wpc_"+u] = info
# side data
for name, u in {"oni": "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
                "nao": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii",
                "pna": "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii"}.items():
    r, info = get(u, None, name + ".txt", 300000); log[name] = info
base = "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/"
r, info = get(base + "ibtracs.NA.list.v04r01.csv", None, "ibtracs_na_head.csv", 30000); log["ibtracs_na"] = info
r, info = get(base + "ibtracs.last3years.list.v04r01.csv", None, "ibtracs_l3_head.csv", 30000); log["ibtracs_l3"] = info
# station wx codes
r, info = get("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py",
              {"station": "ISP", "data": ["wxcodes", "p01i", "tmpf"], "year1": 2025, "month1": 2, "day1": 1, "year2": 2025, "month2": 3, "day2": 1,
               "tz": "Etc/UTC", "format": "onlycomma", "latlon": "no", "missing": "M", "trace": "T", "report_type": [3, 4]}, "isp_wx_2025_02.csv", 2000000)
log["isp_wx"] = info
r, info = get("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py",
              {"station": "ISP", "data": ["wxcodes", "p01i"], "year1": 2003, "month1": 1, "day1": 1, "year2": 2004, "month2": 1, "day2": 1,
               "tz": "Etc/UTC", "format": "onlycomma", "latlon": "no", "missing": "M", "trace": "T", "report_type": [3, 4]}, "isp_wx_2003.csv", 300000)
log["isp_wx_year"] = info
# hourly archive repo
r, info = get("https://josht3113.github.io/contextclimate-hourly-data/", None, "hourly_index.html", 100000); log["hourly_idx"] = info
r, info = get("https://josht3113.github.io/contextclimate-hourly-data/manifest.json", None, "hourly_manifest.json", 500000); log["hourly_manifest"] = info
json.dump(log, open(os.path.join(OUT, "log.json"), "w"), indent=1)
print(json.dumps(log, indent=1)[:5000])
