"""Metadata-only selection with conservative directory-size bounds. No TIFF GETs."""
import concurrent.futures, hashlib, json, re, sys, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/'
STATE=ROOT/'out/wide_discovery'
STATE.mkdir(parents=True,exist_ok=True)
LOW,HIGH=6000,18000
CAP=1_000_000_000

def get(url,head=False):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,method='HEAD' if head else 'GET'),timeout=45) as response:
                return int(response.headers['Content-Length']) if head else response.read()
        except Exception:
            if attempt==2:raise

def metadata(name):
    path=STATE/'meta'/ (name+'.json');path.parent.mkdir(exist_ok=True)
    try:
        if not path.exists():path.write_bytes(get(BASE+name+'/meta.json'))
        d=json.loads(path.read_text())
        b=d.get('bbox')
        hit=bool(b and any(b[0][k]<HIGH and b[1][k]>=LOW for k in [0,2]))
        return {'id':name,'metadata':d,'intersects':hit}
    except Exception as error:return {'id':name,'error':str(error)}

def sizes(row):
    name=row['id'];p=STATE/'sizes'/(name+'.json');p.parent.mkdir(exist_ok=True)
    try:
        if p.exists():return json.loads(p.read_text())
        listing=get(BASE+name+'/').decode()
        names=re.findall('href="([^"]+)"',listing)
        wanted=['meta.json','x.tif','y.tif','z.tif']+(['mask.tif'] if 'mask.tif' in names else [])
        assert all(f in names for f in wanted)
        sizes_html=dict(re.findall(r'href="([^"]+)"[^>]*>.*?</a></td><td class="size">([^<]+)',listing))
        files={}
        for f in wanted:
            value,unit=sizes_html[f].split()
            multiplier={'B':1,'KiB':1024,'MiB':1024**2,'GiB':1024**3}[unit]
            # The listing rounds to one decimal place; add a full last-place
            # unit so this remains an upper bound under either rounding rule.
            from decimal import Decimal,ROUND_CEILING
            increment=Decimal(10)**(-len(value.split('.')[1])) if '.' in value else Decimal(1)
            files[f]=int(((Decimal(value)+increment)*multiplier).to_integral_value(rounding=ROUND_CEILING))
        r={'id':name,'files':files,'bytes':sum(files.values()),'bbox':row['metadata'].get('bbox'),'size_basis':'conservative upper bound from rounded directory listing'}
        p.write_text(json.dumps(r,indent=2));return r
    except Exception as error:return {'id':name,'error':str(error)}

def main(previous):
    previous=Path(previous)
    cached={}
    for filename in ['metadata_selection.json','target_metadata.json','stratified_metadata.json','samewrap_metadata.json']:
        for row in json.loads((previous/filename).read_text()):
            if 'metadata' in row:cached[row['id']]=row['metadata']
    (STATE/'meta').mkdir(exist_ok=True)
    for name,d in cached.items():(STATE/'meta'/(name+'.json')).write_text(json.dumps(d))
    index=(previous/'patch_index.html').read_text()
    auto=[n.rstrip('/') for n in re.findall('href="([^"]+)"',index) if n.startswith('auto_')]
    # Outcome-blind expansion: 1,024 evenly spaced entries of the auto-grown family.
    sampled=[auto[round(i*(len(auto)-1)/1023)] for i in range(1024)]
    names=sorted(set(cached)|set(sampled))
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(64) as pool:
        futures={pool.submit(metadata,n):n for n in names}
        for future in concurrent.futures.as_completed(futures):
            rows.append(future.result())
            if len(rows)%100==0:print('metadata',len(rows),'/',len(names),flush=True)
    rows.sort(key=lambda r:r['id'])
    (STATE/'metadata_census.json').write_text(json.dumps(rows,indent=2))
    hits=[r for r in rows if r.get('intersects')]
    measured=[]
    print('bbox candidates',len(hits),flush=True)
    with concurrent.futures.ThreadPoolExecutor(64) as pool:
        for future in concurrent.futures.as_completed([pool.submit(sizes,r) for r in hits]):
            measured.append(future.result())
            if len(measured)%50==0:print('sized',len(measured),'/',len(hits),flush=True)
    # Preserve the cached named cohort first; then auto-grown IDs, without
    # inspecting any winding inconsistency or changing the z band after outcomes.
    measured.sort(key=lambda r:(r['id'] not in cached,r['id']))
    chosen=[];over_cap=[];total=0
    for row in measured:
        if 'error' in row:continue
        if total+row['bytes']>CAP:over_cap.append(row['id']);continue
        chosen.append(row);total+=row['bytes']
    result={'schema_version':1,'band':[LOW,HIGH],'budget_bytes':CAP,'total_bytes':total,'size_basis':'conservative pre-download upper bound; actual byte counts recorded after download',
            'patch_count':len(chosen),'patches':chosen,'skipped_over_budget':over_cap,
            'metadata_count':len(rows),'metadata_errors':[r for r in rows if 'error' in r],
            'size_errors':[r for r in measured if 'error' in r],
            'selection':'All 689 cached metadata records plus 1024 evenly spaced auto-grown directory entries; conservative bbox intersection, then deterministic 1 GB cap. Not an exhaustive catalogue scan.',
            'catalogue_sha256':hashlib.sha256((previous/'patch_index.html').read_bytes()).hexdigest()}
    destination=ROOT/'results/wide';destination.mkdir(parents=True,exist_ok=True)
    (destination/'download_plan.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PLAN',len(chosen),total,flush=True)
if __name__=='__main__':main(sys.argv[1])
