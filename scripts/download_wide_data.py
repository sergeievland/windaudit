"""Download only the preannounced plan, checking size bounds and hashing files."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import shutil
import ssl
import sys
import threading
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
BASE='https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/'
try:
    import requests
except ImportError:
    requests=None
LOCAL=threading.local()

def transfer(url,temp):
    if requests is None:
        with urllib.request.urlopen(url,timeout=90) as response,temp.open('wb') as out:
            expected=int(response.headers['Content-Length'])
            shutil.copyfileobj(response,out,1024*1024)
    else:
        if not hasattr(LOCAL,'session'):
            LOCAL.session=requests.Session();LOCAL.session.headers['Accept-Encoding']='identity'
            LOCAL.session.verify=ssl.get_default_verify_paths().cafile or True
        with LOCAL.session.get(url,stream=True,timeout=90) as response,temp.open('wb') as out:
            response.raise_for_status();expected=int(response.headers['Content-Length'])
            for chunk in response.iter_content(1024*1024):out.write(chunk)
    assert temp.stat().st_size==expected,('incomplete response',url)

def main(cache_root):
    plan=json.loads((ROOT/'results/wide/download_plan.json').read_text())
    assert plan['total_bytes']<=1_000_000_000
    assert sum(sum(p['files'].values()) for p in plan['patches'])==plan['total_bytes']
    dest=ROOT/'data/scanspace_wide';dest.mkdir(exist_ok=True)
    cache_root=Path(cache_root)
    caches=[cache_root/p for p in ['patches','selected_patches','supplementary_patches']]
    prior_ids=set()
    for filename in ['metadata_selection.json','target_metadata.json','stratified_metadata.json','samewrap_metadata.json']:
        if (cache_root/filename).exists():
            prior_ids.update(r['id'] for r in json.loads((cache_root/filename).read_text()) if 'metadata' in r)
    def fetch(row):
        name=row['id'];folder=dest/name;folder.mkdir(exist_ok=True);files={}
        for filename,size in row['files'].items():
            path=folder/filename
            if not path.exists():
                cached=next((c/name/filename for c in caches if (c/name/filename).exists() and (c/name/filename).stat().st_size<=size),None)
                raw_meta=ROOT/'out/wide_discovery/meta'/(name+'.json')
                if filename=='meta.json' and name not in prior_ids and raw_meta.exists():
                    # New metadata was cached verbatim by plan_wide_data.py.
                    # Prior records were reserialised and must not be used here.
                    cached=raw_meta
                if cached is not None:shutil.copyfile(cached,path)
                else:
                    temp=path.with_suffix(path.suffix+'.partial')
                    for attempt in range(3):
                        try:
                            transfer(BASE+name+'/'+filename,temp)
                            assert 0<temp.stat().st_size<=size,(name,filename,temp.stat().st_size,size)
                            temp.replace(path);break
                        except Exception:
                            if attempt==2:raise
            assert 0<path.stat().st_size<=size
            files[filename]={'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        return {'id':name,'files':files}
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(96) as pool:
        for future in concurrent.futures.as_completed([pool.submit(fetch,row) for row in plan['patches']]):
            rows.append(future.result())
            if len(rows)%25==0:print('downloaded',len(rows),'/',plan['patch_count'],flush=True)
    rows.sort(key=lambda r:r['id'])
    total=sum(f['bytes'] for r in rows for f in r['files'].values())
    (ROOT/'results/wide/dataset_manifest.json').write_text(json.dumps({'base_url':BASE,'patches':rows,'total_bytes':total,'announced_upper_bound':plan['total_bytes']},indent=2)+'\n')
    print('COMPLETE',len(rows),total,flush=True)

if __name__=='__main__':main(sys.argv[1])
