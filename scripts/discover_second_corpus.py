"""Inspect public directory indexes only; never download volumes or tracks."""
import concurrent.futures
import hashlib
import json
import re
import urllib.request
from pathlib import Path
BASE='https://dl.ash2txt.org/datasets/spiral_datasets/'
NEEDED=['relative_windings.json','same_windings.json','abs_winding.json','umbilicus.json']
def listing(url):
    raw=urllib.request.urlopen(url,timeout=40).read()
    return dict(url=url,sha256=hashlib.sha256(raw).hexdigest(),links=re.findall(r'href="([^"]+)"',raw.decode()))
def inspect(name):
    root=listing(BASE+name);rows=[root]
    for d in root['links']:
        if d.endswith('/') and d!='../' and not d.startswith('/'):
            scan=listing(root['url']+d);rows.append(scan)
            for dd in scan['links']:
                if dd=='tracks/':rows.append(listing(scan['url']+dd))
    return dict(dataset=name.rstrip('/'),indexes=rows,complete_input_locations=[x['url'] for x in rows if set(NEEDED)<=set(x['links'])])
def main():
    root=listing(BASE)
    names=[x for x in root['links'] if x.startswith('PHerc') and x.endswith('/') and x!='PHercParis4/']
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:rows=list(ex.map(inspect,names))
    report=dict(checked_utc='2026-09-18',required_files=NEEDED,root_index=root,datasets=rows,
                status='available' if any(r['complete_input_locations'] for r in rows) else 'blocked_no_second_complete_corpus_in_listed_locations',
                scope='Root, scan and tracks indexes for every listed non-Paris dataset. No volumes, tracks or unlisted guessed paths downloaded.')
    p=Path('results/second_corpus/discovery.json');p.parent.mkdir(exist_ok=True);p.write_text(json.dumps(report,indent=1)+'\n')
    print(report['status'],len(rows),flush=True)
if __name__=='__main__':main()
