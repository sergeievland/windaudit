"""Retry only unresolved diagnostics after valid cycle-cut strengthening."""
import concurrent.futures
import json
import multiprocessing
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from diagnostic_order import initialise,trial

def main():
    result=json.loads(Path('results/paris4/order_diagnostic.json').read_text())
    old=[x for x in result['trials'] if x['membership']=='unknown']
    jobs=[]
    # Reconstruct the original site order from the pinned collections.
    from windaudit.pcl import load_point_collections
    from windaudit.controls import _ordered_pids
    for t in old:
        kind,cid=t['collection'].split(':');cid=int(cid)
        cols=load_point_collections('data/paris4/'+('relative_windings.json' if kind=='relative' else 'same_windings.json'))
        jobs.append((t['cohort'],t['defect'],(cid,t['boundary'],_ordered_pids(cols[cid],kind)),t['direction']))
    with concurrent.futures.ProcessPoolExecutor(max_workers=3,initializer=initialise,mp_context=multiprocessing.get_context('spawn')) as pool:
        new=list(pool.map(trial,jobs))
    report=dict(reason='Retry unresolved cases after adding mathematically valid three-node cycle cuts; same data, margins, objectives and 60-second budget.',original=old,retry=new)
    Path('results/paris4/order_retries.json').write_text(json.dumps(report,indent=1)+'\n')
    print('Retried',len(old),'resolved',sum(t['membership']!='unknown' for t in new))
if __name__=='__main__':main()
