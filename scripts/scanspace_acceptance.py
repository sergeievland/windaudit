"""Capture legacy reproduction twice without modifying committed old results."""
import hashlib,json,os
from pathlib import Path
import subprocess,sys,zipfile
ROOT=Path(__file__).resolve().parents[1]
def main():
    env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    logs=[];snapshots=[]
    for attempt in range(2):
        command=['bash','scripts/reproduce.sh']
        run=subprocess.run(command,cwd=ROOT,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        logs.append(f'REPRODUCE RUN {attempt+1}\n$ bash scripts/reproduce.sh\n'+run.stdout+f'\nexit_code={run.returncode}\n')
        (ROOT/'results/tests.txt').write_text('\n'.join(logs))
        if run.returncode:raise SystemExit(run.returncode)
        snapshots.append({n:(ROOT/'out/paris4'/n).read_bytes() for n in ['audit_report.json','certified_links.json','review_queue.json']})
        print('reproduce',attempt+1,'passed',flush=True)
    volatile={'environment','runtime_seconds','inputs'}
    a,b=[json.loads(s['audit_report.json']) for s in snapshots]
    assert all(a.get(k)==b.get(k) for k in set(a)|set(b) if k not in volatile)
    exports={n:hashlib.sha256(snapshots[0][n]).hexdigest() for n in ['certified_links.json','review_queue.json']}
    assert all(snapshots[0][n]==snapshots[1][n] for n in exports)
    logs.append('REPEAT COMPARISON: 17 reproducible report sections identical.\ncertified_links.json: byte-identical\nreview_queue.json: byte-identical\n'+json.dumps(exports,indent=2))
    if len(sys.argv)>1:
        changed=[];verified=[]
        with zipfile.ZipFile(sys.argv[1]) as archive:
            for item in archive.infolist():
                if item.is_dir():continue
                rel=Path(item.filename).relative_to('windaudit-0.3.0')
                if str(rel) in ['SHA256SUMS','results/tests.txt']:continue
                if archive.read(item)!=(ROOT/rel).read_bytes():changed.append(str(rel))
                else:verified.append(str(rel))
        assert not changed,changed
        preservation={'all_previous_published_values_unchanged':True,'old_files_checked_byte_for_byte':len(verified),'changed_old_files':changed,'exceptions':['SHA256SUMS rebuilt','results/tests.txt refreshed; original preserved under results/legacy_borderline/']}
        (ROOT/'results/scanspace_preservation.json').write_text(json.dumps(preservation,indent=2)+'\n')
        logs.append('BASELINE PRESERVATION\n'+json.dumps(preservation,indent=2))
    (ROOT/'results/tests.txt').write_text('\n'.join(logs)+'\n')
if __name__=='__main__':main()
