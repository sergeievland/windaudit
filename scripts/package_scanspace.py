"""Finalise validation logs, manifest and a checked source/data archive."""
from pathlib import Path
import json,re,subprocess,sys,zipfile
ROOT=Path(__file__).resolve().parents[1]
def main(output):
    output=str(Path(output).resolve())
    bad=[]
    for p in ROOT.rglob('*'):
        if not p.is_file() or any(x in {'out','__pycache__','.pytest_cache'} for x in p.parts):continue
        if p.suffix in {'.py','.md','.txt','.json','.patch','.sh','.toml','.yml'}:
            if re.search('[\u0400-\u04ff]',p.read_text()):bad.append(str(p.relative_to(ROOT)))
    assert not bad,bad
    pack=[sys.executable,'scripts/package_release.py',output]
    subprocess.run(pack,cwd=ROOT,check=True)
    result=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    assert result.returncode==0,result.stdout
    with (ROOT/'results/tests.txt').open('a') as f:
        f.write('\nPACKAGE VALIDATION\n$ sha256sum -c SHA256SUMS\n'+result.stdout)
        f.write('All manifest entries verified; manifest rebuilt once more after recording this log.\nCyrillic text scan: clean.\n')
    subprocess.run(pack,cwd=ROOT,check=True)
    result=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    assert result.returncode==0,result.stdout
    with zipfile.ZipFile(output) as z:
        assert not any('__pycache__' in n or n.endswith('.pyc') for n in z.namelist())
        manifest=z.read('windaudit-0.3.0/SHA256SUMS').decode()
        import hashlib
        for line in manifest.splitlines():
            digest,name=line.split('  ',1)
            assert hashlib.sha256(z.read('windaudit-0.3.0/'+name)).hexdigest()==digest,name
    print(json.dumps({'archive':output,'sha256sum_check':'passed','archive_bytes_check':'passed','pycache_entries':0,'cyrillic_text_files':0}))
if __name__=='__main__':main(sys.argv[1])
