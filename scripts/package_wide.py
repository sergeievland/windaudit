"""Build an additive release and assert preservation of every baseline file."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]

def included(path):
    rel=path.relative_to(ROOT)
    return (path.is_file() and not any(p in {'out','__pycache__','.pytest_cache','.git'} or p.endswith('.egg-info') for p in rel.parts)
            and path.suffix not in {'.pyc','.partial'} and not rel.parts[0].startswith('out_'))

def main(baseline,destination):
    baseline=Path(baseline);destination=Path(destination)
    with zipfile.ZipFile(baseline) as old:
        names=[n for n in old.namelist() if not n.endswith('/')]
        for name in names:
            rel=Path(name).relative_to(ROOT.name)
            assert (ROOT/rel).read_bytes()==old.read(name),f'Baseline changed: {rel}'
    files=sorted(p for p in ROOT.rglob('*') if included(p) and p.name!='SHA256SUMS.wide')
    manifest=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix()+'\n' for p in files)
    (ROOT/'SHA256SUMS.wide').write_text(manifest)
    files.append(ROOT/'SHA256SUMS.wide')
    for path in files:
        try:content=path.read_text()
        except (UnicodeDecodeError,ValueError):continue
        assert not any('\u0400'<=c<='\u04ff' for c in content),f'Cyrillic text: {path}'
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as out:
        for path in sorted(files):
            info=zipfile.ZipInfo(ROOT.name+'/'+path.relative_to(ROOT).as_posix(),date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(0o100755 if path.suffix=='.sh' else 0o100644)<<16
            out.writestr(info,path.read_bytes())
    print(json.dumps({'baseline_files_unchanged':len(names),'archive_files':len(files),'archive_bytes':destination.stat().st_size,'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}))

if __name__=='__main__':main(*sys.argv[1:])
