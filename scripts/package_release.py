"""Rebuild a complete SHA-256 manifest and one deterministic source archive."""
import hashlib
from pathlib import Path
import sys
import zipfile
ROOT=Path(__file__).resolve().parents[1]
EXCLUDED={'__pycache__','.pytest_cache','out','.git'}
def included(path):
    parts=path.relative_to(ROOT).parts
    return path.is_file() and not any(p in EXCLUDED or p.endswith('.egg-info') for p in parts) and path.suffix!='.pyc'
def main(output):
    paths=sorted(p for p in ROOT.rglob('*') if included(p) and p!=ROOT/'SHA256SUMS')
    manifest=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix()+'\n' for p in paths)
    (ROOT/'SHA256SUMS').write_text(manifest)
    paths.append(ROOT/'SHA256SUMS')
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in sorted(paths):
            info=zipfile.ZipInfo('windaudit-0.3.0/'+p.relative_to(ROOT).as_posix(),date_time=(2026,9,18,0,0,0))
            info.external_attr=(0o755 if p.suffix=='.sh' else 0o644)<<16
            z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    print(f'Packaged {len(paths)} files in {output}')
if __name__=='__main__':main(sys.argv[1])
