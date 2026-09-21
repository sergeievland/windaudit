"""Compose the existing solver fix and an additive standalone CPU audit tool."""
import difflib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    result=(ROOT/'upstream/solver-candidate.patch').read_text()
    mappings=[('windaudit/scanspace.py','spiral-fitting/scanspace_transport.py'),
              ('windaudit/scanspace_surface.py','spiral-fitting/scanspace_surface.py'),
              ('scripts/scanspace_real_graph.py','spiral-fitting/audit_scanspace_transport.py')]
    for source,target in mappings:
        text=(ROOT/source).read_text()
        if source.startswith('scripts/'):
            text=text.replace('from windaudit.scanspace import','from scanspace_transport import').replace('from windaudit.scanspace_surface import','from scanspace_surface import')
        result+=''.join(difflib.unified_diff([],text.splitlines(keepends=True),fromfile='/dev/null',tofile='b/'+target))
    (ROOT/'upstream/scanspace-combined.patch').write_text(result)
if __name__=='__main__':main()
