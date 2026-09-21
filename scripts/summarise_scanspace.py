"""Build the additive English report from measured evidence."""
import hashlib,json,platform
from pathlib import Path
import numpy,scipy,PIL
ROOT=Path(__file__).resolve().parents[1]
def main():
    transport=json.loads((ROOT/'results/scanspace_transport.json').read_text())
    report=json.loads((ROOT/'results/scanspace/real_graph.json').read_text())
    band=json.loads((ROOT/'results/scanspace/band_transport_gate.json').read_text())
    transport['selected_band_gate']={'patches':len(band),'edges':sum(x['edges'] for x in band),
        'cycle_rank':sum(x['edges']-x['valid_quads']+x['components'] for x in band),
        'inconsistent_edges':sum(x['inconsistent_edges'] for x in band),
        'step_size':1.0,'details':'scanspace/band_transport_gate.json'}
    transport['route_count_note']='distinct_alternative_routes counts routes different from the baseline, not pairwise unique alternatives.'
    transport['independent_synthetic_tests']='tests/test_scanspace.py: analytic turn counts, direct positive-ray intersection oracle, surface routes, plane projection, and a failing annular surface.'
    transport['real_graph_report']='scanspace/real_graph.json'
    transport['real_data_solver_improvement_demonstrated']=report['real_data_improvement_demonstrated']
    transport['same_149_same_151_sheet_identity']='unresolved: target surface coverage is incomplete'
    transport['task_evidence_status']='transport gate passed; real-data non-regression measured; improvement and target sheet identity unproven'
    (ROOT/'results/scanspace_transport.json').write_text(json.dumps(transport,indent=2)+'\n')
    environment={'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,'pillow':PIL.__version__,'device':'CPU','torch_required':False,'checkpoint_required':False}
    (ROOT/'results/scanspace_environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    hashes={}
    for folder in ['scanspace_probe','scanspace_band']:
        for p in sorted((ROOT/'data'/folder).rglob('*')):
            if p.is_file():hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    (ROOT/'results/scanspace_input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    lines=['# Scan-space winding transport','',
        'This additive experiment uses CPU NumPy/SciPy and actual public PHercParis4 tifxyz surfaces. No checkpoint or torch is required. The legacy README, annotation results and exporters are unchanged.','',
        '## Transport gate','',
        f"The first experiment tested {transport['pairs']:,} real-patch endpoint pairs, five weighted-Dijkstra routes per pair, at the pinned upstream default step of 1 grid cell. Disagreeing pairs: {transport['disagreeing_pairs']}; disagreement fraction: {transport['disagreement_fraction']:.6f}. Non-zero closed walks: {transport['nonzero_closed_contours']} of {transport['closed_contours']:,}. Sampling was never tuned.", '',
        'The five weight modes were medial weights 4, 0 and 20, followed by two seeded positive random perturbations of weight 4. Endpoints and the valid 8-neighbour quad graph stayed fixed. Each alternative was compared with the baseline; the alternative plus the reversed baseline formed the closed walk. The saved count of alternatives is the number different from baseline.','',
        f"An additional edge-potential check on {len(band)} band surfaces checked all {sum(x['edges'] for x in band):,} undirected sampled centre edges. Non-zero residuals: {sum(x['inconsistent_edges'] for x in band)}. Its cycle-space rank is {sum(x['edges']-x['valid_quads']+x['components'] for x in band):,}. This certifies all centre-graph closed walks on these surfaces for the fixed sampling rule; unseen surfaces still require the gate.",'',
        'Both bilinear lifting and upstream invalid-sample dropping are preserved. Diagnostic counts explicitly expose dropped samples. Synthetic tests also contain a surface enclosing the axis whose closed walk has winding +1 and must fail the gate.','',
        '## Data selection','',
        'All 102 points in same:149, same:150 and same:151 have z = 12418.1259765625. The half-open band [12417, 12419) contains every one of them. Spatial coverage by an actual surface is checked separately.','',
        'The verified_patches listing contains 89,237 folders and exposes no top-level bounding-box manifest. Individual meta.json files provide bboxes; upstream warns that historical producers used different axis orders. The initial selection inspected 341 named non-auto, non-band-seed, non-same_wrap patches and selected 42 conservatively. Their coordinate planes, masks and metadata total 10,548,497 bytes. The entire patch tree was never downloaded. This is a bounded subset, without a claim of globally minimal or exhaustive selection.','',
        'Additional targeted metadata discovery brought the persisted metadata census to 689 unique patches. Six supplementary candidates totalled 2,623,401 bytes. Their raw surfaces also passed the fixed-step edge-integrability gate before use. The final data directory has 48 candidate patches (13,171,898 bytes); coordinate validation after erosion retains 46. Two bbox candidates fail the actual quad-intersection check. The original 42-patch experiment is preserved under results/scanspace_initial. The supplementary candidates did not supply any surface hits for same:150 or same:151.','',
        'Actual coordinates and masks are checked after download. The CPU experiment retains surfaces whose valid quad extents intersect the band. This explicitly differs from upstream load_tifxyz\'s vertex-only z prefilter, which can miss a two-voxel band between grid rows. Default upstream one-cell binary erosion and per-patch overrides are respected.','',
        '## Real graph and solvers','',
        f"{report['total_source_points']:,} original annotation points were offered to the selected surfaces; {report['attached_points']} attached at tolerance 2.5 voxels. {report['cross_collections']} collections produced {report['reached_patches']} reached patches and {report['measurable_edges']} measurable undirected edges. Unmeasurable edges: {report['unmeasurable_edges']}.",'',
        'Attachment uses the pinned Patch.project triangular fallback on actual valid faces, conservative KD-tree pruning, and largest-area-then-nearest selection. Equivalence to the separate native vc_spiral surface index has not been established. Each retained connected component gets a BFS entry. The upstream build_rel_adjacency function is executed from the pinned source AST with scan-space tour adjustments. The before and after solve_min_edge_fix functions are extracted from the original and the mechanically patched pinned source.','',
        '| Mode | Version | Edges considered | Proposed changes | Retained graph consistent | Seconds | Status |','|---|---|---:|---:|---|---:|---|']
    for mode,key in [('Inconsistent-cycle candidates','solver_comparison'),('All edges editable','unrestricted_solver_comparison')]:
        for version,r in report[key].items():
            lines.append(f"| {mode} | {version} | {r['num_edges_considered']} | {r['num_edges_changed']} | {r['independent_retained_graph_check']['consistent']} | {r['wall_seconds']:.6f} | {r['solver_status']} |")
    lines+=['',f"Real-data improvement demonstrated: **{report['real_data_improvement_demonstrated']}**. This data provides a real-surface non-regression experiment. The CONTRIBUTING requirement for a demonstrated real-data improvement remains open when both solvers retain a consistent graph with no repairs.",'',
        'The independent verifier ignores the MILP labels and traverses every retained edge with exact integer potentials. Suggested detach_or_reattach actions are represented only as removed graph edges; actual annotation editing and re-linking have not been performed.','',
        '## The three traces','',
        '| Trace | Points | Points with an actual surface hit within 2.5 voxels |','|---|---:|---:|']
    for t in report['target_traces']:
        lines.append(f"| {t['collection']} | {len(t['points'])} | {sum(bool(p['hits']) for p in t['points'])} |")
    lines+=['','The per-point hit list, distances and patch IDs are saved in real_graph.json. A missing hit leaves sheet identity unresolved. Merely including a trace\'s z coordinate in the band does not establish surface coverage or decide whether same:149 and same:151 occupy the same sheet.','',
        '## Reproduction and preservation','',
        'Install the existing dev/figures extras (Pillow reads TIFFs), then run:','',
        '```bash','bash scripts/reproduce.sh','bash scripts/reproduce_scanspace.sh','python experiments/scanspace/gate.py','python experiments/scanspace/cycle_probe.py','sha256sum -c SHA256SUMS','```','',
        'The long five-route gate uses the bundled three-patch probe corpus; it writes under out/scanspace_gate. The real graph reproduction runs its exhaustive edge gate before projection or solver calls. All downloaded coordinate files needed for these runs are bundled and hashed.','',
        'results/tests.txt records the full legacy test suite, two reproduction runs, deterministic report/export checks and package checks. results/scanspace_preservation.json compares 101 old files byte for byte against the supplied archive; only SHA256SUMS and the refreshed tests log are exempt, and their originals are retained under results/legacy_borderline.','',
        'upstream/solver-candidate.patch is unchanged. upstream/scanspace-combined.patch adds the standalone CPU files alongside that existing fix. The upstream collection enumeration is already correct; no collection-name collision fix is included or claimed. For the standalone upstream audit, pass the original solver source, patched solver source, output directory and annotation input directory after the patch directory argument.','',
        'Source: ScrollPrize/villa commit 2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7. Public data: https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/ . Detailed source and input hashes accompany the evidence.','']
    (ROOT/'SCANSPACE_TRANSPORT.md').write_text('\n'.join(lines))
if __name__=='__main__':main()
