"""CPU-only theta=0 transport on tifxyz valid-quad paths.

This is an additive experiment, never used by the legacy annotation exporter.
Routing, float32 polyline sampling and bilinear interpolation follow villa
2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7. Scan-space theta is not a claim
of equivalence to every fitted deformation. Test path independence first.
"""
from dataclasses import dataclass, field
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import distance_transform_edt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.csgraph import breadth_first_order


def load_umbilicus(path):
    points = sorted(json.loads(Path(path).read_text())['control_points'],
                    key=lambda p: p['z'])
    return interp1d(np.asarray([p['z'] for p in points], dtype=np.float32),
                    np.asarray([[p['y'], p['x']] for p in points], dtype=np.float32),
                    axis=0, fill_value='extrapolate')


def theta_at(xyz, umbilicus):
    xyz = np.asarray(xyz, dtype=np.float32)
    yx = xyz[..., 1::-1] - umbilicus(xyz[..., 2])
    if not np.isfinite(yx).all() or np.any(np.linalg.norm(yx, axis=-1) == 0):
        raise ValueError('Undefined scan-space angle')
    return np.arctan2(yx[..., 0], yx[..., 1]) % (2 * np.pi)


def crossing_delta(theta):
    theta = np.asarray(theta)
    if not np.isfinite(theta).all():
        raise ValueError('Non-finite angle')
    differences = np.diff(theta)
    return int(np.count_nonzero(differences < -np.pi)
               - np.count_nonzero(differences > np.pi))


def polyline_ijs(points, step_size=1.0):
    if not np.isfinite(step_size) or step_size <= 0:
        raise ValueError('step_size must be finite and positive')
    points = np.asarray(points, dtype=np.float32)
    segments = []
    for a, b in zip(points, points[1:]):
        n = max(1, int(np.ceil(float(np.linalg.norm(b - a)) / step_size))) + 1
        t = np.linspace(0, 1, n, dtype=np.float32)[:, None]
        segment = a[None] * (1 - t) + b[None] * t
        segments.append(segment[1:] if segments else segment)
    return np.concatenate(segments).astype(np.float32) if segments else points.copy()


@dataclass
class ScanPatch:
    xyz: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.xyz = np.asarray(self.xyz, dtype=np.float32)
        if self.xyz.ndim != 3 or self.xyz.shape[-1] != 3:
            raise ValueError('Expected H x W x 3 coordinates')
        valid = np.any(self.xyz != -1, axis=-1)
        self.valid_quad = (valid[:-1, :-1] & valid[1:, :-1]
                           & valid[:-1, 1:] & valid[1:, 1:])

    @classmethod
    def load(cls, directory):
        # Pillow is only needed for the optional real-patch experiment.
        from PIL import Image
        directory = Path(directory)
        xyz = np.stack([np.asarray(Image.open(directory / (c + '.tif')))
                        for c in 'xyz'], axis=-1).astype(np.float32)
        if (directory / 'mask.tif').exists():
            mask = np.asarray(Image.open(directory / 'mask.tif'))
            if mask.ndim == 3:
                mask = mask[..., 0]
            xyz[mask == 0] = -1
        return cls(xyz, json.loads((directory / 'meta.json').read_text()))

    def lift(self, ij):
        ij = np.asarray(ij, dtype=np.float32)
        i, j = np.floor(ij).astype(int).T
        h, w = self.valid_quad.shape
        valid = (i >= 0) & (j >= 0) & (i < h) & (j < w)
        valid[valid] &= self.valid_quad[i[valid], j[valid]]
        result = np.full((len(ij), 3), -1, dtype=np.float32)
        ii, jj = i[valid], j[valid]
        di, dj = (ij[valid] - np.floor(ij[valid])).T
        tl, tr = self.xyz[ii, jj], self.xyz[ii, jj + 1]
        bl, br = self.xyz[ii + 1, jj], self.xyz[ii + 1, jj + 1]
        top = tl + (tr - tl) * dj[:, None]
        bottom = bl + (br - bl) * dj[:, None]
        result[valid] = top + (bottom - top) * di[:, None]
        return result, valid


class PatchGraph:
    def __init__(self, patch, medial_weight=4.0, random_seed=None):
        if medial_weight < 0:
            raise ValueError('Negative medial weight')
        self.patch = patch
        v = patch.valid_quad
        self.ii, self.jj = np.nonzero(v)
        self.n = len(self.ii)
        if not self.n:
            raise ValueError('Patch has no valid quads')
        self.centres = np.column_stack((self.ii, self.jj)).astype(np.float32) + .5
        self.ids = np.full(v.shape, -1, dtype=np.int64)
        self.ids[v] = np.arange(self.n)
        dt = distance_transform_edt(np.pad(v, 1))[1:-1, 1:-1].astype(np.float32)
        rows, cols, data = [], [], []
        rng = np.random.default_rng(random_seed)
        for di, dj in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]:
            ni, nj = self.ii + di, self.jj + dj
            ok = (ni >= 0) & (nj >= 0) & (ni < v.shape[0]) & (nj < v.shape[1])
            src = np.flatnonzero(ok)
            dst = self.ids[ni[ok], nj[ok]]
            ok = dst >= 0
            src, dst = src[ok], dst[ok]
            weights = np.hypot(di, dj) * (1 + medial_weight / np.maximum(dt[self.ii[dst], self.jj[dst]], 1))
            if random_seed is not None:
                weights = weights * np.exp(rng.uniform(-2, 2, len(weights)))
            rows.extend(src); cols.extend(dst); data.extend(weights)
        self.csr = coo_matrix((data, (rows, cols)), shape=(self.n, self.n)).tocsr()
        self.cache = {}

    def nearest_node(self, ij):
        h, w = self.ids.shape
        qi, qj = np.clip(np.floor(ij).astype(int), [0, 0], [h - 1, w - 1])
        if self.ids[qi, qj] >= 0:
            return int(self.ids[qi, qj])
        return int(np.argmin((self.ii - qi) ** 2 + (self.jj - qj) ** 2))

    def route(self, start, end):
        a, b = self.nearest_node(start), self.nearest_node(end)
        if a not in self.cache:
            self.cache[a] = dijkstra(self.csr, directed=True, indices=a,
                                      return_predecessors=True)[1]
        pred = self.cache[a]
        route = [b]
        while route[-1] != a:
            previous = int(pred[route[-1]])
            if previous < 0:
                return None
            route.append(previous)
        return self.centres[route[::-1]]


def path_transport(patch, points, umbilicus, step_size=1.0):
    ij = polyline_ijs(points, step_size)
    xyz, valid = patch.lift(ij)
    if valid.sum() < 2:
        return None
    return {'delta_windings': crossing_delta(theta_at(xyz[valid], umbilicus)),
            'strip_num_points': int(valid.sum()),
            'strip_num_invalid_dropped': int((~valid).sum())}


def strip_winding_delta(graph, from_ij, to_ij, umbilicus, step_size=1.0):
    centres = graph.route(from_ij, to_ij)
    if centres is None:
        return None
    result = path_transport(graph.patch, np.vstack([from_ij, centres, to_ij]),
                            umbilicus, step_size)
    if result is not None:
        result['strip_num_path_quads'] = len(centres)
    return result


def audit_edge_integrability(patch, umbilicus):
    """Check every fundamental cycle of the sampled quad-centre graph.

    Fixed step 1: orthogonal centre edges have endpoints only, diagonals
    include the shared corner. Match upstream's invalid-sample dropping.
    A zero residual on every edge certifies all centre-graph closed walks
    for this sampling rule, including holes and disconnected components.
    """
    graph = PatchGraph(patch)
    centres=graph.centres
    xyz,valid=patch.lift(centres)
    assert valid.all()
    theta=theta_at(xyz,umbilicus)
    row,col=graph.csr.nonzero()
    keep=row<col;row=row[keep];col=col[keep]
    d=theta[col]-theta[row]
    delta=(d < -np.pi).astype(int)-(d > np.pi).astype(int)
    diagonal=np.all(centres[row]!=centres[col],axis=1)
    indices=np.flatnonzero(diagonal)
    mid=(centres[row[indices]]+centres[col[indices]])/2
    xyz,valid=patch.lift(mid)
    k=indices[valid];tm=theta_at(xyz[valid],umbilicus)
    d1=tm-theta[row[k]];d2=theta[col[k]]-tm
    delta[k]=(d1 < -np.pi).astype(int)-(d1 > np.pi).astype(int)+(d2 < -np.pi).astype(int)-(d2 > np.pi).astype(int)
    dg=coo_matrix((np.r_[delta,-delta],(np.r_[row,col],np.r_[col,row])),shape=(graph.n,graph.n)).tocsr()
    potentials=np.full(graph.n, np.iinfo(np.int64).max,dtype=np.int64)
    sentinel=np.iinfo(np.int64).max
    components=0
    for start in range(graph.n):
        if potentials[start]!=sentinel:continue
        components+=1;order,pred=breadth_first_order(graph.csr,start,directed=False)
        potentials[start]=0
        for node in order[1:]:potentials[node]=potentials[pred[node]]+int(dg[pred[node],node])
    residual=potentials[col]-potentials[row]-delta
    return {'valid_quads':graph.n,'components':components,'edges':len(row),
            'inconsistent_edges':int(np.count_nonzero(residual)),
            'max_abs_residual':int(np.max(np.abs(residual))) if len(residual) else 0,
            'invalid_diagonal_midpoints':int((~valid).sum())}
