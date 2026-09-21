import math
import numpy as np
import pytest
from windaudit.scanspace import (ScanPatch, PatchGraph, crossing_delta,
    theta_at, path_transport, strip_winding_delta, polyline_ijs, audit_edge_integrability)
from windaudit.scanspace_surface import Surface


def origin(z):
    return np.zeros(np.shape(z) + (2,))


def independent_ray_count(xyz):
    """Oriented intersections with positive x ray, no atan2 or unwrap."""
    result = 0
    for a, b in zip(xyz, xyz[1:]):
        x0, y0 = map(float, a[:2]); x1, y1 = map(float, b[:2])
        if (y0 < 0 <= y1) or (y1 < 0 <= y0):
            x = x0 + (x1 - x0) * (-y0) / (y1 - y0)
            if x > 0:
                result += 1 if y1 > y0 else -1
    return result


@pytest.mark.parametrize('turns', [-3, -1, 0, 1, 3])
def test_known_turns_against_ray_intersections(turns):
    a = np.linspace(.13, .13 + turns * 2 * math.pi, 601)
    xyz = np.column_stack([10*np.cos(a), 10*np.sin(a), np.ones(len(a))])
    assert crossing_delta(theta_at(xyz, origin)) == independent_ray_count(xyz) == turns


def test_random_small_steps_against_independent_ray():
    rng = np.random.default_rng(72)
    for _ in range(50):
        a = .123 + np.cumsum(rng.uniform(-.2, .2, 200))
        r = rng.uniform(5, 20, len(a))
        xyz = np.column_stack([r*np.cos(a), r*np.sin(a), np.ones(len(a))])
        assert crossing_delta(theta_at(xyz, origin)) == independent_ray_count(xyz)


def test_real_surface_routes_synthetic_known_seam():
    angle = np.linspace(-.5, .5, 20)
    radius = np.arange(10, 15)
    xyz = np.stack(np.broadcast_arrays(radius[:,None]*np.cos(angle),
                      radius[:,None]*np.sin(angle), np.ones((5,20))), axis=-1)
    patch = ScanPatch(xyz)
    for weight in [0, 4, 20]:
        graph = PatchGraph(patch, weight)
        assert strip_winding_delta(graph, [1.2,1.2], [2.2,17.2], origin)['delta_windings'] == 1
        assert strip_winding_delta(graph, [2.2,17.2], [1.2,1.2], origin)['delta_windings'] == -1


def test_closed_loop_can_fail_gate_on_surface_around_axis():
    y, x = np.mgrid[-4:5, -4:5]
    patch = ScanPatch(np.stack([x, y, np.ones_like(x)], axis=-1))
    loop = [[1.5,1.5],[1.5,6.5],[6.5,6.5],[6.5,1.5],[1.5,1.5]]
    result = path_transport(patch, loop, origin)
    assert result['delta_windings'] == 1
    lifted, valid = patch.lift(polyline_ijs(loop))
    assert valid.all()
    assert independent_ray_count(lifted) == 1


def test_disconnected_and_invalid_samples():
    y,x = np.mgrid[:6,:8]
    xyz = np.stack([x+10,y+10,np.ones_like(x)],axis=-1).astype(float)
    xyz[:,3:5] = -1
    patch = ScanPatch(xyz)
    assert strip_winding_delta(PatchGraph(patch), [1.2,1.2], [1.2,5.2], origin) is None
    _, valid = patch.lift([[-1,0],[1,3],[1,1]])
    assert valid.tolist() == [False,False,True]


def test_bilinear_lift_and_sampling():
    y,x = np.mgrid[:4,:4]
    patch = ScanPatch(np.stack([x, y, x*y], axis=-1))
    xyz,valid = patch.lift([[1.25,2.5]])
    np.testing.assert_allclose(xyz, [[2.5,1.25,3.125]])
    assert valid.all()
    pts=polyline_ijs([[0,0],[1,1],[2,1]])
    np.testing.assert_array_equal(pts, [[0,0],[.5,.5],[1,1],[2,1]])
    with pytest.raises(ValueError): polyline_ijs([[0,0],[1,1]],0)
    with pytest.raises(ValueError): theta_at([[0,0,1]], origin)


def test_triangle_projection_independent_plane_solution():
    y,x=np.mgrid[:6,:7]
    patch=ScanPatch(np.stack([x,y,np.full_like(x,10)],axis=-1))
    surface=Surface(patch)
    rng=np.random.default_rng(3)
    for _ in range(60):
        point=rng.uniform([.1,.1,8],[5.9,4.9,12])
        hit=surface.project(point,2.5)
        assert hit is not None
        np.testing.assert_allclose(hit['ij'],point[1::-1],atol=1e-5)
        assert hit['distance']==pytest.approx(abs(point[2]-10),abs=2e-6)
    assert surface.project([2,2,14],2.5) is None


def test_edge_gate_detects_hole_holonomy():
    y,x=np.mgrid[-5:6,-5:6]
    xyz=np.stack([x,y,np.ones_like(x)],axis=-1).astype(float)
    xyz[4:7,4:7]=-1
    r=audit_edge_integrability(ScanPatch(xyz),origin)
    assert r['inconsistent_edges']>0
    assert r['max_abs_residual']==1


def test_edge_gate_accepts_simply_connected_off_axis_patch():
    y,x=np.mgrid[:7,:9]
    patch=ScanPatch(np.stack([x+20,y-3,np.ones_like(x)],axis=-1))
    assert audit_edge_integrability(patch,origin)['inconsistent_edges']==0
