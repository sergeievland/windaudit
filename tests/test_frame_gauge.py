"""The identity behind frame invariance of winding contradictions.

For a path sampled finely enough that every step turns by less than pi in both
frames, the signed branch-crossing count in a rotated frame equals the count in
the original frame plus k(end) - k(start), where k is an integer function of
position alone. Closed walks therefore keep their count, and graph equations
change only by gauge differences.
"""

import numpy as np
import pytest

from windaudit.scanspace import crossing_delta


def admissible_twist(rng):
    a, b, c = rng.uniform(-2, 2), rng.uniform(-2, 2), rng.uniform(-0.8, 0.8)
    ph = rng.uniform(0, 2 * np.pi)

    def alpha(theta, r, z):  # |d alpha / d theta| < 1: a homeomorphism of each circle
        return a * np.sin(z / 50 + ph) + b * np.cos(r / 7) + c * np.sin(theta)
    return alpha


def sample_path(rng, n):
    """A smooth random curve in (r, theta, z) staying away from the axis."""
    t = np.linspace(0, 1, n)
    theta = rng.uniform(0, 2 * np.pi) + rng.uniform(-9, 9) * t + 0.6 * np.sin(rng.uniform(1, 6) * t)
    r = 20 + 5 * np.sin(rng.uniform(1, 5) * t + rng.uniform(0, 6))
    z = 100 * t + 10 * np.cos(3 * t)
    return theta, r, z


@pytest.mark.parametrize("seed", range(40))
def test_rotated_count_differs_by_endpoint_gauge(seed):
    rng = np.random.default_rng(seed)
    alpha = admissible_twist(rng)
    theta, r, z = sample_path(rng, 4000)
    base = np.mod(theta, 2 * np.pi)
    lifted = theta + alpha(np.mod(theta, 2 * np.pi), r, z)
    rotated = np.mod(lifted, 2 * np.pi)
    # k(p) = floor((theta(p) + alpha(p)) / 2 pi): an integer function of the point alone
    k = np.floor((base + alpha(base, r, z)) / (2 * np.pi)).astype(int)
    assert crossing_delta(rotated) == crossing_delta(base) + k[-1] - k[0]


@pytest.mark.parametrize("seed", range(40))
def test_closed_walk_count_is_frame_independent(seed):
    rng = np.random.default_rng(1000 + seed)
    alpha = admissible_twist(rng)
    theta, r, z = sample_path(rng, 4000)
    # close the walk: return to the start point along the reversed path
    theta = np.concatenate([theta, theta[::-1]])
    r = np.concatenate([r, r[::-1]])
    z = np.concatenate([z, z[::-1]])
    base = np.mod(theta, 2 * np.pi)
    rotated = np.mod(theta + alpha(base, r, z), 2 * np.pi)
    assert crossing_delta(rotated) == crossing_delta(base) == 0


def test_untransported_gap_breaks_closure():
    """The attachment gap: an annotation point at theta = 6.28 whose surface
    attachment sits at theta = 0.003, just across the branch ray. Transported as
    one closed walk the count is zero. Counting the two pieces separately — the
    tour at the annotation point, the strip from the attachment — drops the step
    across the ray and leaves a spurious one-winding residual."""
    loop = np.array([6.20, 6.28, 0.003, 0.05, 6.27, 6.20])
    assert crossing_delta(loop) == 0
    tour_piece = loop[:2]    # ends at the annotation point
    strip_piece = loop[2:]   # starts at the attachment
    assert crossing_delta(tour_piece) + crossing_delta(strip_piece) == -1
