"""R2 GO-manifolds and R3 cycle geometry."""

from __future__ import annotations

import numpy as np
import pytest

from rookfields import PAPER
from rookfields import r2_manifolds as R2
from rookfields import r3_manifolds as R3
from rookfields.networks import INTRO_PERIODIC, RAMP_SYSTEMS, VAN_DER_POL
from rookfields.ramp import RampSystem

SYSTEMS = list(RAMP_SYSTEMS.values())
IDS = [s.name for s in SYSTEMS]
ADMISSIBLE_WIDTH = 0.005


def _system(spec):
    return RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))


# ---------------------------------------------------------------------------
# R2
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_go_pairs_exist(spec):
    pairs, _stg = R2.go_pairs(_system(spec), spec=PAPER)
    assert pairs, "expected codimension-two indecisive-drift pairs with GO-pairs"
    assert any(p.external for p in pairs)


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_control_region_is_degenerate_in_n_o(spec):
    """`Q^delta(xi, xi')` lies in the `n_o`-hyperplane.

    `n_o in J_i(xi) cap J_e(xi')`, so the ramp window and the adjacent plateau
    meet in exactly one point (GlobalDynR2v9.tex:196).  That is what makes the
    base `Nul_{n_o}^delta` have dimension `N-2` and its flow-out codimension one.
    """
    system = _system(spec)
    delta = 0.25 * R2.max_admissible_delta(system)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    for pair in pairs:
        region = R2.control_region(system, pair, delta)
        lo, hi = region[pair.n_opaque]
        assert abs(hi - lo) < 1e-12, f"pair {pair.face},{pair.coface}: {region}"


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_base_lies_on_the_nullcline(spec):
    system = _system(spec)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    for pair in pairs[:12]:
        manifold = R2.build_go_manifold(system, pair, base_samples=5)
        assert R2.base_is_on_nullcline(system, manifold) < 1e-9


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_admissible_epsilon_is_positive_and_needed(spec):
    """`defn:admissible-GO-perturbation` requires `r_{n_g} F_{eps,n_g} > 0`.

    Beyond the supremum the perturbation reverses the `n_g`-motion and the sweep
    never reaches the directed wall, so the bound is not cosmetic.
    """
    system = _system(spec)
    delta = 0.25 * R2.max_admissible_delta(system)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    external = [p for p in pairs if p.external]
    assert external
    for pair in external[:6]:
        bound = R2.max_admissible_epsilon(system, pair, delta)
        assert bound > 0

        good = R2.build_go_manifold(system, pair, delta=delta, epsilon=0.5 * bound)
        assert good.reached_wall

        bad = R2.build_go_manifold(
            system, pair, delta=delta, epsilon=4.0 * bound, max_time=5.0
        )
        assert not bad.reached_wall


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_external_go_manifolds_reach_the_directed_wall(spec):
    """`lem:GO-wall-hitting`: the flow-out attains the directed `n_g`-wall.

    The proof uses H2 to bound the crossing time; the width used here satisfies
    H2 (see the hypotheses report).
    """
    system = _system(spec)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    for pair in [p for p in pairs if p.external][:8]:
        manifold = R2.build_go_manifold(system, pair, base_samples=5)
        assert manifold.reached_wall
        assert np.all(
            np.abs(manifold.terminal[:, pair.n_grad] - manifold.directed_wall) < 1e-6
        )


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_go_crossing_sign(spec):
    """`prop:GO-crossing-sign`: `<F, z> = eps x_{n_g} |z_{n_g}| > 0` off the base.

    Equivalently the deviation `F - F_eps` is exactly `r eps x_{n_g} e^{(n_g)}`
    in the `(n_o, n_g)` block, which is the identity the proof uses.
    """
    system = _system(spec)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    for pair in pairs[:10]:
        manifold = R2.build_go_manifold(system, pair, base_samples=5)
        verdict = R2.crossing_sign(system, manifold)
        assert verdict["holds"], verdict


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_go_product_structure(spec):
    """`lem:GO-product`: the transverse directions are tangent to the sheet.

    In codimension two every direction outside `{n_o, n_g}` sits in a gap on
    both sides of the shared wall, so the base and its flow-out are products in
    those coordinates.
    """
    system = _system(spec)
    pairs, _ = R2.go_pairs(system, spec=PAPER)
    for pair in pairs[:10]:
        manifold = R2.build_go_manifold(system, pair, base_samples=5)
        verdict = R2.product_structure(system, manifold)
        assert verdict["holds"], verdict


# ---------------------------------------------------------------------------
# R3
# ---------------------------------------------------------------------------


def test_h3_bound_does_not_force_a_complex_pair():
    """A recorded defect, reproduced.

    `issues/source-defects.md` states that the positive-3-cycle `H_3` bound does
    not force a complex stable pair, with the explicit witness
    `(gamma_1,gamma_2,gamma_3) = (1,2,100)`, `P = 31000`.  `eq:F3-bounds` is
    equivalent to `|P| > D(gamma)` for the signed slope product, and that
    witness satisfies it while giving three real eigenvalues.
    """
    verdict = R3.check_h3_forces_complex_pair((1.0, 2.0, 100.0), 31000.0)
    assert verdict.satisfies_h3
    assert not verdict.complex_pair
    assert np.all(np.abs(np.imag(verdict.eigenvalues)) < 1e-9)


def test_real_spectrum_counterexamples_are_not_isolated():
    hits = R3.search_real_spectrum_counterexamples(limit=5)
    assert len(hits) >= 3
    for verdict in hits:
        assert verdict.satisfies_h3 and not verdict.complex_pair


def test_cyclic_feedback_matrix_matches_the_product_equation():
    """`det(lam I - J) = 0` iff `prod (lam + gamma_n) = P`."""
    gamma = (1.0, 2.0, 3.0)
    slopes = (2.0, -3.0, -5.0)
    cycle = (0, 1, 2)
    J = R3.cyclic_feedback_matrix(gamma, slopes, cycle)
    from_matrix = np.sort_complex(np.linalg.eigvals(J))
    P = float(np.prod(slopes))
    from_product = np.sort_complex(R3.eigenvalues_from_product(gamma, P))
    assert np.allclose(from_matrix, from_product, atol=1e-8)


def test_swept_surface_tangent_detection():
    """A generating direction of a swept surface is not transverse to it."""
    generator = lambda s: np.array([np.cos(2 * np.pi * s), np.sin(2 * np.pi * s), 0.0])
    flow = lambda _x: np.array([0.0, 0.0, 1.0])
    surface = R3.swept_surface(generator, flow, samples=16, steps=16)

    assert not R3.is_transverse(surface, [0.0, 0.0, 1.0])["transverse"]
    assert R3.is_transverse(surface, [1.0, 0.0, 0.0])["transverse"]


@pytest.mark.parametrize(
    "name,index,expected_lengths",
    [
        ("N3_B", 2_472_287, {3}),
        ("N3_A", 52_718_681_992, {2, 3}),
    ],
)
def test_cycle_cells_are_found(name, index, expected_lengths):
    """The R3 construction acts only on semi-opaque cells with a nontrivial cycle."""
    cells, _stg = R3.cycle_cells_for_parameter(name, index, spec=PAPER)
    assert cells, "expected semi-opaque cells with a nontrivial regulation cycle"
    assert {len(c.cycle) for c in cells} == expected_lengths
    for cell in cells:
        assert cell.signature in (-1, 1)
        assert cell.feedback in ("positive", "negative")


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_published_ramp_systems_have_no_regulation_cycles(spec):
    """Neither directly specified ramp system triggers the `F_3` cycle rule.

    A finding, not a limitation of this code: with no semi-opaque cell carrying
    a nontrivial cycle, Conditions 3.1 and 3.2 are vacuous, so `F_3 = F_2` for
    these parameters.  `ex:ramp_van_der_pol` is presented as an `F_3`
    computation, but the cycle rule never fires there.
    """
    from rookfields.blowup import SpecCubicalBlowupGraph

    system = _system(spec)
    assert R3.cycle_cells(system, spec=PAPER) == []

    labelling, K = system.wall_labelling()
    f2 = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=K, spec=PAPER, level=2
    )
    f3 = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=K, spec=PAPER, level=3
    )
    assert set(f2.digraph.edges()) == set(f3.digraph.edges())


# ---------------------------------------------------------------------------
# orange manifolds
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,index", [("N2_b", 974), ("N3_B", 2_472_287)])
def test_d3_triples_and_orange_manifolds(name, index):
    """`defn:orange-manifold` for the triples of `D_3(rook)`.

    A triple pairs an externally pruned coface with an internally pruned one
    over the same lower cell; the connector is the affine interpolation between
    the two anchors `Sigma + r_{n_g} eta e^{(n_g)}` and
    `Sigma + p_{n_o} alpha e^{(n_o)}`.
    """
    from rookfields.blowup import SpecCubicalBlowupGraph
    from rookfields.networks import realises_dsgrn_labelling
    from rookfields.ramp import RampSystem as _RS

    spec, realised = realises_dsgrn_labelling(name, index, 0.05)
    assert realised
    system = _RS(spec)
    labelling, K = system.wall_labelling()
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=K, spec=PAPER, level=2
    )
    triples = R2.d3_triples(stg)
    assert triples, "expected D_3 triples at this parameter node"

    for triple in triples[:6]:
        # the two cofaces really are pruned in opposite directions
        assert not stg.has_edge(triple.face, triple.external_coface)
        assert not stg.has_edge(triple.internal_coface, triple.face)

        orange = R2.build_orange_manifold(system, triple)
        assert orange.sigma.size
        assert orange.eta > 0 and orange.alpha > 0

        # Sigma lies on the codimension-two flat Pi(xi, xi'')
        for point in orange.sigma:
            assert np.isclose(
                point[triple.n_grad],
                R2.band_endpoint(system, triple.face_coords, triple.face_shape,
                                 triple.n_grad, triple.r_grad),
            )
            assert np.isclose(
                point[triple.n_opaque],
                R2.band_endpoint(system, triple.face_coords, triple.face_shape,
                                 triple.n_opaque, triple.p_opaque),
            )

        # the anchors are displaced by exactly eta and alpha along their axes
        d_g = orange.anchor_g - orange.sigma
        d_o = orange.anchor_o - orange.sigma
        assert np.allclose(d_g[:, triple.n_grad], triple.r_grad * orange.eta)
        assert np.allclose(d_o[:, triple.n_opaque], triple.p_opaque * orange.alpha)
        for n in range(system.dim):
            if n != triple.n_grad:
                assert np.allclose(d_g[:, n], 0.0)
            if n != triple.n_opaque:
                assert np.allclose(d_o[:, n], 0.0)

        # the sheet interpolates linearly between them
        sheet = orange.points(samples=5)
        assert np.allclose(sheet[:, 0, :], orange.anchor_g)
        assert np.allclose(sheet[:, -1, :], orange.anchor_o)
