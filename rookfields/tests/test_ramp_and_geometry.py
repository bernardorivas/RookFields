"""Ramp systems, width conditions, geometrization, and alignment."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from rookfields import LEGACY, PAPER, conley_morse_graph
from rookfields import hypotheses as H
from rookfields.alignment import check_alignment
from rookfields.geometrization import RectangularGeometrization, JanusComplex
from rookfields.networks import INTRO_PERIODIC, RAMP_SYSTEMS, VAN_DER_POL
from rookfields.ramp import RampSystem
from rookfields.wall_labeling import WallLabeling

SYSTEMS = list(RAMP_SYSTEMS.values())
IDS = [s.name for s in SYSTEMS]

#: A width satisfying H_0, H_1, H_2 and H_3 for both published systems.
ADMISSIBLE_WIDTH = 0.005


# ---------------------------------------------------------------------------
# the ramp cell complex reproduces the index sets printed in the text
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_index_set_matches_the_text(spec):
    assert spec.claimed_index_set is not None
    assert spec.index_set == spec.claimed_index_set


# ---------------------------------------------------------------------------
# defn:ramp-wall-labeling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_wall_labelling_is_valid_and_dissipative(spec):
    """A ramp induced wall labelling must satisfy `thm:ramp-wall-labeling`.

    `def:wall_labeling` validity plus strong dissipativity, the standing
    assumption from `CombinatorialDynamics.tex:4` onward.
    """
    system = RampSystem(spec)
    labelling, K = system.wall_labelling()
    wl = WallLabeling(labelling, K)
    assert wl.is_valid()
    assert wl.is_strongly_dissipative()


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
@pytest.mark.parametrize("h", [0.1, 0.05, 0.02, 0.005, 0.002])
def test_wall_labelling_is_constant_in_h(spec, h):
    """`prop:wall-labeling-const` (RampSystemsv4.tex:490).

    The labelling does not depend on `h` inside `H_0`.  This holds for the
    `GB_n` sentinel of `defn:ramp-wall-labeling`; it fails for the
    `theta + 10h` surrogate in `WallLabelling.py:44`, which is what
    `test_legacy_sentinel_breaks_dissipativity` pins down.
    """
    reference, K0 = RampSystem(spec).wall_labelling()
    shrunk, K1 = RampSystem(spec.with_uniform_width(h)).wall_labelling()
    assert K0 == K1
    assert shrunk == reference


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_legacy_sentinel_breaks_dissipativity_at_small_width(spec):
    """The `theta + 10h` outer sentinel loses strong dissipativity as `h -> 0`.

    `defn:ramp-wall-labeling` sets the outermost sentinel threshold to `GB_n`.
    `DSGRN_utils.WallLabelling.ramp_system_wall_labelling` (:44) uses
    `theta_{m_K} + 10 h_{m_K}` instead, which shrinks with `h` and eventually
    falls inside the region it is supposed to bound.  The failure is invisible
    at the published widths, which is presumably why it survived.
    """
    published = RampSystem(spec)
    labelling, K = published.wall_labelling(outer="legacy")
    assert WallLabeling(labelling, K).is_strongly_dissipative()

    small = RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))
    legacy_labelling, legacy_K = small.wall_labelling(outer="legacy")
    assert not WallLabeling(legacy_labelling, legacy_K).is_strongly_dissipative()

    correct_labelling, correct_K = small.wall_labelling(outer="global_bound")
    assert WallLabeling(correct_labelling, correct_K).is_strongly_dissipative()


# ---------------------------------------------------------------------------
# the width conditions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_published_widths_violate_H1(spec):
    """The published ramp parameters do not satisfy `H_1`.

    Recorded as a finding, not a defect of this code: `defn:H1` forbids the
    focal value `E_n(D_v)/gamma_n` from lying in a ramp half-window, and at the
    published widths it does.  `H_2` and `H_3` are subsets of `H_1`, so they
    fail too, and with them the hypotheses of `thm:R1ABlattice` and
    `thm:R3ABlattice`.
    """
    system = RampSystem(spec)
    assert H.check_admissible_theta_h(system).holds, "H_0 should hold"
    assert H.check_lambda_R(system).holds, "Lambda(R) should hold"
    assert not H.check_H1(system).holds


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_small_width_satisfies_every_condition(spec):
    system = RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))
    for verdict in (
        H.check_admissible_theta_h(system),
        H.check_lambda_R(system),
        H.check_H1(system),
        H.check_H2(system),
        H.check_H3(system),
    ):
        assert verdict.holds, str(verdict)


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_shrinking_h_preserves_the_morse_graph(spec):
    """Repairing the width does not disturb the published combinatorics.

    So the computed Morse graphs stand; only the stated widths need changing.
    """
    published = RampSystem(spec)
    repaired = RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))
    a = conley_morse_graph(
        labelling=published.wall_labelling()[0],
        num_thresholds=published.num_thresholds,
        spec=PAPER,
        level=3,
    )
    b = conley_morse_graph(
        labelling=repaired.wall_labelling()[0],
        num_thresholds=repaired.num_thresholds,
        spec=PAPER,
        level=3,
    )
    assert a.conley_indices == b.conley_indices
    assert a.edges == b.edges


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_prop_smallh_bound_really_lands_in_H1(spec):
    """`prop:smallh`'s constructive bound gives a width inside `H_1`."""
    system = RampSystem(spec)
    bound = H.uniform_width_bound(system, 1)
    assert bound > 0
    shrunk = RampSystem(spec.with_uniform_width(0.99 * bound))
    assert H.check_admissible_theta_h(shrunk).holds
    assert H.check_H1(shrunk).holds


# ---------------------------------------------------------------------------
# rectangular geometrization
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_rectangular_geometrization_fills_the_bounding_box(spec):
    """`eq:rectangular_geo_boundary`: the embedded complex is exactly `prod [0, GB_n]`."""
    geo = RectangularGeometrization(RampSystem(spec))
    assert geo.covers_bounding_box()
    assert geo.intervals_are_ordered()


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_blowup_cell_image_matches_I_n(spec):
    """`rem:g(b(xi))` agrees with `I_n(xi)` away from the sentinel indices."""
    system = RampSystem(spec)
    geo = RectangularGeometrization(system)
    compared = 0
    for v in itertools.product(*[range(1, k + 2) for k in system.num_thresholds]):
        for w in itertools.product(*[[0, 1]] * system.dim):
            if geo.uses_sentinel(list(v), list(w)):
                continue
            compared += 1
            a = geo.blowup_cell_image_of(list(v), list(w))
            b = [system.interval(list(v), list(w), n) for n in range(system.dim)]
            assert np.allclose(a, b), (v, w, a, b)
    assert compared > 0


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_chart_restricts_to_faces(spec):
    """`a_zeta|_{K(u)} = a_{[v,u]}` for `u <_Z w`, the inductive requirement."""
    geo = RectangularGeometrization(RampSystem(spec))
    dim = geo.system.dim
    coords = [2] * dim
    full = [1] * dim
    a_full = geo.chart(coords, full)
    for drop in range(dim):
        sub = list(full)
        sub[drop] = 0
        a_sub = geo.chart(coords, sub)
        for y in itertools.product([0.0, 0.5, 1.0], repeat=dim - 1):
            # embed the sub-face parameter into the full cube at y_drop = 0
            full_y = list(y)
            full_y.insert(drop, 0.0)
            assert np.allclose(a_full(full_y), a_sub(y))


# ---------------------------------------------------------------------------
# Janus complex
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_janus_adjacent_top_slices_are_disjoint(spec):
    """Regression for the corrected `eq:barS`.

    The earlier carrier assigned the fine top cell based at `8(v+1)` to both
    adjacent coarse top cells, which is what made `pi_J` and the top-cell branch
    of `pi_2` ill defined.
    """
    janus = JanusComplex(RampSystem(spec))
    assert janus.adjacent_top_slices_are_disjoint(limit=200)


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_janus_subdivision_endpoint(spec):
    """A coarse edge contains exactly the eight fine edges based at `8v .. 8v+7`."""
    janus = JanusComplex(RampSystem(spec))
    assert janus.subdivision_chain_endpoint_ok()
    system = RampSystem(spec)
    assert janus.vertex_range == [
        8 * (2 * (k + 1) + 2) for k in system.num_thresholds
    ]


# ---------------------------------------------------------------------------
# alignment -- thm:R1ABlattice
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
@pytest.mark.parametrize("algo", [LEGACY, PAPER], ids=lambda s: s.name)
def test_rectangular_geometrization_is_aligned_with_F1(spec, algo):
    """`thm:R1ABlattice` (GlobalDynR1v3.tex:15-21), checked numerically.

    Every rectangular geometrization is aligned with the ramp system over all
    `N in N(F_1)`, once the width condition holds.
    """
    system = RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))
    report = check_alignment(system, spec=algo, level=1, samples_per_axis=3)
    assert report.walls, "no oriented walls to check"
    assert report.ok, str(report) + f"\nfirst failures: {report.failures[:3]}"


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_F1_alignment_needs_the_width_condition(spec):
    """The hypothesis of `thm:R1ABlattice` is doing real work.

    At the published widths -- which violate `H_1` -- the rectangular
    geometrization is *not* aligned.
    """
    report = check_alignment(RampSystem(spec), spec=PAPER, level=1, samples_per_axis=3)
    assert not report.ok


@pytest.mark.parametrize("spec", SYSTEMS, ids=IDS)
def test_F2_misalignments_are_exactly_the_newly_oriented_walls(spec):
    """Where the rectangular geometry stops sufficing.

    `thm:R1ABlattice` covers `F_1` only; `F_2`/`F_3` orient walls that `F_1`
    left as double edges, and on precisely those the rectangular face has the
    wrong crossing direction.  Those are the faces the R2 GO-manifold and R3
    cycle-surface constructions are meant to replace.
    """
    from rookfields.blowup import SpecCubicalBlowupGraph

    system = RampSystem(spec.with_uniform_width(ADMISSIBLE_WIDTH))
    labelling, K = system.wall_labelling()
    f1 = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=K, spec=PAPER, level=1
    )
    f1_edges = set(f1.digraph.edges())

    report = check_alignment(system, spec=PAPER, level=2, samples_per_axis=3)
    assert report.failures, "expected the rectangular geometry to fall short for F_2"
    for wall in report.failures:
        assert (wall.source, wall.target) in f1_edges
        assert (wall.target, wall.source) in f1_edges
