"""Combinatorial claims from Parts II of the monograph, tested by computation.

These are statements the manuscript proves; a failure here is a finding about
the mathematics or the implementation, not a broken test fixture.
"""

from __future__ import annotations

import pytest

from rookfields import LEGACY, PAPER, SpecCubicalBlowupGraph, conley_morse_graph
from rookfields import networks
from rookfields.wall_labeling import enumerate_wall_labelings

SPECS = [LEGACY, PAPER]


def _double_edges(stg):
    """Unordered pairs of distinct non-fringe cells joined in both directions."""
    edges = set(stg.digraph.edges())
    return {
        frozenset((a, b))
        for (a, b) in edges
        if a != b
        and (b, a) in edges
        and not stg.blowup_complex.rightfringe(a)
        and not stg.blowup_complex.rightfringe(b)
    }


# ---------------------------------------------------------------------------
# thm:2dDone -- CombinatorialDynamics.tex:1135
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
@pytest.mark.parametrize("num_thresholds", [[1, 1], [1, 2]])
def test_2d_done_exhaustive(spec, num_thresholds):
    """`F_3` has no double edges, for *every* abstract 2D wall labeling.

    The theorem is stated for an arbitrary wall labeling on a two-dimensional
    cubical complex, so it is tested against `def:wall_labeling`-valid,
    strongly dissipative labelings enumerated directly rather than against
    DSGRN-derived ones only.
    """
    labelings = list(enumerate_wall_labelings(num_thresholds))
    assert labelings, "enumerator produced nothing"
    for wl in labelings:
        stg = SpecCubicalBlowupGraph(
            labelling=wl.labelling, num_thresholds=num_thresholds, spec=spec, level=3
        )
        assert not _double_edges(stg), f"double edge for labelling {wl.labelling}"


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
def test_2d_done_sampled_larger(spec):
    for wl in enumerate_wall_labelings([2, 2], seed=17, limit=120):
        stg = SpecCubicalBlowupGraph(
            labelling=wl.labelling, num_thresholds=[2, 2], spec=spec, level=3
        )
        assert not _double_edges(stg), f"double edge for labelling {wl.labelling}"


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
@pytest.mark.parametrize("name", ["N2_a", "N2_b", "toggle"])
def test_2d_done_dsgrn_derived(spec, name):
    """The same theorem over every DSGRN parameter of the 2D example networks."""
    total = networks.parameter_graph(name).size()
    for index in range(total):
        stg = SpecCubicalBlowupGraph(
            networks.parameter(name, index), spec=spec, level=3
        )
        assert not _double_edges(stg), f"{name}[{index}]"


# ---------------------------------------------------------------------------
# Pointwise nonemptiness: prop:F1-welldefined, cor:F2-well-defined,
# prop:F3-welldefined  (CombinatorialDynamics.tex:139, ~880, ~1110)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
@pytest.mark.parametrize("level", [1, 2, 3, 4])
@pytest.mark.parametrize("name", ["toggle", "repressilator", "cycle3"])
def test_maps_are_pointwise_nonempty(spec, level, name):
    """`F_i(xi) != empty` for every cell.

    Every map in the monograph is required to be nonempty pointwise; this is
    what makes `SCC(F_i)` and the induced grading well defined.
    """
    total = networks.parameter_graph(name).size()
    indices = range(total) if total <= 30 else range(0, total, max(1, total // 30))
    for index in indices:
        stg = SpecCubicalBlowupGraph(
            networks.parameter(name, index), spec=spec, level=level
        )
        for cell in stg.blowup_complex(stg.dim):
            if stg.blowup_complex.rightfringe(cell):
                continue
            assert stg.digraph.adjacencies(cell), (
                f"{name}[{index}] level={level} spec={spec.name}: "
                f"F_{level}(cell {cell}) is empty"
            )


# ---------------------------------------------------------------------------
# Refinement order (def:refinement, CombinatorialDynamics.tex:16)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
@pytest.mark.parametrize("name", ["toggle", "repressilator", "cycle3"])
def test_f2_refines_f1(spec, name):
    """`F_2 = F_1 cap F_{2.1}` is a refinement, so its edges are a subset."""
    total = networks.parameter_graph(name).size()
    indices = range(total) if total <= 20 else range(0, total, max(1, total // 20))
    for index in indices:
        p = networks.parameter(name, index)
        f1 = set(SpecCubicalBlowupGraph(p, spec=spec, level=1).digraph.edges())
        f2 = set(SpecCubicalBlowupGraph(p, spec=spec, level=2).digraph.edges())
        assert f2 <= f1, f"{name}[{index}]: F_2 is not a refinement of F_1"


@pytest.mark.parametrize("name", ["toggle", "repressilator", "cycle3"])
def test_f3_is_f2_refinement_plus_unstable_cells(name):
    """`F_3 = (F_2 cap F_{3.1}) u U`, so `F_3 \\ F_2` consists of U-edges only.

    `defn:Rule3` adds cells, so `F_3` is *not* a pointwise submap of `F_2`;
    every extra edge must come from `U(xi)` (`defn:Rule3.2`).
    """
    spec = PAPER
    total = networks.parameter_graph(name).size()
    indices = range(total) if total <= 20 else range(0, total, max(1, total // 20))
    for index in indices:
        p = networks.parameter(name, index)
        stg2 = SpecCubicalBlowupGraph(p, spec=spec, level=2)
        stg3 = SpecCubicalBlowupGraph(p, spec=spec, level=3)
        extra = set(stg3.digraph.edges()) - set(stg2.digraph.edges())
        allowed = set()
        for cc_cell in stg3.cubical_complex:
            if stg3.cubical_complex.rightfringe(cc_cell):
                continue
            targets = stg3.unstable_targets(cc_cell)
            if not targets:
                continue
            source = stg3.cubical2blowup(cc_cell)
            allowed.update((source, stg3.cubical2blowup(t)) for t in targets)
        assert extra <= allowed, (
            f"{name}[{index}]: F_3 gained edges not accounted for by U(xi)"
        )


# ---------------------------------------------------------------------------
# prop:CH=0 -- LatticeStructures13.tex:1918
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.name)
@pytest.mark.parametrize("name", ["toggle", "repressilator", "cycle3"])
@pytest.mark.parametrize("level", [2, 4])
def test_common_gradient_direction_implies_trivial_index(spec, name, level):
    """`p in GRC(F_i) => CH_*(p) = 0`.

    The manuscript records that the only offered proof of `prop:CH=0` uses a
    later ramp theorem with extra hypotheses, so this is worth checking
    directly.  `MorseGraph.nontrivial_scc` relies on it: it drops SCCs with a
    common gradient direction, which would discard a nonzero index if the
    proposition failed.
    """
    total = networks.parameter_graph(name).size()
    indices = range(total) if total <= 20 else range(0, total, max(1, total // 20))
    for index in indices:
        result = conley_morse_graph(
            networks.parameter(name, index), spec=spec, level=level
        )
        stg, gc = result.stg, result.graded_complex
        fringe_grade = gc.value(stg.complex().size() - 1)

        cells_by_grade: dict[int, list[int]] = {}
        for cell in stg.digraph.vertices():
            value = gc.value(cell)
            if value == fringe_grade:
                continue
            cells_by_grade.setdefault(value, []).append(cell)

        counts = result.connection_matrix.count()
        for grade, cells in cells_by_grade.items():
            common = None
            for cell in cells:
                dirs = set(stg.gradient_directions(stg.blowup2cubical(cell)))
                common = dirs if common is None else (common & dirs)
                if not common:
                    break
            if not common:
                continue
            index_vector = counts.get(grade)
            assert not (index_vector and any(index_vector)), (
                f"{name}[{index}] level={level} spec={spec.name}: grade {grade} "
                f"has common gradient directions {sorted(common)} but "
                f"Conley index {index_vector} -- prop:CH=0 would be false"
            )


# ---------------------------------------------------------------------------
# Why the conditions of the intersection can be evaluated selectively
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["repressilator", "cycle3", "cycle4", "N3_B"])
def test_condition_2_1_never_contradicts_F1(name):
    """Condition 2.1 can only fire where `F_1` left a double edge.

    `cor:F2-well-defined` argues that every GO-pair produces an `F_1` double
    edge, so a pair `F_1` has already oriented is not in `D(Phi)` and the
    decision wall cannot speak about it.  `F_2 = F_1 cap F_{2.1}` may therefore
    skip the decision wall wherever `F_1` returned a verdict -- which is what
    makes the intersection affordable in six dimensions.

    `strict_intersection=True` disables the skip and raises if the two ever
    disagree; this test drives that path.
    """
    from rookfields.blowup import SpecCubicalBlowupGraph

    total = networks.parameter_graph(name).size()
    indices = range(total) if total <= 20 else range(0, total, max(1, total // 20))
    for index in indices:
        p = networks.parameter(name, index)
        strict = SpecCubicalBlowupGraph(
            p, spec=PAPER, level=2, strict_intersection=True
        )
        fast = SpecCubicalBlowupGraph(p, spec=PAPER, level=2)
        assert set(strict.digraph.edges()) == set(fast.digraph.edges())
