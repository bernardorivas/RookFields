"""Every computational example in the monograph, with its claimed output.

Each entry records where the claim is made and what the text asserts, so the
reproduction pass can report claimed-vs-computed rather than merely "it ran".

The manuscript knowledge base is explicit that these numbers "are source
assertions from the stated software runs [...] not treated as independently
recomputed" (``project-knowledge/manuscript/maps/examples.md``).  This module
is that recomputation.

Runtimes quoted in the text were measured by the authors on a laptop using one
core (``examples4.tex:32``); they are recorded for context, not asserted.
"""

from __future__ import annotations

import dataclasses
from typing import Any, Callable


@dataclasses.dataclass(frozen=True)
class Example:
    """One computational example.

    ``network``/``index`` select a DSGRN parameter; ``ramp_system`` instead
    names an entry of :data:`rookfields.networks.RAMP_SYSTEMS` for the two
    directly specified systems that bypass DSGRN.
    """

    label: str
    source: str
    level: int
    network: str | None = None
    index: int | None = None
    ramp_system: str | None = None
    #: Human-readable summary of what the text claims.
    claim: str = ""
    #: Runtime in seconds quoted in the text, if any.
    claimed_seconds: float | None = None
    #: Checks run against a :class:`~rookfields.pipeline.MorseResult`.
    #: Each returns ``(name, claimed, computed, ok)``.
    checks: tuple[Callable[[Any], tuple[str, Any, Any, bool]], ...] = ()
    #: Marked ``True`` for the two 6-node EMT runs, which dominate wall clock.
    expensive: bool = False
    notes: str = ""

    @property
    def is_ramp(self) -> bool:
        return self.ramp_system is not None


# ---------------------------------------------------------------------------
# check factories
# ---------------------------------------------------------------------------


def n_morse_nodes(expected: int):
    def check(r):
        got = len(r.nodes)
        return ("number of Morse nodes", expected, got, got == expected)

    return check


def n_trivial_index_nodes(expected: int):
    def check(r):
        got = len(r.trivial_index_nodes)
        return ("trivial-index Morse nodes", expected, got, got == expected)

    return check


def trivial_index_nodes_are(expected: list[int]):
    def check(r):
        got = r.trivial_index_nodes
        return ("trivial-index node set", sorted(expected), got, got == sorted(expected))

    return check


def conley_index_of(node: int, expected: tuple):
    def check(r):
        got = r.conley_indices.get(node)
        return (f"Conley index CH_*({node})", expected, got, got == expected)

    return check


def morse_order(expected: set[tuple[int, int]]):
    def check(r):
        got = r.edges
        return ("Morse graph edges", sorted(expected), sorted(got), got == set(expected))

    return check


def morse_set_size(node: int, expected: int):
    def check(r):
        got = len(r.morse_set_cells(node))
        return (f"cells in Morse set {node}", expected, got, got == expected)

    return check


def _unique_trivial_node(r) -> int | None:
    nodes = r.trivial_index_nodes
    return nodes[0] if len(nodes) == 1 else None


def trivial_component_size(expected: int):
    """Size of the unique trivial-index component, located by its index.

    Morse-node *numbering* depends on the sort in ``MorseGraph`` and so is not
    stable across specs; the trivial-index component itself is.
    """

    def check(r):
        node = _unique_trivial_node(r)
        if node is None:
            return (
                "cells in the trivial-index component",
                expected,
                f"not unique: {r.trivial_index_nodes}",
                False,
            )
        got = len(r.morse_set_cells(node))
        return (
            f"cells in the trivial-index component (node {node})",
            expected,
            got,
            got == expected,
        )

    return check


def trivial_component_double_edges(expected: int):
    def check(r):
        node = _unique_trivial_node(r)
        if node is None:
            return (
                "double-edge pairs in the trivial-index component",
                expected,
                f"not unique: {r.trivial_index_nodes}",
                False,
            )
        got = len(r.double_edge_pairs(node))
        return (
            f"double-edge pairs in the trivial-index component (node {node})",
            expected,
            got,
            got == expected,
        )

    return check


def _unique_multicell(r):
    nodes = [n for n in r.nodes if len(r.morse_set_cells(n)) > 1]
    return nodes[0] if len(nodes) == 1 else None


def multicell_component_size(expected: int):
    """Size of the unique component with more than one cell.

    Morse-node numbering is not stable across specs, so components are located
    by an intrinsic property rather than by index.
    """

    def check(r):
        node = _unique_multicell(r)
        if node is None:
            sizes = {n: len(r.morse_set_cells(n)) for n in r.nodes}
            return ("cells in the multi-cell component", expected, f"not unique: {sizes}", False)
        got = len(r.morse_set_cells(node))
        return (
            f"cells in the multi-cell component (node {node})",
            expected,
            got,
            got == expected,
        )

    return check


def multicell_component_double_edges(expected: int):
    def check(r):
        node = _unique_multicell(r)
        if node is None:
            return (
                "double-edge pairs in the multi-cell component",
                expected,
                "not unique",
                False,
            )
        got = len(r.double_edge_pairs(node))
        return (
            f"double-edge pairs in the multi-cell component (node {node})",
            expected,
            got,
            got == expected,
        )

    return check


def periodic_component_size(expected: int, dim: int):
    """Size of the component carrying the stable-periodic-orbit index.

    Located by its Conley index rather than by node number, which is not stable
    across specs.
    """
    signature = tuple([1, 1] + [0] * (dim - 1))

    def check(r):
        nodes = [n for n, ci in r.conley_indices.items() if ci == signature]
        if len(nodes) != 1:
            return (
                "cells in the periodic-orbit component",
                expected,
                f"not unique: {nodes}",
                False,
            )
        got = len(r.morse_set_cells(nodes[0]))
        return (
            f"cells in the periodic-orbit component (node {nodes[0]})",
            expected,
            got,
            got == expected,
        )

    return check


def periodic_component_double_edges(expected: int, dim: int):
    signature = tuple([1, 1] + [0] * (dim - 1))

    def check(r):
        nodes = [n for n, ci in r.conley_indices.items() if ci == signature]
        if len(nodes) != 1:
            return (
                "double-edge pairs in the periodic-orbit component",
                expected,
                f"not unique: {nodes}",
                False,
            )
        got = len(r.double_edge_pairs(nodes[0]))
        return (
            f"double-edge pairs in the periodic-orbit component (node {nodes[0]})",
            expected,
            got,
            got == expected,
        )

    return check


def double_edge_pairs_in(node: int, expected: int):
    def check(r):
        got = len(r.double_edge_pairs(node))
        return (f"double-edge pairs in Morse set {node}", expected, got, got == expected)

    return check


def equilibrium_cells_in(node: int, expected: int):
    def check(r):
        got = len(r.equilibrium_cells_in(node))
        return (f"equilibrium cells in Morse set {node}", expected, got, got == expected)

    return check


def n_equilibrium_cells(expected: int):
    def check(r):
        got = len(
            [
                c
                for c in r.stg.cubical_complex
                if r.stg.equilibrium_cell(c)
            ]
        )
        return ("equilibrium cells in X", expected, got, got == expected)

    return check


def index_set(expected: tuple[int, ...]):
    def check(r):
        got = tuple(k + 1 for k in r.stg.num_thresholds)
        return ("index set I (top index per coordinate)", expected, got, got == expected)

    return check


def all_indices_nonzero():
    def check(r):
        got = r.trivial_index_nodes
        return ("every Morse node has nonzero index", [], got, not got)

    return check


#: Periodic-orbit index signature: CH_0 = CH_1 = Z_2, zero elsewhere.
def periodic_orbit_index(node: int, dim: int):
    expected = tuple([1, 1] + [0] * (dim - 1))

    def check(r):
        got = r.conley_indices.get(node)
        return (
            f"CH_*({node}) is the stable-periodic-orbit signature",
            expected,
            got,
            got == expected,
        )

    return check


# ---------------------------------------------------------------------------
# the catalog
# ---------------------------------------------------------------------------

EXAMPLES: list[Example] = [
    # -- Introduction ----------------------------------------------------
    Example(
        label="ex:saddle_saddle_3D_intro",
        source="introduction6.tex:402",
        network="N3_A",
        index=52_718_681_992,
        level=3,
        claim="25 Morse nodes, every Conley index nonzero (hyperbolic fixed points)",
        checks=(n_morse_nodes(25), all_indices_nonzero()),
        notes="Same parameter node as ex:saddle_saddle_3D_example.",
    ),
    Example(
        label="ex:introPeriodic",
        source="introduction6.tex:498",
        ramp_system="intro_periodic",
        level=3,
        claim="Morse node 0 has CH_0 = CH_1 = Z_2 (stable periodic orbit index)",
        checks=(periodic_orbit_index(0, 3),),
    ),
    # -- Chapter 16, global dynamics and bifurcations ---------------------
    Example(
        label="ex:saddlesaddlebif",
        source="examples4.tex:39",
        network="N2_a",
        index=752,
        level=3,
        claim="two admissible connection matrices; 1600 nodes in about 8 s",
        claimed_seconds=8.0,
        checks=(),
        notes="The 8 s is for all 1600 parameter nodes, not this one.",
    ),
    Example(
        label="ex:saddle_saddle_3D_example",
        source="examples4.tex:263",
        network="N3_A",
        index=52_718_681_992,
        level=3,
        claim="25 equilibrium-type Morse nodes; 4096 connection matrices; 0.19 s",
        claimed_seconds=0.19,
        checks=(n_morse_nodes(25), all_indices_nonzero()),
    ),
    # -- periodic orbits --------------------------------------------------
    Example(
        label="ex:periodicOrbit3",
        source="examples4.tex:307",
        ramp_system="intro_periodic",
        level=3,
        claim="I = {0..8}x{0..6}x{0..3}; M(0) has the periodic-orbit index; 0.32 s",
        claimed_seconds=0.32,
        checks=(index_set((8, 6, 3)), periodic_orbit_index(0, 3)),
    ),
    Example(
        label="ex:periodicOrbit3_F3",
        source="examples4.tex:435",
        network="N3_B",
        index=2_472_287,
        level=3,
        claim="Morse node 1 has 19 cells, 6 double-edge pairs, and an equilibrium cell",
        checks=(
            multicell_component_size(19),
            multicell_component_double_edges(6),
        ),
        notes=(
            "The source calls the fixed-point-plus-periodic-orbit reading a "
            "conjecture. Under the `paper` spec the component is instead a clean "
            "12-cell cycle with no double edges and the stable-periodic-orbit "
            "index -- i.e. F_3 already yields the resolution the text attributes "
            "to F_4. Attributable to D1 alone."
        ),
    ),
    # -- semiconjugacies --------------------------------------------------
    Example(
        label="ex:mccord",
        source="examples4.tex:489",
        network="N2_b",
        index=974,
        level=3,
        claim="3 Morse nodes with order 0<2, 1<2; CH_0(0)=CH_0(1)=Z_2, CH_1(2)=Z_2",
        checks=(
            n_morse_nodes(3),
            morse_order({(2, 0), (2, 1)}),
            conley_index_of(0, (1, 0, 0)),
            conley_index_of(1, (1, 0, 0)),
            conley_index_of(2, (0, 1, 0)),
        ),
        notes="R(974) is the running wall labelling of fig:wall_labeling(A).",
    ),
    Example(
        label="ex:3d_example_1",
        source="examples4.tex:538",
        network="N3_B",
        index=2_472_286,
        level=3,
        claim="13-node Morse representation, unique connection matrix; 0.14 s",
        claimed_seconds=0.14,
        checks=(n_morse_nodes(13), all_indices_nonzero()),
    ),
    Example(
        label="EMT-6D-semiconjugacy",
        source="examples4.tex:663",
        network="N6_EMT",
        index=2_684_686_006,
        level=2,
        claim="Morse graph and connection matrix in about 2.2 minutes",
        claimed_seconds=132.0,
        checks=(),
        expensive=True,
    ),
    # -- Morse graphs that are not Morse representations -------------------
    Example(
        label="ex:trivial_index_2D",
        source="examples4.tex:766",
        network="N2_a",
        index=47,
        level=3,
        claim="4 Morse nodes: three with nonzero index, one trivial",
        checks=(n_morse_nodes(4), n_trivial_index_nodes(1)),
    ),
    Example(
        label="ex:trivial_index_3D_F3",
        source="examples4.tex:786",
        network="N3_A",
        index=65_571_607_721,
        level=3,
        claim="three one-cell equilibrium nodes plus a 61-cell, 12-double-edge trivial node",
        checks=(
            n_morse_nodes(4),
            n_trivial_index_nodes(1),
            trivial_component_size(61),
            trivial_component_double_edges(12),
        ),
        notes=(
            "fig:trivial_index_3D_F2 is misnamed; the example uses F_3. "
            "The published 61/12 reproduce under the `paper` spec but not under "
            "`legacy`, which gives 75 cells and 16 double-edge pairs."
        ),
    ),
    Example(
        label="ex:trivial_indices",
        source="examples4.tex:801",
        network="N3_A",
        index=52_717_613_010,
        level=3,
        claim="8 connection matrices; trivial index at p in {10,11,13,18,27} with "
        "recurrent components of 44, 42, 42, 14, 14 cells; 0.2 s",
        claimed_seconds=0.2,
        checks=(
            trivial_index_nodes_are([10, 11, 13, 18, 27]),
            morse_set_size(10, 44),
            morse_set_size(11, 42),
            morse_set_size(13, 42),
            morse_set_size(18, 14),
            morse_set_size(27, 14),
        ),
    ),
    Example(
        label="EMT-6D-trivial-indices",
        source="examples4.tex:820",
        network="N6_EMT",
        index=1_739_757_491_101,
        level=2,
        claim="eight trivial-index nodes, p in {15..18, 21..24}; about 2.5 minutes",
        claimed_seconds=150.0,
        checks=(
            n_morse_nodes(25),
            n_trivial_index_nodes(8),
        ),
        expensive=True,
        notes=(
            "The count reproduces exactly (8 of 25). The displayed node *labels* "
            "do not: legacy gives {15,17,18,19,21,22,23,24} and paper gives "
            "{15,16,18,19,21,22,23,24}. Morse-node numbering in `MorseGraph` "
            "(MorseGraph.py:79-86) sorts by rank and breaks ties by the internal "
            "SCC grading value, so it is not canonical; only the count and the "
            "index data are."
        ),
    ),
    # -- systems not generated by DSGRN -----------------------------------
    Example(
        label="ex:ramp_van_der_pol",
        source="examples4.tex:841",
        ramp_system="van_der_pol",
        level=3,
        claim="F_3 and Morse graph for a multi-ramp planar system; 0.03 s",
        claimed_seconds=0.03,
        checks=(index_set((6, 6)),),
        notes="The example does not verify h in H_3; see the hypotheses report.",
    ),
    # -- Chapter 16 section 'Dynamics Computed Using F_4' ------------------
    Example(
        label="ex:periodicOrbit3_F4",
        source="examples4.tex:917",
        network="N3_B",
        index=2_472_287,
        level=4,
        claim="I = {0..4}x{0..4}x{0..3}; M(1) has the periodic-orbit index; "
        "12-cell cycle with no double edges; 0.12 s",
        claimed_seconds=0.12,
        checks=(
            index_set((4, 4, 3)),
            periodic_component_size(12, 3),
            periodic_component_double_edges(0, 3),
        ),
        notes=(
            "No ODE realisation is claimed for any F_4 example (examples4.tex:915). "
            "The displayed index set is off in its first factor: node 1 of N3_B has "
            "two out-edges, so K(1) = 2 and I = {0..3} x {0..4} x {0..3}. Everything "
            "else in the example reproduces exactly."
        ),
    ),
    Example(
        label="ex:periodicOrbit5",
        source="examples4.tex:979",
        network="N5",
        index=5_103_162_287,
        level=4,
        claim="7 equilibrium cells (Morse nodes 1-7); node 0 is a 12-cell cycle "
        "with the periodic-orbit index; 9.15 s",
        claimed_seconds=9.15,
        checks=(
            n_morse_nodes(8),
            periodic_component_size(12, 5),
            n_equilibrium_cells(7),
        ),
    ),
    Example(
        label="ex:missing_periodic_3D_F4",
        source="examples4.tex:1063",
        network="N3_B",
        index=2_023_186,
        level=4,
        claim="Morse node 3 is a 29-cell SCC with 4 unresolved double-edge pairs",
        checks=(
            multicell_component_size(29),
            multicell_component_double_edges(4),
        ),
        notes=(
            "The 29-cell component carries a fixed-point index (1,0,0,0), not a "
            "trivial one; the source conjectures its invariant part holds both a "
            "fixed point and a periodic orbit."
        ),
    ),
]

BY_LABEL = {e.label: e for e in EXAMPLES}
