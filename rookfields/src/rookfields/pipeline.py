"""Wall labelling -> F_i -> SCC -> Morse graph -> Conley complex.

Mirrors ``DSGRN_utils.ComputeMorseGraph.ConleyMorseGraph`` (ComputeMorseGraph.py:8)
but is spec-aware and returns the ``scc_dag`` and ``connection_matrix`` that the
original computes and discards -- the reason
``SaveDatabaseJSON_CubicalBlowup.py:311-312`` recomputes them.

The step sequence is the one stated in ``prelude.tex:190-296``:

    omega  ->  Phi  ->  F_i  ->  SCC(F_i), pi  ->  pi_b on X_b  ->  MG, Delta
"""

from __future__ import annotations

import dataclasses
import time
from collections import Counter

import pychomp
from DSGRN_utils.MorseGraph import MorseGraph

from .blowup import SpecCubicalBlowupGraph
from .spec import PAPER, Spec


@dataclasses.dataclass
class MorseResult:
    """Everything the pipeline produces for one wall labelling."""

    morse_graph: object
    stg: SpecCubicalBlowupGraph
    graded_complex: object
    scc_dag: object
    connection_matrix: object
    spec: Spec
    level: int
    seconds: float

    # -- convenience accessors used by the example checks and the audit ----

    @property
    def nodes(self) -> list[int]:
        """Morse-node indices, in the Morse graph's own sorted order."""
        return sorted(self.morse_graph.vertex_label(v)[0] for v in self.morse_graph.vertices())

    @property
    def conley_indices(self) -> dict[int, tuple]:
        """``{morse node: Betti vector}``, Z_2 coefficients."""
        return {
            self.morse_graph.vertex_label(v)[0]: tuple(self.morse_graph.vertex_label(v)[2])
            for v in self.morse_graph.vertices()
        }

    @property
    def top_cell_counts(self) -> dict[int, int]:
        """``{morse node: number of top-dimensional cells}``."""
        return {
            self.morse_graph.vertex_label(v)[0]: self.morse_graph.vertex_label(v)[1]
            for v in self.morse_graph.vertices()
        }

    @property
    def edges(self) -> set[tuple[int, int]]:
        """Morse-graph edges as index pairs (already a transitive reduction)."""
        label = {v: self.morse_graph.vertex_label(v)[0] for v in self.morse_graph.vertices()}
        return {
            (label[v], label[u])
            for v in self.morse_graph.vertices()
            for u in self.morse_graph.adjacencies(v)
        }

    @property
    def stable_nodes(self) -> list[int]:
        """Morse nodes with no outgoing edge (attractors)."""
        return sorted(
            self.morse_graph.vertex_label(v)[0]
            for v in self.morse_graph.vertices()
            if not self.morse_graph.adjacencies(v)
        )

    @property
    def trivial_index_nodes(self) -> list[int]:
        """Morse nodes whose Conley index vanishes in every degree."""
        return sorted(n for n, ci in self.conley_indices.items() if not any(ci))

    def _cells_by_node(self) -> dict[int, set[int]]:
        """Invert the grading once; the STG has O(10^5) vertices in six dimensions."""
        cached = getattr(self, "_cells_by_node_cache", None)
        if cached is None:
            grade_of_node = {
                self.morse_graph.vertex_label(v)[0]: v
                for v in self.morse_graph.vertices()
            }
            node_of_grade = {v: n for n, v in grade_of_node.items()}
            cached = {n: set() for n in grade_of_node}
            for cell in self.stg.digraph.vertices():
                node = node_of_grade.get(self.graded_complex.value(cell))
                if node is not None:
                    cached[node].add(cell)
            self._cells_by_node_cache = cached
        return cached

    def morse_set_cells(self, node: int) -> set[int]:
        """Blowup cells carrying Morse node ``node``."""
        cells = self._cells_by_node()
        if node not in cells:
            raise KeyError(f"no Morse node {node}")
        return cells[node]

    def double_edge_pairs(self, node: int) -> set[frozenset]:
        """Unordered pairs of cells in one Morse set joined by arrows both ways."""
        cells = self.morse_set_cells(node)
        pairs = set()
        for c in cells:
            for d in self.stg.digraph.adjacencies(c):
                if d != c and d in cells and c in self.stg.digraph.adjacencies(d):
                    pairs.add(frozenset((c, d)))
        return pairs

    def equilibrium_cells_in(self, node: int) -> set[int]:
        """Cells of Morse set ``node`` that are equilibrium cells (``def:eqcell``)."""
        return {
            cell
            for cell in self.morse_set_cells(node)
            if self.stg.equilibrium_cell(self.stg.blowup2cubical(cell))
        }

    @property
    def diagnostics(self) -> Counter:
        return self.stg.diagnostics


def conley_morse_graph(
    parameter=None,
    *,
    labelling=None,
    num_thresholds=None,
    spec: Spec = PAPER,
    level: int = 4,
) -> MorseResult:
    """Run the full pipeline for one DSGRN parameter or one raw wall labelling."""
    if parameter is None and (labelling is None or num_thresholds is None):
        raise ValueError("provide either parameter, or labelling and num_thresholds")
    if parameter is not None and (labelling is not None or num_thresholds is not None):
        raise ValueError("provide parameter, or labelling and num_thresholds, not both")

    start = time.perf_counter()
    stg = SpecCubicalBlowupGraph(
        parameter=parameter,
        labelling=labelling,
        num_thresholds=num_thresholds,
        spec=spec,
        level=level,
    )
    scc_dag, graded_complex = pychomp.FlowGradedComplex(stg.complex(), stg.adjacencies())
    connection_matrix = pychomp.ConnectionMatrix(graded_complex)
    morse_graph = MorseGraph(
        stg, scc_dag, graded_complex, connection_matrix, prune_grad=spec.prune_grad
    )
    elapsed = time.perf_counter() - start

    return MorseResult(
        morse_graph=morse_graph,
        stg=stg,
        graded_complex=graded_complex,
        scc_dag=scc_dag,
        connection_matrix=connection_matrix,
        spec=spec,
        level=level,
        seconds=elapsed,
    )


def morse_node_correspondence(a: MorseResult, b: MorseResult) -> dict | None:
    """Match the Morse nodes of two results by their cell sets.

    Returns ``{node in a: node in b}`` when every Morse set of ``a`` occurs
    verbatim in ``b``, otherwise ``None``.

    Morse-node *numbering* is not canonical: ``MorseGraph`` sorts by rank and
    breaks ties by the internal SCC grading value (MorseGraph.py:79-86), so two
    runs that produce the same decomposition can number it differently.  Any
    comparison between specs has to go through the cell sets.
    """
    cells_a = {n: frozenset(a.morse_set_cells(n)) for n in a.nodes}
    cells_b = {n: frozenset(b.morse_set_cells(n)) for n in b.nodes}
    if len(set(cells_b.values())) != len(cells_b):
        return None  # two nodes with identical cell sets: cannot match uniquely
    inverse = {v: k for k, v in cells_b.items()}
    mapping = {}
    for node, cells in cells_a.items():
        target = inverse.get(cells)
        if target is None:
            return None
        mapping[node] = target
    return mapping


def same_morse_decomposition(a: MorseResult, b: MorseResult) -> bool:
    """Do two results give the same Morse decomposition, up to renumbering?

    Compares the Morse sets as collections of cells, the Conley index carried by
    each, and the reachability order -- but not the node numbers.
    """
    mapping = morse_node_correspondence(a, b)
    if mapping is None:
        return False
    if any(a.conley_indices[n] != b.conley_indices[m] for n, m in mapping.items()):
        return False
    return {(mapping[u], mapping[v]) for (u, v) in a.edges} == b.edges


#: Verdicts returned by :func:`compare_results`, from weakest to strongest
#: difference.
IDENTICAL = "identical"
SAME_INVARIANTS = "same invariants, different cells"
DIFFERENT = "different"


def compare_results(a: MorseResult, b: MorseResult) -> str:
    """Classify how two results differ.

    ``IDENTICAL``
        the same Morse decomposition, possibly with the nodes renumbered;
    ``SAME_INVARIANTS``
        the same number of Morse nodes and the same multiset of Conley indices
        -- so every homological quantity the monograph reports is unchanged --
        but at least one Morse set is a different collection of cells;
    ``DIFFERENT``
        the node count or the Conley indices themselves differ.

    Morse-node numbering is not canonical (``MorseGraph.py:79-86`` breaks rank
    ties by the internal SCC grading value), so nothing here compares node
    labels.
    """
    if same_morse_decomposition(a, b):
        return IDENTICAL
    same_count = len(a.nodes) == len(b.nodes)
    same_indices = sorted(a.conley_indices.values()) == sorted(b.conley_indices.values())
    if same_count and same_indices:
        return SAME_INVARIANTS
    return DIFFERENT
