"""Spec-parameterised construction of the multivalued maps F_0 .. F_4.

``SpecCubicalBlowupGraph`` subclasses ``DSGRN_utils.CubicalBlowupGraph`` and
overrides exactly the methods that the paper/code divergences touch.  Everything
else -- the rook field, wall labels, flow directions, lap numbers, unstable
cells, and the blowup/cubical index translation -- is reused unchanged, because
those already match the manuscript.

Manuscript sources (``paper/Rook_Field_Paper_v2``):

* ``def:rookfield``            RookFields6.tex:630-668
* ``def:active_regulation``    RookFields6.tex:738-758
* ``def:regulationmap``        RookFields6.tex:863-875
* ``defn:opaque``              RookFields6.tex:1087
* ``def:eqcell``               RookFields6.tex:1222
* ``def:Rule0``                CombinatorialDynamics.tex:25-33
* ``def:Rule1.1/1.2/Rule1``    CombinatorialDynamics.tex:65-105
* ``defn:GO-pair``             CombinatorialDynamics.tex:221-239
* ``defn:back-walls``          CombinatorialDynamics.tex:356-366
* ``defn:indecisive``          CombinatorialDynamics.tex:367-391
* ``defn:F2`` (Condition 2.1)  CombinatorialDynamics.tex:806-822
* ``def:Rule2``                CombinatorialDynamics.tex:834-841
* ``defn:partially_opaque``    CombinatorialDynamics.tex:999-1013
* ``defn:Rule3.1``             CombinatorialDynamics.tex:1064-1072
* ``defn:Rule3.2`` (U(xi))     CombinatorialDynamics.tex:1073-1101
* ``defn:Rule3``               CombinatorialDynamics.tex:1102-1108
* ``defn:F4`` (Condition 4.1)  CombinatorialDynamics.tex:1184-1199
* ``def:Rule4``                CombinatorialDynamics.tex:1200-1209
"""

from __future__ import annotations

import itertools
from collections import Counter

from DSGRN_utils.CubicalBlowupGraph import CubicalBlowupGraph

from .spec import PAPER, Spec


class BackWallOutOfRangeError(ValueError):
    """Raised under ``back_wall_out_of_range='error'``."""


class SpecCubicalBlowupGraph(CubicalBlowupGraph):
    """A ``CubicalBlowupGraph`` whose algorithm follows an explicit ``Spec``.

    Diagnostics are accumulated in ``self.diagnostics`` (a ``Counter``) during
    construction so the audit can attribute differences without re-running the
    computation.
    """

    def __init__(
        self,
        parameter=None,
        *,
        spec: Spec = PAPER,
        labelling=None,
        num_thresholds=None,
        level: int = 4,
        strict_intersection: bool = False,
    ):
        self.spec = spec
        self.diagnostics: Counter = Counter()
        #: (face, coface) pairs where two conditions demanded opposite
        #: orientations, so the definition leaves no edge between them.
        self.conflicting_pairs: list = []
        #: (face, coface) pairs whose back-wall set is not contained in X.
        self.out_of_range_pairs: list = []
        #: cells where prop:back_wall_well_defined fails empirically.
        self.back_wall_disagreements: list = []
        super().__init__(
            parameter=parameter,
            labelling=labelling,
            num_thresholds=num_thresholds,
            self_edges=spec.self_edges,
            level=level,
            # Every leaf method the specs touch is overridden below, so the base
            # class flag only matters for code paths we do not override.
            legacy=(spec.regulation_codomain == "all"),
            strict_intersection=strict_intersection,
        )

    # ------------------------------------------------------------------
    # D1 -- the regulation map o_xi
    # ------------------------------------------------------------------

    def compute_active_regulation_map(self, cc_cell):
        """The regulation map ``o_xi`` at ``cc_cell``.

        Overrides the base computation; the base class memoises it.

        Live ``def:active_regulation`` types this ``Act(xi) -> J_i(xi)``;
        ``DSGRN_utils`` lets the target range over all N directions, which is
        the retired ``def:active_regulation-old``.
        """
        coords = self.cubical_complex.coordinates(cc_cell)
        inessential = self.inessential_directions(cc_cell)

        if self.spec.regulation_domain == "interior":
            domain = [n for n in inessential if 0 < coords[n] < self.limits[n]]
        else:
            domain = list(inessential)

        if self.spec.regulation_codomain == "inessential":
            targets = inessential
        else:
            targets = range(self.dim)

        reg_map = {}
        for n in domain:
            hits = [m for m in targets if self.active_regulation(cc_cell, n, m)]
            if not hits:
                continue
            if len(hits) > 1:
                # `def:regulationmap` presents o_xi as a map; a multivalued
                # entry would contradict that.  Record rather than silently
                # take the last, which is what the DSGRN_utils dict
                # comprehension (CubicalBlowupGraph.py:308) does.
                self.diagnostics["regulation_map_multivalued"] += 1
            reg_map[n] = hits[-1]
        return reg_map

    # ------------------------------------------------------------------
    # boundary / semi-opaque helpers
    # ------------------------------------------------------------------

    def in_boundary(self, cc_cell) -> bool:
        """``xi in bdy(X)``: some inessential direction sits at an outer end."""
        coords = self.cubical_complex.coordinates(cc_cell)
        return any(
            coords[n] == 0 or coords[n] == self.limits[n]
            for n in self.inessential_directions(cc_cell)
        )

    def compute_semi_opaque_cell(self, cc_cell) -> bool:
        """``defn:partially_opaque``.

        The manuscript requires ``xi in X \\ (X^(N) u bdy(X))`` with ``o_xi`` a
        bijection of ``Act(xi)``.  ``DSGRN_utils`` checks only the bijection.
        Under ``semi_opaque_guard='explicit'`` the two excluded classes are
        enforced as well.
        """
        if self.spec.semi_opaque_guard == "explicit":
            if self.cubical_complex.cell_dim(cc_cell) == self.dim:
                return False
            if self.in_boundary(cc_cell):
                return False
        reg_map = self.active_regulation_map(cc_cell)
        return set(reg_map.keys()) == set(reg_map.values())

    def nontrivial_cycles(self, cc_cell):
        """Cycles of length >= 2 in the cycle decomposition of ``o_xi``."""
        if not self.semi_opaque_cell(cc_cell):
            return []
        reg_map = self.active_regulation_map(cc_cell)
        return [c for c in self.cycle_decomposition(reg_map) if len(c) > 1]

    def unstable_targets(self, cc_cell):
        """``U(xi)`` of ``defn:Rule3.2``, computed cellwise and unconditionally."""
        cycles = self.nontrivial_cycles(cc_cell)
        if not cycles:
            return []
        return self.unstable_cells(cc_cell, cycles)

    # ------------------------------------------------------------------
    # D2 -- back walls
    # ------------------------------------------------------------------

    def rook_value_set(self, cc_cell, n):
        """``R_n(xi) = { Phi_n(xi, mu) : mu in Top_X(xi) }`` (eq:Rn)."""
        return {
            self.rook_field_component(cc_cell, top, n)
            for top in self.top_star(cc_cell)
        }

    def back_wall_range_ok(self, cc_coface) -> bool:
        """``1 <= v'_j <= K(j)`` for every ``j in J_i(xi')``.

        The condition under which every member of ``Back(xi, xi')`` is a cell of
        ``X`` (``issues/source-defects.md``, critical path B3).
        """
        coords = self.cubical_complex.coordinates(cc_coface)
        return all(
            1 <= coords[j] <= self.num_thresholds[j]
            for j in self.inessential_directions(cc_coface)
        )

    def back_walls(self, cc_face, cc_coface, n_opaque):
        """Every ``(top_cell, n_o, side)`` in ``Back(xi, xi')``.

        ``defn:back-walls`` indexes the set by independent choices
        ``r_n in R_n(xi)`` for ``n in J_i(xi')``; ``DSGRN_utils`` uses the single
        coherent choice realised by ``top_star(xi)[0]``.
        """
        face_coords = self.cubical_complex.coordinates(cc_face)
        coface_coords = self.cubical_complex.coordinates(cc_coface)
        side = -1 if face_coords[n_opaque] == coface_coords[n_opaque] else 1
        directions = self.inessential_directions(cc_coface)

        if self.spec.back_wall_choice == "first":
            top0 = self.top_star(cc_face)[0]
            choices = [
                tuple(self.rook_field_component(cc_face, top0, n) for n in directions)
            ]
        else:
            per_direction = [sorted(self.rook_value_set(cc_face, n)) for n in directions]
            choices = list(itertools.product(*per_direction)) if per_direction else [()]

        walls = []
        for choice in choices:
            new_coords = list(coface_coords)
            for n, r_n in zip(directions, choice):
                # J_i(xi') subset J_i(xi), so Phi_n(xi, .) is a wall label and
                # r_n in {-1, +1}; 0 would make (1 + r_n)/2 non-integral.
                assert r_n in (-1, 1), f"rook value {r_n} on an inessential direction"
                if r_n == 1:
                    new_coords[n] -= 1
            walls.append((new_coords, side))
        return walls

    def decision_wall(self, cc_face, cc_coface):
        """One decision wall for ``(xi, xi')``, or ``()``.

        Kept for compatibility with ``DSGRN_utils``; the spec-aware logic lives
        in :meth:`decision_wall_direction`.
        """
        walls = self._decision_walls(cc_face, cc_coface)
        if not walls:
            return tuple()
        n_opaque, side, coords_list = walls
        return (
            self.cubical_complex.cell_index(coords_list[0], self.top_shape),
            n_opaque,
            side,
        )

    def _decision_walls(self, cc_face, cc_coface):
        """``(n_opaque, side, [coords, ...])`` for an indecisive-drift pair.

        Returns ``()`` when there is no GO-pair, when the drift is not
        indecisive, or when the spec rejects an out-of-range back-wall set.
        """
        go_pair = self.gradient_opaque_pair(cc_face, cc_coface)
        if not go_pair:
            return ()
        _n_grad, n_opaque = go_pair

        # defn:indecisive (ii): every other opaque inessential target has only
        # itself in its regulation-map fibre.
        act_map = self.active_regulation_map(cc_face)
        face_inessential = set(self.inessential_directions(cc_face))
        opaque_others = (
            set(self.opaque_directions(cc_face)) & face_inessential
        ) - {n_opaque}
        if any(act_map[n] in opaque_others and act_map[n] != n for n in act_map):
            return ()

        policy = self.spec.back_wall_out_of_range
        if policy != "wrap" and not self.back_wall_range_ok(cc_coface):
            self.diagnostics["back_wall_out_of_range_strict"] += 1
            self.out_of_range_pairs.append((cc_face, cc_coface))
            if policy == "error":
                raise BackWallOutOfRangeError(
                    f"Back({cc_face}, {cc_coface}) is not contained in X"
                )
            if policy == "discard_strict":
                return ()

        walls = self.back_walls(cc_face, cc_coface, n_opaque)
        side = walls[0][1]
        coords_list = [coords for coords, _ in walls]

        if any(c < 0 or c > self.limits[j] for coords in coords_list for j, c in enumerate(coords)):
            self.diagnostics["back_wall_shift_out_of_range"] += 1
            if policy in ("discard_shifted", "discard_strict"):
                return ()
            if policy == "error":
                raise BackWallOutOfRangeError(
                    f"back-wall shift for ({cc_face}, {cc_coface}) leaves X"
                )
            # "wrap": fall through, reproducing the DSGRN_utils behaviour of
            # handing negative coordinates to cell_index.

        return (n_opaque, side, coords_list)

    def decision_wall_direction(self, cc_cell1, cc_cell2):
        """Condition 2.1 (``defn:F2``) or Condition 4.1 (``defn:F4``).

        Returns ``-1`` to remove ``face -> coface``, ``+1`` to remove
        ``coface -> face``, ``0`` for no constraint.
        """
        if cc_cell1 < cc_cell2:
            cc_face, cc_coface, face_sign = cc_cell1, cc_cell2, 1
        else:
            cc_face, cc_coface, face_sign = cc_cell2, cc_cell1, -1

        walls = self._decision_walls(cc_face, cc_coface)
        if not walls:
            return 0
        n_opaque, side, coords_list = walls

        labels = set()
        for coords in coords_list:
            top = self.cubical_complex.cell_index(coords, self.top_shape)
            labels.add(self.wall_label(top, n_opaque, side))
        if len(labels) > 1:
            # prop:back_wall_well_defined asserts this cannot happen.
            self.diagnostics["back_wall_disagreement"] += 1
            self.back_wall_disagreements.append((cc_face, cc_coface))
            return 0
        label = labels.pop()

        if label == side:
            # xi_hat in E^-(xi_hat'):  xi' not in F_{2.1}(xi)
            return -face_sign

        # xi_hat in E^+(xi_hat').  Condition 2.1 additionally requires
        #     n_o not in Act(xi)   or   (n_o in Act(xi) and o_xi(n_o) = n_o).
        # Condition 4.1 drops that requirement.
        if self.level < 4:
            act_map = self.active_regulation_map(cc_face)
            if n_opaque in act_map and act_map[n_opaque] != n_opaque:
                return 0
        return face_sign

    # ------------------------------------------------------------------
    # Condition 3.1
    # ------------------------------------------------------------------

    def condition_3_1_direction(self, cc_cell1, cc_cell2):
        """``defn:Rule3.1``: remove ``xi -> xi'`` when ``Ex(xi, xi') subset S_sigma``."""
        if cc_cell1 < cc_cell2:
            cc_face, cc_coface, face_sign = cc_cell1, cc_cell2, 1
        else:
            cc_face, cc_coface, face_sign = cc_cell2, cc_cell1, -1

        cycles = self.nontrivial_cycles(cc_face)
        if not cycles:
            return 0
        ext = set(self.extension_directions(cc_face, cc_coface))
        if any(ext <= set(cycle) for cycle in cycles):
            return -face_sign
        return 0

    # ------------------------------------------------------------------
    # combinatorial data consumed by the analytic hypotheses
    # ------------------------------------------------------------------

    def indecisive_drift_pairs(self, codim: tuple[int, int] | None = None):
        """``D(Phi)``: pairs ``(xi, xi')`` exhibiting indecisive drift.

        Yields ``(cc_face, cc_coface, n_grad, n_opaque)``.  ``codim`` optionally
        restricts to a pair of cell dimensions, e.g. ``(N-2, N-1)`` for the
        pairs H2 quantifies over (``defn:H2``, RampSystemsv4.tex:572).
        """
        cc = self.cubical_complex
        for cc_coface in cc:
            if cc.rightfringe(cc_coface):
                continue
            coface_dim = cc.cell_dim(cc_coface)
            if codim is not None and coface_dim != codim[1]:
                continue
            for cc_face in cc.boundary({cc_coface}):
                if cc.rightfringe(cc_face):
                    continue
                if codim is not None and cc.cell_dim(cc_face) != codim[0]:
                    continue
                go_pair = self.gradient_opaque_pair(cc_face, cc_coface)
                if not go_pair:
                    continue
                if not self._decision_walls(cc_face, cc_coface):
                    continue
                n_grad, n_opaque = go_pair
                yield (cc_face, cc_coface, n_grad, n_opaque)

    def has_edge(self, cc_source, cc_target) -> bool:
        """``target in F(source)``, in cubical-cell terms."""
        return self.cubical2blowup(cc_target) in set(
            self.digraph.adjacencies(self.cubical2blowup(cc_source))
        )

    def three_cycle_vertices(self):
        """Vertices whose regulation map is a 3-cycle (``defn:H3``)."""
        cc = self.cubical_complex
        out = []
        for cell in cc(0):
            if cc.rightfringe(cell):
                continue
            cycles = self.nontrivial_cycles(cell)
            if len(cycles) == 1 and len(cycles[0]) == 3:
                out.append((cell, cycles[0]))
        return out

    # ------------------------------------------------------------------
    # graph assembly
    # ------------------------------------------------------------------

    def compute_multivalued_map(self):
        if self.spec.f3_composition == "cascade":
            return super().compute_multivalued_map()
        return self._compute_intersect_union()

    def _compute_intersect_union(self):
        """``F_i`` as an intersection of refinements, then a union with ``U``.

        ``def:Rule2``   F_2 = F_1 cap F_{2.1}
        ``defn:Rule3``  F_3 = ( F_2 cap F_{3.1} ) u U
        ``def:Rule4``   F_4 = ( F_1 cap F_{3.1} cap F_{4.1} ) u U
        """
        if self.level == 0:
            return self.trivial_multivalued_map()

        use_decision_wall = self.level > 1
        use_cycles = self.level > 2

        for cell1 in self.blowup_complex(self.dim):
            cc_cell1 = self.blowup2cubical(cell1)

            # Condition 1.1: the self-arrow survives only at equilibrium cells.
            if self.spec.self_edges and self.equilibrium_cell(cc_cell1):
                self.digraph.add_edge(cell1, cell1)

            for cell2 in self.parallel_neighbors(cell1):
                fringe1 = self.blowup_complex.rightfringe(cell1)
                fringe2 = self.blowup_complex.rightfringe(cell2)
                if fringe1:
                    self.digraph.add_edge(cell1, cell2)
                if fringe2:
                    self.digraph.add_edge(cell2, cell1)
                if fringe1 or fringe2:
                    continue

                cc_cell2 = self.blowup2cubical(cell2)

                flow = self.flow_direction(cc_cell1, cc_cell2)
                verdicts = [flow]
                # See CubicalBlowupGraph.compute_paper_multivalued_map: every
                # GO-pair leaves an F_1 double edge, so Condition 2.1 cannot
                # fire on a pair F_1 has already oriented.
                if use_decision_wall and (flow == 0 or self.strict_intersection):
                    decision = self.decision_wall_direction(cc_cell1, cc_cell2)
                    if flow != 0 and decision not in (0, flow):
                        raise AssertionError(
                            "Condition 2.1 contradicted F_1 at "
                            f"({cc_cell1}, {cc_cell2}): {flow} vs {decision}"
                        )
                    verdicts.append(decision)
                if use_cycles:
                    verdicts.append(self.condition_3_1_direction(cc_cell1, cc_cell2))

                # A verdict of -1 removes cell1 -> cell2; +1 removes
                # cell2 -> cell1.  Refinements intersect, so a removal by any
                # single condition is final.
                remove_forward = any(v == -1 for v in verdicts)
                remove_backward = any(v == 1 for v in verdicts)

                if remove_forward and remove_backward:
                    self.diagnostics["conflicting_pair"] += 1
                    self.conflicting_pairs.append((cc_cell1, cc_cell2))
                    continue
                if not remove_forward:
                    self.digraph.add_edge(cell1, cell2)
                if not remove_backward:
                    self.digraph.add_edge(cell2, cell1)

        if use_cycles:
            self._add_unstable_cells()

    def _add_unstable_cells(self):
        """Union with ``U(xi)``, applied to every cell rather than per pair."""
        for cc_cell in self.cubical_complex:
            if self.cubical_complex.rightfringe(cc_cell):
                continue
            targets = self.unstable_targets(cc_cell)
            if not targets:
                continue
            source = self.cubical2blowup(cc_cell)
            for cc_target in targets:
                self.digraph.add_edge(source, self.cubical2blowup(cc_target))
                self.diagnostics["unstable_edge"] += 1
