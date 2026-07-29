"""Algorithm specifications for the combinatorial multivalued maps F_0 .. F_4.

The implementation in ``DSGRN_utils.CubicalBlowupGraph`` and the definitions in
``paper/Rook_Field_Paper_v2`` have drifted apart.  Rather than pick a winner, we
parameterise every point of disagreement so both semantics can be computed and
diffed.  ``LEGACY`` reproduces ``DSGRN_utils`` exactly; ``PAPER`` implements the
live manuscript definitions.

Each field documents the manuscript definition it tracks and the DSGRN_utils
line that currently implements the other reading.
"""

from __future__ import annotations

import dataclasses
from typing import Literal

# --- D1: codomain of the regulation map o_xi ------------------------------
#
# Live `def:active_regulation` (RookFields6.tex:738-758) reads
#
#     Let xi in X and let n, m in J_i(xi).  We say that n actively regulates m
#     at xi if there exist n-adjacent top cells mu, mu' in Top_X(xi) whose
#     common n-wall contains xi, such that Phi_m(xi, mu) != Phi_m(xi, mu').
#
# so the regulation map is typed  o_xi : Act(xi) -> J_i(xi).  The superseded
# `def:active_regulation-old` (RookFields6.tex:726-736) let k range over all of
# {1..N}.  CubicalBlowupGraph.active_regulation_map (:298-309) ranges k over
# `range(self.dim)` -- i.e. the old definition.
RegulationCodomain = Literal["all", "inessential"]

# --- D1b: domain of the regulation map ------------------------------------
#
# DSGRN_utils restricts n to *interior* inessential directions
# (`0 < coords[n] < limits[n]`, CubicalBlowupGraph.py:302).  The manuscript
# states no such restriction, but strong dissipativity makes boundary walls
# absorbing.  Exposed so the restriction can be audited; both presets keep it.
RegulationDomain = Literal["interior", "all"]

# --- D2: back walls that leave the complex --------------------------------
#
# `defn:back-walls` (CombinatorialDynamics.tex:356-366) shifts by
#
#     v_hat = sum_{n in J_i(xi')} (1 + r_n)/2 * 0^{(n)},   r_n in R_n(xi)
#
# which is not a cell of X when a direction inessential at xi' sits at an outer
# endpoint.  `issues/source-defects.md` records the safe condition as
# 1 <= v'_j <= K(j) for every such j, and leaves the algorithm's response to the
# authors.  CubicalBlowupGraph.decision_wall (:427-434) applies the shift with
# no bounds check and hands negative coordinates to `cell_index`, which wraps to
# an unrelated cell.
#
#   "wrap"            -- current DSGRN_utils behaviour (no check)
#   "discard_shifted" -- drop the pair when the *chosen* shift leaves X
#   "discard_strict"  -- drop the pair unless 1 <= v'_j <= K(j) for all
#                        j in J_i(xi'), i.e. unless every back wall in the set
#                        is a cell of X
#   "error"           -- raise, for use in tests
BackWallOutOfRange = Literal["wrap", "discard_shifted", "discard_strict", "error"]

# --- D2b: which back wall to read -----------------------------------------
#
# `Back(xi, xi')` is a *set*, indexed by the choices r_n in R_n(xi).
# `prop:back_wall_well_defined` (CombinatorialDynamics.tex:762) asserts the
# n_o-rook value does not depend on the choice.  DSGRN_utils picks the first
# top cell of the star (`face_top_star[0]`, :426).  "all_must_agree" evaluates
# every choice and reports disagreement -- an empirical test of that
# proposition.
BackWallChoice = Literal["first", "all_must_agree"]

# --- D3: how F_3 / F_4 are composed ---------------------------------------
#
# `defn:Rule3` (CombinatorialDynamics.tex:1102-1108) is
#
#     F_3(xi) = ( F_2(xi) cap F_{3.1}(xi) ) cup U(xi)
#
# an intersection of independently-computed refinements, followed by an
# unconditional union.  CubicalBlowupGraph.compute_multivalued_map (:621-682)
# is a first-match cascade: Condition 3.1 and U(xi) are only reached when both
# F_1 and the decision wall returned 0 for that pair.
F3Composition = Literal["cascade", "intersect_union"]

# --- D4: the two excluded classes in defn:partially_opaque -----------------
#
# `defn:partially_opaque` (CombinatorialDynamics.tex:999-1013) requires a
# semi-opaque cell to lie in X \ (X^(N) u bdy(X)) *and* have o_xi a bijection of
# Act(xi).  CubicalBlowupGraph.semi_opaque_cell (:470-475) checks only the
# bijection; a top cell has an empty regulation map and so passes vacuously, and
# a cell on bdy(X) via an *essential* direction is not excluded either.
#
# Both presets default to "implicit" (today's behaviour) because the effect has
# not been measured; the audit reports on it.
SemiOpaqueGuard = Literal["implicit", "explicit"]


@dataclasses.dataclass(frozen=True)
class Spec:
    """A choice of reading for every known paper/code divergence."""

    name: str
    regulation_codomain: RegulationCodomain = "all"
    regulation_domain: RegulationDomain = "interior"
    back_wall_out_of_range: BackWallOutOfRange = "wrap"
    back_wall_choice: BackWallChoice = "first"
    f3_composition: F3Composition = "cascade"
    semi_opaque_guard: SemiOpaqueGuard = "implicit"
    self_edges: bool = True
    prune_grad: bool = True
    #: Homology coefficients.  The connection-matrix enumeration of
    #: LatticeStructures13.tex:1968-1976 is only tractable over Z_2.
    coefficients: int = 2

    def replace(self, **changes) -> "Spec":
        changes.setdefault("name", f"{self.name}+" + ",".join(sorted(changes)))
        return dataclasses.replace(self, **changes)

    def __str__(self) -> str:  # pragma: no cover - display only
        return self.name


#: Exactly today's ``DSGRN_utils`` behaviour.  Every previously published figure
#: must reproduce bit-for-bit under this spec.
LEGACY = Spec(name="legacy")

#: The live manuscript definitions.
PAPER = Spec(
    name="paper",
    regulation_codomain="inessential",
    back_wall_out_of_range="discard_strict",
    f3_composition="intersect_union",
    semi_opaque_guard="explicit",
)


#: One-toggle-at-a-time variants, so the audit can attribute each difference to
#: a single divergence rather than to the bundle.
ISOLATED = {
    "D1-regulation-codomain": LEGACY.replace(
        name="D1-regulation-codomain", regulation_codomain="inessential"
    ),
    "D2-back-wall-strict": LEGACY.replace(
        name="D2-back-wall-strict", back_wall_out_of_range="discard_strict"
    ),
    "D2-back-wall-shifted": LEGACY.replace(
        name="D2-back-wall-shifted", back_wall_out_of_range="discard_shifted"
    ),
    "D3-f3-composition": LEGACY.replace(
        name="D3-f3-composition", f3_composition="intersect_union"
    ),
    "D4-semi-opaque-guard": LEGACY.replace(
        name="D4-semi-opaque-guard", semi_opaque_guard="explicit"
    ),
    "D2b-back-wall-all": LEGACY.replace(
        name="D2b-back-wall-all", back_wall_choice="all_must_agree"
    ),
    "D1b-regulation-domain": LEGACY.replace(
        name="D1b-regulation-domain", regulation_domain="all"
    ),
}


def by_name(name: str) -> Spec:
    """Resolve a spec name from the CLI."""
    if name in ("legacy", "LEGACY"):
        return LEGACY
    if name in ("paper", "PAPER"):
        return PAPER
    if name in ISOLATED:
        return ISOLATED[name]
    raise KeyError(
        f"unknown spec {name!r}; known: legacy, paper, " + ", ".join(sorted(ISOLATED))
    )


ALL_SPECS = {"legacy": LEGACY, "paper": PAPER, **ISOLATED}
