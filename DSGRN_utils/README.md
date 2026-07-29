# DSGRN_utils

DSGRN utilities.

## This copy diverges from upstream

Upstream is <https://github.com/marciogameiro/DSGRN_utils>. This copy is vendored
into the RookFields repository and **implements the definitions of
`Rook_Field_Paper_v2`**, which upstream does not. Passing `legacy=True` to
`CubicalBlowupGraph`, `ConleyMorseGraph`, or `save_morse_graph_database_json`
restores the previous behaviour, so figures produced before the alignment still
reproduce exactly.

| | file | change |
|---|---|---|
| C1 | `CubicalBlowupGraph.py` `active_regulation_map` | targets of `o_xi` restricted to `J_i(xi)`, per the live `def:active_regulation`; the previous code used the retired `def:active_regulation-old`, which ranged over all `N` directions |
| C2 | `CubicalBlowupGraph.py` `decision_wall` | back-wall pairs outside `1 <= v'_j <= K(j)` are discarded instead of handing negative coordinates to `cell_index`, which wrapped to an unrelated cell |
| C3 | `WallLabelling.py` | new `global_bound=` argument gives the `GB_n` sentinel of `defn:ramp-wall-labeling`; the previous `theta + 10 h` surrogate depends on `h` and loses strong dissipativity for small `h` |
| C4 | `CubicalBlowupGraph.py` `semi_opaque_cell` | excludes top cells and `bdy(X)`, as `defn:partially_opaque` requires |
| C5 | `IsomorphismQuery.py` | rewritten around a Morse-graph signature; it previously crashed against current HEAD, and used `continue` where `break` was meant |
| — | `CubicalBlowupGraph.py` `compute_multivalued_map` | assembles `F_i` as `defn:Rule3` states — an intersection of refinements plus an unconditional union with `U(xi)` — rather than as a first-match cascade |

Two performance changes keep the faithful assembly affordable: the cell-indexed
quantities are memoised, and Condition 2.1 is evaluated only where `F_1` gave no
verdict. The second is not an approximation — `cor:F2-well-defined` argues every
GO-pair leaves an `F_1` double edge, so a pair `F_1` has already oriented is not
in `D(Phi)`. `strict_intersection=True` disables the skip and raises on any
disagreement.

Measured effect of each change, and the evidence behind it, is in
`../rookfields/reports/FINDINGS.md`. The differences are exercised by
`../rookfields/tests/test_legacy_equivalence.py`, which pins both readings.

The upstream git history for this directory was moved to
`../../.dsgrn_utils-upstream.git` when it was vendored; merges from upstream are
now manual.
