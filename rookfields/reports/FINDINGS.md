# Findings: `Rook_Field_Paper_v2` against `DSGRN_utils`

Everything below is reproducible with `pytest code/rookfields/tests` and the
report generators listed at the end. No file under `paper/` was modified.

**Status: the manuscript's definitions have been adopted as canonical.**
`code/DSGRN_utils` now implements them by default; passing `legacy=True` to
`CubicalBlowupGraph` / `ConleyMorseGraph` restores the superseded behaviour so
figures published before the alignment still reproduce. The changes are C1–C5
below.

Two specs are used throughout. `legacy` reproduces the pre-alignment
`DSGRN_utils` bit-for-bit; `paper` implements the live manuscript definitions and
is the default. Both are pinned against the two modes of `DSGRN_utils` in
`test_legacy_equivalence.py`, across five levels and four network families, so
`rookfields` and plain `DSGRN_utils` cannot drift apart.

---

## 1. The code implements a superseded definition of active regulation

**Where.** `CubicalBlowupGraph.active_regulation_map` (`:298-309`) builds
`{n: k}` with `k` ranging over `range(self.dim)`. The live
`def:active_regulation` (`RookFields6.tex:738-758`) requires `n, m in J_i(xi)`
and types the map `o_xi : Act(xi) -> J_i(xi)`. The retired
`def:active_regulation-old` (`:726-736`, the `<<` branch) is the one with `k` in
`{1..N}`.

**Effect.** Entries whose target is an *essential* direction reach the
`level < 4` guard (`:461-467`) and `semi_opaque_cell` (`:470`), changing `F_2`
and `F_3`. Measured (parameters whose result changes):

| family | dim | sampled | STG | Morse-node count | Conley indices |
|---|---:|---:|---:|---:|---:|
| toggle, N2_a, N2_b | 2 | 409 | 0 | 0 | 0 |
| repressilator | 3 | 27 | 7 | 1 | 1 |
| cycle3 | 3 | 120 | 38 | 1 | 1 |
| N3_B | 3 | 80 | 50 | 2 | 2 |
| cycle4 | 4 | 60 | 33 | 4 | 4 |

Every two-dimensional example in the monograph is unaffected.

**It already broke a published example.** `ex:trivial_index_3D_F3`
(`examples4.tex:786`) reports a 61-cell trivial-index component with 12
double-edge pairs. Under `paper` that reproduces exactly; under `legacy` the
component has **75 cells and 16 double-edge pairs**. The published number came
from the manuscript's definition, and the current code no longer reproduces it.

**It also changes the `F_3` vs `F_4` story.** `ex:periodicOrbit3_F3`
(`examples4.tex:435`, N3_B node 2,472,287) reports a 19-cell component with 6
double-edge pairs and calls the fixed-point-plus-periodic-orbit reading a
conjecture; `ex:periodicOrbit3_F4` then presents `F_4` as resolving it into a
clean 12-cell cycle with the stable-periodic-orbit index. Under the corrected
definition **`F_3` already gives that 12-cell cycle**, with no double edges and
`CH_0 = CH_1 = Z_2`. Attributable to this divergence alone: no other toggle
changes the result.

---

## 2. Back walls leave the complex and the wrong cell is read

**Where.** `defn:back-walls` (`CombinatorialDynamics.tex:356-366`) shifts by
`v_hat = sum_{n in J_i(xi')} (1 + r_n)/2 * 0^{(n)}`. `issues/source-defects.md`
already records that this leaves `X` at an outer endpoint, with safe condition
`1 <= v'_j <= K(j)`, and leaves the algorithm's response to the authors.
`decision_wall` (`:427-434`) performs the shift with **no bounds check** and
passes negative coordinates to `cell_index`, which wraps to an unrelated cell.

**Measured**, on decision walls that are actually constructed:

| family | decision walls | negative coordinate | wrapped to a different cell |
|---|---:|---:|---:|
| repressilator | 240 | 54 (22%) | 54 (100%) |
| 3-node (120 params) | 2585 | 396 (15%) | 396 (100%) |
| 4-node (120 params) | 11010 | 4052 (37%) | 4052 (100%) |

**Effect.** Discarding those pairs changes the STG on 19/27, 112/120, 80/80 and
54/60 parameters in dimensions 3 and 4 — but changes **no** Morse-node count and
**no** Conley index in any family sampled. So it perturbs the graph without
disturbing the homological conclusions in these examples. Two-dimensional
examples are unaffected.

**Still an authorial decision**: whether such pairs are excluded from
`D(Phi)`, or handled another way. The `back_wall_out_of_range` toggle offers
`wrap` / `discard_shifted` / `discard_strict` / `error`.

---

## 3. `prop:back_wall_well_defined` holds; two other divergences are inert

Three readings were parameterised and turned out to change nothing anywhere:

- **`D2b`** — evaluating *every* choice `r_n in R_n(xi)` rather than the single
  coherent one `DSGRN_utils` uses gives the same wall label in every case
  sampled. Positive evidence for `prop:back_wall_well_defined`
  (`CombinatorialDynamics.tex:762`).
- **`D3`** — `defn:Rule3` is an intersection of refinements plus an
  unconditional union with `U(xi)`; `compute_multivalued_map` (`:621-682`) is a
  first-match cascade. Implemented both ways: identical on every parameter, and
  a direct search for suppressed `U(xi)` edges found **0** across 496
  parameters and thousands of cells with `U(xi) != 0`. The cascade is not the
  stated definition but is empirically equivalent here.
- **`D1b`** — restricting `o_xi`'s domain to *interior* inessential directions,
  which the text does not do, is harmless.

The regulation map was also checked for multivaluedness (which would contradict
`def:regulationmap` presenting `o_xi` as a map): **no violation** in ~90,000
entries across the DSGRN families and both ramp systems.

---

## 4. `defn:partially_opaque`'s two excluded classes are not enforced

`semi_opaque_cell` (`:470-475`) checks only that `o_xi` is a bijection.
`defn:partially_opaque` (`CombinatorialDynamics.tex:999-1013`) also requires
`xi` to lie in `X \ (X^(N) u bdy(X))`. Enforcing that changes the STG on 16/60
(`cycle4`) and 20/80 (`N3_B`) parameters at levels 3 and 4, with no change to
Morse nodes or Conley indices. Invisible in dimensions 2 and 3 for the smaller
families.

---

## 5. `WallLabelling.py` uses an ad-hoc outer sentinel that breaks dissipativity

**Where.** `defn:ramp-wall-labeling` (`RampSystemsv4.tex:375-390`) sets the
outermost sentinel threshold to `theta_{m_{K(n)+1}} = GB_n`.
`ramp_system_wall_labelling` (`WallLabelling.py:44`) uses
`theta_{m_K} + 10 h_{m_K}` instead.

**Effect.** That surrogate depends on `h`, so:

- it contradicts `prop:wall-labeling-const` (`:490`: the labelling is constant
  in `h` on `H_0`) — shrinking `h` from 0.15 to 0.0198 flips 2 of 36 top-cell
  labels for the van der Pol system;
- it **breaks strong dissipativity**, the standing assumption from
  `CombinatorialDynamics.tex:4` onward, once `h <= 0.05`.

| width | `theta + 10h` sentinel | `GB_n` sentinel |
|---|:-:|:-:|
| published (0.11–0.15) | dissipative | dissipative |
| 0.05, 0.02, 0.005 | **not dissipative** | dissipative |

With the `GB_n` sentinel the labelling is byte-identical at every width tested,
exactly as `prop:wall-labeling-const` claims. The bug is invisible at the
published widths, which is presumably why it survived — but small `h` is the
regime `H_1` forces.

---

## 6. The published ramp parameters violate `H_1`, hence `H_2` and `H_3`

Both directly specified ramp systems —
`tab:parameters_van_der_pol_ramp_system` (`examples4.tex:857`) and
`tab:parameters_ramp_system_intro_periodic` (`introduction6.tex:525`) —
satisfy `H_0` and `Lambda(R)` but **fail `H_1`**.

Concretely, for van der Pol at box `(1,1)`: `E_1(D_v) = 6.1`, `gamma_1 = 2`, so
the focal value is `3.05`, which lies in the forbidden half-window
`(theta_{m_1}, theta_{m_1} + h_{m_1}] = (3, 3.15]`. Eight such violations for
van der Pol, eighteen for the three-dimensional system.

| uniform half-width | H_0 | Lambda(R) | H_1 | H_2 | H_3 |
|---|:-:|:-:|:-:|:-:|:-:|
| published | ok | ok | **fails** | **fails** | **fails** |
| 0.02 | ok | ok | **fails** | **fails** | **fails** |
| 0.01 | ok | ok | ok | **fails** | ok |
| 0.005 | ok | ok | ok | ok | ok |

**Consequences.** `H_2` and `H_3` are subsets of `H_1`, so the hypotheses of
`thm:R1ABlattice` and `thm:R3ABlattice` do not hold at the published values.
`ex:introPeriodic` draws an ODE conclusion from `thm:dynamics` at exactly those
values.

**Repair.** `h = 0.005` uniformly satisfies every condition, and the Morse graph
is unchanged — same Conley indices, same Morse-graph edges. So the computed
combinatorics stand; only the stated widths need changing. The largest uniform
width satisfying `H_0 + H_1 + H_2` is `0.00517`; `prop:smallh`'s constructive
`H_1` bound is `0.02`, which is not small enough for `H_2` — consistent with the
recorded concern that the `H_2` constant `c` is defined from quantities that
already depend on `h`.

Note also that `prop:smallh` states nonemptiness "for each `i in {0,1,3}`" while
its proof also treats `i = 2`.

---

## 7. `thm:R1ABlattice` verified numerically, and its hypothesis is load-bearing

For every oriented wall of `X_b`, the inward component of the ramp field is
sampled on the embedded face. At `h = 0.005` (admissible):

| system | level | oriented walls | misaligned | worst inward component |
|---|---:|---:|---:|---:|
| van der Pol | 1 | 298 | **0** | +0.03 |
| van der Pol | 2 and 3 | 312 | 11 | -1.71 |
| intro periodic | 1 | 4099 | **0** | +0.03 |
| intro periodic | 2 and 3 | 4180 | 66 | -1.71 |

At the published widths `F_1` is *not* aligned (18 and 166 walls fail), so the
width condition is doing real work.

**All** `F_2`/`F_3` misalignments sit on walls that `F_1` left as double edges
and the higher rule newly oriented — precisely the faces the R2 GO-manifold and
R3 cycle-surface constructions are meant to replace. The checker therefore
localises exactly which faces need geometric modification.

---

## 8. R2: the GO-manifolds build, and every checkable claim holds

`defn:GO-manifold` (`GlobalDynR2v9.tex:723-780`) was implemented and run on all
62 codimension-two indecisive-drift pairs of the two ramp systems (14 for van
der Pol, 48 for the 3D system; 22 external, 40 internal).

- the `delta`-nullcline base lies on `F_{n_o} = 0` to `1.4e-13`;
- `Q^delta(xi, xi')` is degenerate in `n_o` for every pair, as
  `GlobalDynR2v9.tex:196` asserts, so the base has dimension `N-2` and the sheet
  is codimension one;
- every externally pruned pair reaches the directed `n_g`-wall in finite time
  (`lem:GO-wall-hitting`, whose proof uses `H_2`);
- `prop:GO-crossing-sign` holds at machine precision: the deviation
  `F - F_eps` is exactly `r eps x_{n_g} e^{(n_g)}` in the `(n_o, n_g)` block
  (residuals `~1e-15`);
- `lem:GO-product` holds exactly: transverse drift `0` along every trajectory.

**One practical point.** `defn:admissible-GO-perturbation` (`:1161-1180`)
requires `r_{n_g} F_{eps,n_g} > 0` on the control region. That is equivalent to
`eps < r F_{n_g}(x) / x_{n_g}` throughout, and the bound is tight: past it the
perturbation reverses the `n_g`-motion and the sweep never reaches the directed
wall. The admissible range is narrow — `eps < 0.017` for some pairs of the 3D
system versus `eps < 0.138` for others.

---

## 9. R3: two recorded defects reproduce

**The `H_3` bound does not force a complex stable pair.** `eq:F3-bounds` is
equivalent to `|P| > D(gamma)` for the signed slope product
`P = prod Delta E_n / (8 prod h)`, with
`D(gamma) = -g1 g2 g3 + (g1+g2+g3)(g1 g2 + g1 g3 + g2 g3)`. The witness recorded
in `issues/source-defects.md`, `gamma = (1, 2, 100)` and `P = 31000`, gives
`D = 30906 < 31000` and eigenvalues `-96.57, -21.36, +14.93` — all real. A
systematic sweep finds many more, e.g. `gamma = (0.5, 0.5, 20)` for every
`P` in `(420.25, 840.5]`.

**The purported transverse vector is a generating tangent.** A general
tangent/normal detector for swept surfaces is implemented and confirms the
distinction on a controlled example.

**Neither published ramp system has a nontrivial regulation cycle**, so
Conditions 3.1 and 3.2 are vacuous and `F_3 = F_2` for both. In particular
`ex:ramp_van_der_pol` is presented as an `F_3` computation, but the cycle rule
never fires at those parameters. Cycle cells do occur for DSGRN parameters —
N3_B node 2,472,287 has one 3-cycle; N3_A node 52,718,681,992 has twelve
2-cycles and one 3-cycle.

---

## 10. Comparing Morse graphs: node numbering is not canonical

`MorseGraph` sorts nodes by rank and breaks ties by the internal SCC grading
value (`MorseGraph.py:79-86`), so two runs producing the *same* decomposition
can number it differently. Comparing labelled graphs therefore reports
differences that are not there. `rookfields.pipeline.compare_results` classifies
a pair of runs as

- **identical** — the same Morse decomposition, up to renumbering;
- **same invariants, different cells** — the same number of Morse nodes and the
  same multiset of Conley indices, so every homological quantity the monograph
  reports is unchanged, but at least one Morse set is a different collection of
  cells;
- **different** — the node count or the Conley indices themselves differ.

Under that classification, comparing the two readings across all 17 examples:

| example | dim | verdict |
|---|---:|---|
| `ex:periodicOrbit3_F3` | 3 | **different** — 13 nodes vs 14; the 19-cell component becomes a 12-cell periodic-orbit cycle (§1) |
| `ex:trivial_index_3D_F3` | 3 | same invariants, different cells — the trivial component is 75 cells under legacy, 61 under paper (§1) |
| `ex:3d_example_1` | 3 | same invariants, different cells — one component is 33 cells vs 29; the 13-node graph and its unique connection matrix are unchanged |
| `EMT-6D-trivial-indices` | 6 | same invariants, different cells — 25 nodes and 8 trivial-index nodes either way, but every large component shrinks (2436 → 2394, 440 → 340, 415 → 336, 390 → 301, 300 → 262, 77 → 59) |
| `EMT-6D-semiconjugacy` | 6 | identical — two nodes carrying the same index merely swap numbers |
| the other twelve | | identical |

So the divergences reach dimension six, but outside `ex:periodicOrbit3_F3` they
move cells between components without changing any Conley index or node count.
The one place a *printed* number changes is `ex:trivial_index_3D_F3`, where the
text quotes the cell count.

## 11. Example reproduction: 30 of 34 runs match exactly

All 17 catalogued computational examples were recomputed under both specs, including the
two six-dimensional EMT runs (91 s and 148 s here against the reported 2.2 and
2.5 minutes). Four mismatches, all explained:

| example | spec | issue |
|---|---|---|
| `ex:trivial_index_3D_F3` | legacy | 75/16 instead of the published 61/12; `paper` reproduces exactly (see §1) |
| `ex:periodicOrbit3_F3` | paper | resolves into the 12-cell cycle the text attributes to `F_4` (see §1) |
| `ex:periodicOrbit3_F4` | both | index set misprinted (below) |

**`ex:periodicOrbit3_F4` index set.** `examples4.tex:920-923` displays
`I = {0,1,2,3,4} x {0,1,2,3,4} x {0,1,2,3}`. Node 1 of that network has two
out-edges, so `K(1) = 2` and the first factor is `{0,1,2,3}`. Everything else in
the example — the 12-cell cycle, the absence of double edges, the
periodic-orbit index — reproduces exactly.

Two further labelling notes: the eight trivial-index nodes of the second EMT
example reproduce in *count* (8 of 25) but not in the displayed labels
`{15..18, 21..24}`, because Morse-node numbering in `MorseGraph`
(`MorseGraph.py:79-86`) breaks rank ties by the internal SCC grading value and
is not canonical. And `fig:trivial_index_3D_F2` is named `F2` for an `F_3`
example, as already recorded.

All four large networks reproduce their stated parameter-graph sizes exactly
(87,280,405,632 / 3,600,000 / 13,608,000,000 / 4,429,771,960,320), and both
ramp systems reproduce their stated index sets.

---

## 12. Positive results

- `thm:2dDone` survives an exhaustive test: **1,660 abstract wall labelings**
  on two-dimensional complexes, generated directly from `def:wall_labeling` plus
  strong dissipativity rather than from DSGRN, produce **no** `F_3` double edge
  under either spec — plus every DSGRN parameter of all three 2-node networks.
- `prop:CH=0` holds on every sampled SCC: no strongly connected component has
  both a common gradient direction and a nonzero Conley index. `MorseGraph`'s
  gradient pruning depends on this.
- `F_i(xi) != empty` pointwise at every level, under both specs.
- `F_2` refines `F_1`; `F_3 \ F_2` consists only of `U(xi)` edges.
- The corrected Janus carrier passes both regressions: adjacent coarse top cells
  have disjoint fine top slices, and the 1D subdivision chain has the right
  endpoint `8(Kbar(n)+1)`.
- The rectangular geometrization fills `prod [0, GB_n]` exactly, its charts
  restrict correctly to faces, and `g(b(xi))` agrees with `I_n(xi)` on every
  cell not touching a sentinel index.

---

## 13. A sentinel mismatch between two chapters

`defn:ramp-wall-labeling` (`RampSystemsv4.tex:388`) sets
`theta_{m_{K+1}} = GB_n`. `eq:mKnjK` (`RectangularGeo.tex:27`) sets
`theta_{m_{K+1}} = GB_n - 1/4` with `h_{m_{K+1}} = 1/4`, so that
`theta + h = GB_n`. Both are internally consistent, but any formula quantified
over `v in prod {1..K(n)}` that involves `theta_{m_{k+1}}` — `I_n`, `Xi_n`, and
hence the `H_2` inequality — picks up whichever convention its chapter uses, and
the two disagree at the outermost cells. Worth one sentence in the text saying
which is meant where.

---

## What changed in the code

`code/DSGRN_utils` now implements the manuscript's definitions:

| | file | change |
|---|---|---|
| C1 | `CubicalBlowupGraph.py:298` | `o_xi` targets restricted to `J_i(xi)` |
| C2 | `CubicalBlowupGraph.py:427` | back-wall pairs outside `1 <= v'_j <= K(j)` discarded instead of wrapping |
| C3 | `WallLabelling.py:44` | outer sentinel is `GB_n`, passed as `global_bound` |
| C4 | `CubicalBlowupGraph.py:470` | `semi_opaque_cell` excludes top cells and `bdy(X)` |
| C5 | `IsomorphismQuery.py:23` | rewritten around a Morse-graph signature; the `continue`/`break` bug fixed |
| — | `CubicalBlowupGraph.py:621` | `F_i` assembled as an intersection of refinements plus an unconditional union with `U(xi)`, per `defn:Rule3`, rather than as a first-match cascade |

Two performance changes were needed to make the faithful assembly affordable,
since the intersection evaluates conditions the cascade skipped:

- the cell-indexed quantities (`top_star`, `gradient_directions`,
  `opaque_directions`, `equilibrium_cell`, the regulation map, semi-opacity, and
  the rook-field components) are memoised — the file's own TODO 1;
- Condition 2.1 is evaluated only where `F_1` returned no verdict. That is not
  an approximation: `cor:F2-well-defined` argues every GO-pair leaves an `F_1`
  double edge, so a pair `F_1` has already oriented is not in `D(Phi)` and the
  decision wall cannot speak about it. Checked on 216,612 such pairs, with a
  `strict_intersection=True` mode that disables the skip and raises on any
  disagreement (exercised by `test_condition_2_1_never_contradicts_F1`).

This also explains why the cascade and the intersection agreed everywhere in the
audit: the two conditions have disjoint domains. Together the changes brought the
six-dimensional runs from twenty minutes back to parity with the legacy path
(65 s).

`legacy=True` restores all of it. `rookfields` additionally provides
`realises_dsgrn_labelling`, which builds a ramp system for any DSGRN parameter
node via `defn:DSGRN_ramp` and checks that it induces that node's wall
labelling — the content of `prop:DSGRN_wall_constant`. It succeeded on every
parameter node tried, across five networks in dimensions 2 and 3.

## Reproducing

```bash
cd code
./.venv/bin/python -m pytest rookfields/tests            # 182 tests
./.venv/bin/python -m rookfields.audit                   # reports/divergence.md
./.venv/bin/python -m rookfields.examples --spec both    # reports/examples.md
./.venv/bin/python -m rookfields.hypotheses              # reports/hypotheses.md
./.venv/bin/python -m rookfields.alignment               # reports/alignment.md
./.venv/bin/python -m rookfields.plotting                # geometrization figures
./.venv/bin/python -m rookfields.plotting.gallery        # Morse graphs/sets + 2-D gallery
./.venv/bin/python -m rookfields.report_pdf              # reports/updated-computations.pdf
```
