# rookfields

Companion code for `paper/Rook_Field_Paper_v2`. Three jobs:

1. compute the multivalued maps `F_0 .. F_4` under both the `DSGRN_utils`
   implementation and the live manuscript definitions, and diff them;
2. recompute every example in the monograph and compare against the claimed
   output;
3. construct, render, and numerically validate the geometrizations of Part III.

Read [`reports/FINDINGS.md`](reports/FINDINGS.md) first — it is the summary of
what all of this turned up.

`paper/` is never modified. Regenerated figures land in `figures/` for the
authors to copy across if they want them.

## Status

`code/DSGRN_utils` now implements the manuscript's definitions by default.
`CubicalBlowupGraph(..., legacy=True)` and `ConleyMorseGraph(..., legacy=True)`
restore the superseded behaviour, so figures published before the alignment
still reproduce exactly. The changes and their measured effect are in
[`reports/FINDINGS.md`](reports/FINDINGS.md); the recomputed examples, with
Morse graphs and Morse sets under both readings, are in
[`reports/updated-computations.pdf`](reports/updated-computations.pdf).

## Environment

The venv at `code/.venv` was built by a pyenv interpreter that no longer exists;
it has been re-pointed at Homebrew's Python 3.13 (the site-packages tree,
including the DSGRN and pyCHomP2 C extensions, was intact and needed no rebuild).

```bash
cd code
./.venv/bin/python -c "import DSGRN, pychomp, DSGRN_utils, rookfields; print('ok')"
```

`DSGRN_utils` is installed editable from `code/DSGRN_utils` (the current fork).
The stale duplicate that the old editable install pointed at has been moved to
`code/_archive/src-DSGRN_utils-stale/`, along with the legacy `setup.py` that
packaged it.

## Layout

| module | what it does |
|---|---|
| `spec.py` | every paper/code divergence as an independent toggle; `LEGACY` and `PAPER` presets |
| `blowup.py` | `SpecCubicalBlowupGraph`, subclassing `DSGRN_utils.CubicalBlowupGraph` and overriding only what the specs touch |
| `pipeline.py` | wall labelling → `F_i` → SCC → Morse graph → Conley complex, returning the `scc_dag` and connection matrix that `ConleyMorseGraph` discards |
| `wall_labeling.py` | `def:wall_labeling` validity, strong dissipativity, and enumeration of *abstract* wall labelings |
| `networks.py` | the monograph's networks and ramp systems, cross-checked against the parameter-graph sizes in the text |
| `ramp.py` | ramp vector field, flow, global bound, `I_n`/`Xi_n`/`L_n`/`U_n`, and `defn:ramp-wall-labeling` |
| `hypotheses.py` | `Lambda(R)`, `H_0`, `H_1`, `H_2`, `H_3`, and the `prop:smallh` bounds |
| `geometrization.py` | rectangular geometrization of `X_b` and the corrected Janus subdivision |
| `alignment.py` | numerical inward-crossing check — the content of `thm:R1ABlattice` |
| `r2_manifolds.py` | GO-manifolds of `defn:GO-manifold`, with the R2 sign and product checks |
| `r3_manifolds.py` | cycle cells, the cyclic-feedback spectrum, swept-surface tangent/normal |
| `audit.py` | LEGACY vs PAPER diff across networks and levels |
| `examples/` | the 17 computational examples with their claimed outputs |
| `plotting/` | Morse graphs and Morse sets (DSGRN style), 2D/3D rendering in real ramp coordinates, and the geometrization gallery |
| `report_pdf.py` | the PDF report of recomputed examples and open problems |

## Usage

```python
from rookfields import conley_morse_graph, LEGACY, PAPER
from rookfields import networks

p = networks.parameter("N3_B", 2_472_287)
result = conley_morse_graph(p, spec=PAPER, level=4)
result.nodes, result.conley_indices, result.morse_set_cells(1)
```

`LEGACY` reproduces `DSGRN_utils` bit-for-bit, so anything computed with it
matches the previously published figures exactly. `PAPER` applies the live
manuscript definitions; `rookfields.spec.ISOLATED` holds one-toggle-at-a-time
variants so a difference can be attributed to a single divergence.

Directly specified ramp systems bypass DSGRN:

```python
from rookfields.networks import VAN_DER_POL
from rookfields.ramp import RampSystem
from rookfields import hypotheses as H

system = RampSystem(VAN_DER_POL)
for verdict in H.report(system):
    print(verdict)
```

## Commands

```bash
./.venv/bin/python -m pytest rookfields/tests            # 182 tests
./.venv/bin/python -m rookfields.audit                   # reports/divergence.md
./.venv/bin/python -m rookfields.examples --spec both    # reports/examples.md
./.venv/bin/python -m rookfields.hypotheses              # reports/hypotheses.md
./.venv/bin/python -m rookfields.alignment               # reports/alignment.md
./.venv/bin/python -m rookfields.plotting                # geometrization figures
./.venv/bin/python -m rookfields.plotting.gallery        # Morse graphs/sets + 2-D gallery
./.venv/bin/python -m rookfields.report_pdf              # reports/updated-computations.pdf
```

`rookfields.plotting.gallery` writes two families of figure: `figures/examples/`
holds a Morse graph and Morse set comparison for every catalogued example under
both readings, and `figures/gallery2d/` holds the two-dimensional
geometrizations with their GO-manifolds and orange connectors.

Any DSGRN parameter node can be turned into an explicit ramp system, which is
what makes the geometrization pictures possible:

```python
from rookfields.networks import realises_dsgrn_labelling
from rookfields.ramp import RampSystem

spec, ok = realises_dsgrn_labelling("N2_b", 974, h_bar=0.05)
system = RampSystem(spec)   # ok is True when it induces that node's wall labelling
```

`--skip-expensive` on the examples runner drops the two six-dimensional EMT
runs, which dominate wall clock at roughly 90 s and 150 s each per spec.

## Tests

The suite is written against the manuscript, not against the implementation:
each test names the statement it exercises. Notable ones are the exhaustive
`thm:2dDone` sweep over abstract wall labelings, the `prop:CH=0` check that
`MorseGraph`'s gradient pruning depends on, and the R2 checks of
`prop:GO-crossing-sign` and `lem:GO-product`.

Tests that encode a *finding* rather than a passing claim say so in their
docstring — for instance `test_published_widths_violate_H1` and
`test_h3_bound_does_not_force_a_complex_pair`.
