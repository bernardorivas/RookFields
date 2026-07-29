"""Diff the DSGRN_utils reading of F_i against the live manuscript reading.

Run ``python -m rookfields.audit`` to regenerate ``reports/divergence.md``.

Each toggle in :mod:`rookfields.spec` is measured on its own, so a difference
can be attributed to a single divergence rather than to the bundle.  For every
(network, level, spec) we record how many parameters change their state
transition graph, their Morse-graph node count, or their Conley indices.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import random
import time
from collections import Counter
from pathlib import Path

from . import networks
from .blowup import SpecCubicalBlowupGraph
from .pipeline import conley_morse_graph
from .spec import ISOLATED, LEGACY, PAPER, Spec

REPORTS = Path(__file__).resolve().parents[2] / "reports"

#: (network, sample size).  ``None`` means sweep the whole parameter graph.
DEFAULT_FAMILIES: list[tuple[str, int | None]] = [
    ("toggle", None),
    ("repressilator", None),
    ("cycle3", 120),
    ("cycle4", 60),
    ("N2_a", 200),
    ("N2_b", 200),
    ("N3_B", 80),
]

DEFAULT_LEVELS = (1, 2, 3, 4)


@dataclasses.dataclass
class Divergence:
    network: str
    dim: int
    level: int
    spec: str
    sampled: int
    stg_changed: int
    morse_nodes_changed: int
    conley_changed: int
    diagnostics: Counter

    def as_row(self) -> dict:
        d = dataclasses.asdict(self)
        d["diagnostics"] = dict(self.diagnostics)
        return d


def sample_indices(name: str, size: int | None, seed: int) -> list[int]:
    total = networks.parameter_graph(name).size()
    if size is None or total <= size:
        return list(range(total))
    return sorted(random.Random(seed).sample(range(total), size))


def _signature(result) -> tuple:
    """What we compare: STG edges, Morse nodes, Conley indices."""
    return (
        frozenset(result.stg.digraph.edges()),
        tuple(result.nodes),
        tuple(sorted(result.conley_indices.items())),
    )


def compare(
    name: str,
    indices: list[int],
    level: int,
    spec: Spec,
    baseline: dict[int, tuple],
) -> Divergence:
    net = networks.network(name)
    diagnostics: Counter = Counter()
    stg_changed = nodes_changed = conley_changed = 0

    for i in indices:
        result = conley_morse_graph(
            networks.parameter(name, i), spec=spec, level=level
        )
        diagnostics.update(result.diagnostics)
        edges, nodes, conley = _signature(result)
        b_edges, b_nodes, b_conley = baseline[i]
        if edges != b_edges:
            stg_changed += 1
        if nodes != b_nodes:
            nodes_changed += 1
        if conley != b_conley:
            conley_changed += 1

    return Divergence(
        network=name,
        dim=net.size(),
        level=level,
        spec=spec.name,
        sampled=len(indices),
        stg_changed=stg_changed,
        morse_nodes_changed=nodes_changed,
        conley_changed=conley_changed,
        diagnostics=diagnostics,
    )


def unstable_edge_suppression(name: str, indices: list[int], level: int) -> dict:
    """D3: does the cascade ever drop a ``U(xi)`` edge the definition demands?

    ``defn:Rule3`` unions ``U(xi)`` in unconditionally, but
    ``compute_multivalued_map`` only reaches it for pairs that both F_1 and the
    decision wall left unoriented.
    """
    missing = 0
    cells_with_u = 0
    affected_parameters = 0
    for i in indices:
        stg = SpecCubicalBlowupGraph(
            networks.parameter(name, i), spec=LEGACY, level=level
        )
        local = 0
        for cc_cell in stg.cubical_complex:
            if stg.cubical_complex.rightfringe(cc_cell):
                continue
            targets = stg.unstable_targets(cc_cell)
            if not targets:
                continue
            cells_with_u += 1
            source = stg.cubical2blowup(cc_cell)
            adjacent = set(stg.digraph.adjacencies(source))
            local += sum(
                1 for t in targets if stg.cubical2blowup(t) not in adjacent
            )
        missing += local
        if local:
            affected_parameters += 1
    return {
        "network": name,
        "level": level,
        "sampled": len(indices),
        "cells_with_nonempty_U": cells_with_u,
        "missing_U_edges": missing,
        "parameters_affected": affected_parameters,
    }


def run(
    families=DEFAULT_FAMILIES,
    levels=DEFAULT_LEVELS,
    specs=None,
    seed: int = 20260728,
) -> dict:
    specs = specs or {**ISOLATED, "paper": PAPER}
    rows: list[Divergence] = []
    suppression: list[dict] = []
    started = time.time()

    for name, size in families:
        indices = sample_indices(name, size, seed)
        print(f"[audit] {name}: {len(indices)} parameters, dim {networks.network(name).size()}")
        for level in levels:
            baseline = {}
            for i in indices:
                baseline[i] = _signature(
                    conley_morse_graph(networks.parameter(name, i), spec=LEGACY, level=level)
                )
            for spec in specs.values():
                row = compare(name, indices, level, spec, baseline)
                rows.append(row)
                print(
                    f"    level={level} {spec.name:24s} "
                    f"stg={row.stg_changed:4d} nodes={row.morse_nodes_changed:4d} "
                    f"conley={row.conley_changed:4d}"
                )
            if level >= 3:
                suppression.append(unstable_edge_suppression(name, indices, level))

    return {
        "generated_seconds": time.time() - started,
        "seed": seed,
        "rows": [r.as_row() for r in rows],
        "unstable_edge_suppression": suppression,
    }


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

_TOGGLE_NOTES = {
    "D1-regulation-codomain": (
        "`o_xi` targets restricted to `J_i(xi)` per live `def:active_regulation` "
        "(RookFields6.tex:738); DSGRN_utils ranges over all N directions, which is "
        "the retired `def:active_regulation-old` (:726)."
    ),
    "D1b-regulation-domain": (
        "`o_xi` domain widened from interior-inessential to all inessential "
        "directions (CubicalBlowupGraph.py:302 has no counterpart in the text)."
    ),
    "D2-back-wall-strict": (
        "Pairs discarded unless `1 <= v'_j <= K(j)` for all `j in J_i(xi')`, the "
        "condition under which every member of `Back(xi,xi')` (`defn:back-walls`, "
        "CombinatorialDynamics.tex:356) is a cell of X."
    ),
    "D2-back-wall-shifted": (
        "Weaker variant: discard only when the chosen shift itself leaves X."
    ),
    "D2b-back-wall-all": (
        "Evaluate every `r_n in R_n(xi)` choice rather than the single coherent "
        "one DSGRN_utils uses -- an empirical test of "
        "`prop:back_wall_well_defined` (CombinatorialDynamics.tex:762)."
    ),
    "D3-f3-composition": (
        "`F_3 = (F_2 cap F_3.1) u U` computed as an intersection plus an "
        "unconditional union (`defn:Rule3`, CombinatorialDynamics.tex:1102) "
        "rather than as the first-match cascade of "
        "`compute_multivalued_map` (:621)."
    ),
    "D4-semi-opaque-guard": (
        "`defn:partially_opaque` (CombinatorialDynamics.tex:999) excludes top "
        "cells and `bdy(X)`; `semi_opaque_cell` (:470) checks only that `o_xi` "
        "is a bijection."
    ),
    "paper": "All live-manuscript readings at once.",
}


def render(data: dict) -> str:
    rows = data["rows"]
    out: list[str] = []
    out.append("# F_i divergence audit: DSGRN_utils vs. the live manuscript\n")
    out.append(
        "Generated by `python -m rookfields.audit`. Each row counts parameters "
        "whose result changes relative to the `legacy` spec, which reproduces "
        "`DSGRN_utils` bit-for-bit.\n"
    )
    out.append(
        f"Sampling seed `{data['seed']}`; total runtime "
        f"{data['generated_seconds']:.1f} s.\n"
    )

    families = []
    for r in rows:
        key = (r["network"], r["dim"], r["sampled"])
        if key not in families:
            families.append(key)

    # -- what actually differs -------------------------------------------
    nonzero = [
        r
        for r in rows
        if r["stg_changed"] or r["morse_nodes_changed"] or r["conley_changed"]
    ]
    out.append("\n## Everything that differs\n")
    if not nonzero:
        out.append("No toggle changed any result on any sampled parameter.\n")
    else:
        out.append("| network | dim | parameters | level | toggle | STG | Morse nodes | Conley indices |")
        out.append("|---|---:|---:|---:|---|---:|---:|---:|")
        for r in nonzero:
            out.append(
                f"| {r['network']} | {r['dim']} | {r['sampled']} | {r['level']} | "
                f"`{r['spec']}` | {r['stg_changed']} | {r['morse_nodes_changed']} | "
                f"{r['conley_changed']} |"
            )
    silent = sorted({r["spec"] for r in rows} - {r["spec"] for r in nonzero})
    if silent:
        out.append(
            "\nToggles that changed nothing anywhere: "
            + ", ".join(f"`{s}`" for s in silent)
            + ".\n"
        )

    out.append("\n## Full grid\n")
    for name, dim, sampled in families:
        out.append(f"\n### {name} (dim {dim}, {sampled} parameters)\n")
        out.append("| toggle | level | STG changed | Morse nodes changed | Conley indices changed |")
        out.append("|---|---:|---:|---:|---:|")
        for r in rows:
            if r["network"] != name:
                continue
            out.append(
                f"| `{r['spec']}` | {r['level']} | {r['stg_changed']} | "
                f"{r['morse_nodes_changed']} | {r['conley_changed']} |"
            )

    out.append("\n## What each toggle changes\n")
    for k, v in _TOGGLE_NOTES.items():
        out.append(f"- **`{k}`** — {v}")

    out.append("\n## D3: unconditional union with U(xi)\n")
    out.append(
        "`defn:Rule3` unions `U(xi)` into `F_3` for every semi-opaque cell with a "
        "nontrivial regulation cycle. The DSGRN_utils cascade only reaches that "
        "step for pairs left unoriented by both F_1 and the decision wall. A "
        "nonzero `missing_U_edges` would be a genuine omission.\n"
    )
    out.append("| network | level | parameters | cells with U(xi) != 0 | missing U edges | parameters affected |")
    out.append("|---|---:|---:|---:|---:|---:|")
    for s in data["unstable_edge_suppression"]:
        out.append(
            f"| {s['network']} | {s['level']} | {s['sampled']} | "
            f"{s['cells_with_nonempty_U']} | {s['missing_U_edges']} | "
            f"{s['parameters_affected']} |"
        )

    out.append("\n## Diagnostics\n")
    out.append("| network | level | toggle | counter | total |")
    out.append("|---|---:|---|---|---:|")
    for r in rows:
        for k, v in sorted(r["diagnostics"].items()):
            if not v:
                continue
            out.append(f"| {r['network']} | {r['level']} | `{r['spec']}` | `{k}` | {v} |")

    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", action="append", help="restrict to these networks")
    parser.add_argument("--sample", type=int, default=None, help="override sample size")
    parser.add_argument("--level", action="append", type=int, help="restrict to these levels")
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument("--out", type=Path, default=REPORTS / "divergence.md")
    parser.add_argument("--json", type=Path, default=REPORTS / "divergence.json")
    args = parser.parse_args(argv)

    families = DEFAULT_FAMILIES
    if args.network:
        families = [(n, args.sample) for n in args.network]
    elif args.sample is not None:
        families = [(n, args.sample) for n, _ in DEFAULT_FAMILIES]
    levels = tuple(args.level) if args.level else DEFAULT_LEVELS

    data = run(families=families, levels=levels, seed=args.seed)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render(data))
    args.json.write_text(json.dumps(data, indent=2, default=str))
    print(f"\nwrote {args.out}\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
