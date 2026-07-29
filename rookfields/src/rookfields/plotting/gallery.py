"""Generate the full figure set: example comparisons and the 2-D gallery.

    python -m rookfields.plotting.gallery

Two families of figure:

``examples/``
    for each catalogued example, the Morse graph (with Conley indices and the
    FP/PO/T/M classification) and the Morse sets with the state transition
    graph, under both algorithm specifications side by side;

``gallery2d/``
    for a variety of two-node DSGRN parameter nodes, the rectangular
    geometrization drawn in real ramp coordinates with the vector field, the
    nullclines, the GO-manifolds and the orange connectors, next to the Morse
    graph and Morse sets for the same node.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .. import networks  # noqa: E402
from ..examples.catalog import EXAMPLES  # noqa: E402
from ..pipeline import IDENTICAL, compare_results, conley_morse_graph  # noqa: E402
from ..ramp import RampSystem  # noqa: E402
from ..spec import LEGACY, PAPER  # noqa: E402
from .morse import comparison_figure, draw_morse_sets, render_morse_graph  # noqa: E402
from .phase2d import plot_phase_portrait_2d  # noqa: E402

FIGURES = Path(__file__).resolve().parents[3] / "figures"

#: Two-node DSGRN parameter nodes for the geometrization gallery.
GALLERY_2D = [
    ("N2_b", 974, "ex:mccord -- the running wall labelling of fig:wall_labeling(A)"),
    ("N2_a", 752, "ex:saddlesaddlebif -- two admissible connection matrices"),
    ("N2_a", 47, "ex:trivial_index_2D -- one trivial-index node"),
    ("toggle", 4, "toggle switch -- bistable"),
    ("N2_b", 1376, "the richest GO structure found: 6 pairs, 3 triples"),
    ("N2_a", 968, "several GO pairs and a three-node Morse graph"),
]

RAMP_WIDTH = 0.05


def _result_pair(example, level=None):
    level = level or example.level
    out = {}
    for spec in (LEGACY, PAPER):
        if example.is_ramp:
            system = networks.RAMP_SYSTEMS[example.ramp_system]
            labelling, K = system.wall_labelling()
            out[spec.name] = conley_morse_graph(
                labelling=labelling, num_thresholds=K, spec=spec, level=level
            )
        else:
            out[spec.name] = conley_morse_graph(
                networks.parameter(example.network, example.index),
                spec=spec,
                level=level,
            )
    return out


def generate_example_figures(out_dir: Path, *, include_expensive: bool = False) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for example in EXAMPLES:
        if example.expensive and not include_expensive:
            print(f"  skipping {example.label} (expensive)")
            continue
        print(f"  {example.label} ...", flush=True)
        results = _result_pair(example)
        target = out_dir / f"{example.label.replace(':', '_')}.png"
        subtitle = (
            f"{example.network} node {example.index:,}"
            if example.network
            else f"ramp system {example.ramp_system}"
        )
        comparison_figure(
            results,
            target,
            title=f"{example.label}   ({example.source})",
            subtitle=subtitle + f",  $\\mathcal{{F}}_{{{example.level}}}$",
        )
        verdict = compare_results(results["legacy"], results["paper"])
        differs = verdict != IDENTICAL
        written.append(
            {
                "label": example.label,
                "file": target.name,
                "source": example.source,
                "level": example.level,
                "dim": results["paper"].stg.dim,
                "differs": differs,
                "verdict": verdict,
                "legacy_nodes": len(results["legacy"].nodes),
                "paper_nodes": len(results["paper"].nodes),
            }
        )
    return written


def generate_2d_gallery(out_dir: Path, cases=None) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cases = cases or GALLERY_2D
    written = []
    for name, index, caption in cases:
        print(f"  {name}[{index}] ...", flush=True)
        spec, realised = networks.realises_dsgrn_labelling(name, index, RAMP_WIDTH)
        if not realised:
            print("    no ramp realisation found for this region; skipped")
            continue
        system = RampSystem(spec)

        from .. import r2_manifolds as R2

        labelling, K = system.wall_labelling()
        result = conley_morse_graph(
            labelling=labelling, num_thresholds=K, spec=PAPER, level=3
        )
        pairs, stg = R2.go_pairs(system, spec=PAPER)
        triples = R2.d3_triples(stg)

        fig = plt.figure(figsize=(13.5, 12.6))
        gs = fig.add_gridspec(2, 2, hspace=0.22, wspace=0.20)

        ax = fig.add_subplot(gs[0, 0])
        plot_phase_portrait_2d(system, ax=ax, spec=PAPER, level=3)
        ax.set_title(
            "Rectangular geometrization with the ramp field\n"
            "black: GO-manifolds (solid external, dashed internal); "
            "orange: connectors",
            fontsize=10,
        )

        # a detail view around the GO / orange cluster
        ax = fig.add_subplot(gs[0, 1])
        # Zoom on a single codimension-two corner: the first triple if there is
        # one (so the orange connector is in frame), otherwise the first pair.
        points = []
        if triples:
            focus = triples[0]
            orange = R2.build_orange_manifold(system, focus)
            points.extend(orange.anchor_g[:, :2])
            points.extend(orange.anchor_o[:, :2])
            points.extend(orange.sigma[:, :2])
            for pair in pairs:
                if pair.face == focus.face:
                    manifold = R2.build_go_manifold(system, pair, base_samples=3)
                    for traj in manifold.trajectories:
                        points.extend(traj[:, :2])
        elif pairs:
            manifold = R2.build_go_manifold(system, pairs[0], base_samples=3)
            for traj in manifold.trajectories:
                points.extend(traj[:, :2])
        if points:
            import numpy as np

            pts = np.array(points)
            lo, hi = pts.min(axis=0), pts.max(axis=0)
            pad = np.maximum((hi - lo) * 1.6, 0.04)
            plot_phase_portrait_2d(
                system, ax=ax, spec=PAPER, level=3, field_grid=16,
                xlim=(lo[0] - pad[0], hi[0] + pad[0]),
                ylim=(lo[1] - pad[1], hi[1] + pad[1]),
            )
            ax.set_title(
                "Detail: the codimension-two corners where the\n"
                "GO-manifolds and orange connectors sit",
                fontsize=10,
            )
        else:
            ax.text(0.5, 0.5, "no GO pairs at this parameter node",
                    ha="center", va="center", fontsize=11)
            ax.axis("off")

        ax = fig.add_subplot(gs[1, 0])
        draw_morse_sets(result, ax)
        ax.set_title("Morse sets and the state transition graph", fontsize=11)

        ax = fig.add_subplot(gs[1, 1])
        import matplotlib.image as mpimg
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            png = render_morse_graph(result, Path(tmp) / "mg.png")
            ax.imshow(mpimg.imread(png))
        ax.axis("off")
        ax.set_title("Morse graph with Conley indices", fontsize=11)

        fig.suptitle(
            f"{name} node {index:,} --- {caption}\n"
            f"{len(pairs)} codimension-2 GO pairs "
            f"({sum(1 for p in pairs if p.external)} external), "
            f"{len(triples)} $\\mathcal{{D}}_3$ triples,  "
            f"$h={RAMP_WIDTH}$",
            fontsize=13,
        )
        target = out_dir / f"gallery_{name}_{index}.png"
        fig.savefig(target, dpi=115, bbox_inches="tight")
        plt.close(fig)
        written.append(
            {
                "network": name,
                "index": index,
                "caption": caption,
                "file": target.name,
                "go_pairs": len(pairs),
                "external": sum(1 for p in pairs if p.external),
                "d3_triples": len(triples),
                "morse_nodes": len(result.nodes),
            }
        )
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=FIGURES)
    parser.add_argument("--include-expensive", action="store_true")
    parser.add_argument("--only", choices=["examples", "gallery"], default=None)
    args = parser.parse_args(argv)

    manifest = {}
    if args.only in (None, "examples"):
        print("example comparisons (legacy vs paper)")
        manifest["examples"] = generate_example_figures(
            args.out / "examples", include_expensive=args.include_expensive
        )
    if args.only in (None, "gallery"):
        print("two-dimensional geometrization gallery")
        manifest["gallery2d"] = generate_2d_gallery(args.out / "gallery2d")

    path = args.out / "manifest.json"
    existing = json.loads(path.read_text()) if path.exists() else {}
    existing.update(manifest)
    path.write_text(json.dumps(existing, indent=2))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
