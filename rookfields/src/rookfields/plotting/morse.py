"""DSGRN-style Morse graph and Morse set rendering, and legacy/paper comparisons.

``DSGRN_utils.PlotMorseGraph`` and ``DSGRN_utils.PlotMorseSets`` take exactly the
objects a :class:`~rookfields.pipeline.MorseResult` carries, so they are used
directly rather than reimplemented: the plots are the ones the monograph's
figures were made with.

The comparison figures put the two algorithm specifications side by side so a
difference in the Morse graph, the Conley indices, or the Morse sets themselves
is visible at a glance.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from DSGRN_utils.PlotMorseGraph import PlotMorseGraph  # noqa: E402
from DSGRN_utils.PlotMorseSets import PlotMorseSets  # noqa: E402

from ..pipeline import IDENTICAL, MorseResult, compare_results  # noqa: E402

#: `PlotMorseGraph` node annotations: 'c' Conley index, 's' size,
#: 'a' FP/PO/T/M classification.  Colons concatenate them.
DEFAULT_LABEL = "c:a"


def render_morse_graph(result: MorseResult, path: Path, *, label: str = DEFAULT_LABEL) -> Path:
    """Render the Morse graph to PNG via graphviz.

    Node labels carry the Conley index and the FP/PO/T/M classification that
    ``PlotMorseGraph.annotation`` derives from it (fixed point, periodic orbit,
    trivial, or a more complex invariant set).
    """
    source = PlotMorseGraph(result.morse_graph, label=label)
    path = Path(path)
    with tempfile.TemporaryDirectory() as tmp:
        dot = Path(tmp) / "mg.dot"
        dot.write_text(source.source)
        subprocess.run(
            ["dot", "-Tpng", "-Gdpi=170", str(dot), "-o", str(path)],
            check=True,
            capture_output=True,
        )
    return path


def draw_morse_sets(result: MorseResult, ax, *, proj_dims=None, proj_slice=None, **kw):
    """Draw the Morse sets and the state transition graph on ``ax``.

    Delegates to ``DSGRN_utils.PlotMorseSets``.  For dimension greater than two
    a projection is drawn; without a slice the cells collapse onto the chosen
    plane and each projected box takes the colour of the lowest Morse node it
    meets, which is that function's documented behaviour.
    """
    if proj_dims is None and result.stg.dim > 2:
        proj_dims = (0, 1)
    PlotMorseSets(
        result.morse_graph,
        result.stg,
        result.graded_complex,
        proj_dims=proj_dims,
        proj_slice=proj_slice,
        ax=ax,
        **kw,
    )
    return ax


def _conley_caption(result: MorseResult) -> str:
    parts = []
    for node in result.nodes:
        index = result.conley_indices[node]
        cells = result.top_cell_counts.get(node, 0)
        parts.append(f"{node}: {tuple(index)} [{cells}]")
    return ", ".join(parts)


#: Above this dimension the Morse set panel is dropped.  ``PlotMorseSets``
#: iterates every top cell of the blowup complex and its adjacencies, which for
#: a six-node network is a few hundred thousand cells; the resulting projection
#: onto two coordinates is unreadable anyway.
MAX_MORSE_SET_DIM = 4


def comparison_figure(
    results: dict[str, MorseResult],
    path: Path,
    *,
    title: str = "",
    subtitle: str = "",
    proj_dims=None,
    proj_slice=None,
) -> Path:
    """Morse graph and Morse sets for each spec, laid out for comparison.

    ``results`` maps a spec name to its :class:`MorseResult`.  Above
    :data:`MAX_MORSE_SET_DIM` only the Morse graphs are drawn.
    """
    path = Path(path)
    specs = list(results)
    n = len(specs)
    dim = next(iter(results.values())).stg.dim
    with_sets = dim <= MAX_MORSE_SET_DIM

    graph_pngs = {}
    with tempfile.TemporaryDirectory() as tmp:
        for spec in specs:
            graph_pngs[spec] = render_morse_graph(
                results[spec], Path(tmp) / f"{spec}.png"
            )

        rows = 2 if with_sets else 1
        fig, axes = plt.subplots(
            rows, n,
            figsize=(7.2 * n, 12.0 if with_sets else 6.0),
            gridspec_kw={"height_ratios": [1.05, 1.4]} if with_sets else None,
            squeeze=False,
        )
        if not with_sets:
            axes = axes.reshape(1, n)

        for col, spec in enumerate(specs):
            result = results[spec]

            ax = axes[0, col]
            ax.imshow(mpimg.imread(graph_pngs[spec]))
            ax.axis("off")
            differs = ""
            if n == 2 and col == 1:
                verdict = compare_results(results[specs[0]], result)
                if verdict != IDENTICAL:
                    differs = f"  ({verdict})"
            ax.set_title(
                f"$\\mathcal{{F}}_{{{result.level}}}$ Morse graph --- "
                f"{spec}{differs}",
                fontsize=13,
            )

            if not with_sets:
                continue

            ax = axes[1, col]
            try:
                draw_morse_sets(result, ax, proj_dims=proj_dims, proj_slice=proj_slice)
            except Exception as exc:  # pragma: no cover - defensive
                ax.text(0.5, 0.5, f"Morse sets unavailable:\n{exc}",
                        ha="center", va="center", fontsize=9)
                ax.axis("off")
            projected = result.stg.dim > 2
            ax.set_title(
                "Morse sets and the state transition graph"
                + (f"\nprojected to coordinates {proj_dims or (0, 1)}" if projected else ""),
                fontsize=11,
            )

        caption = "\n".join(
            f"{spec}:  " + _conley_caption(results[spec]) for spec in specs
        )
        fig.suptitle(
            (title or "") + (f"\n{subtitle}" if subtitle else ""),
            fontsize=15,
            y=0.995,
        )
        fig.text(
            0.5, 0.012,
            "Morse node: (Conley index) [cells]\n" + caption,
            ha="center", va="bottom", fontsize=9, family="monospace",
        )
        fig.tight_layout(rect=(0, 0.055 if with_sets else 0.10, 1, 0.965))
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return path


def morse_sets_figure(
    result: MorseResult,
    path: Path,
    *,
    title: str = "",
    proj_dims=None,
    proj_slice=None,
    fig_size: float = 8.0,
) -> Path:
    """A single Morse set plot, DSGRN style."""
    path = Path(path)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    draw_morse_sets(result, ax, proj_dims=proj_dims, proj_slice=proj_slice)
    if title:
        ax.set_title(title, fontsize=12)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def graphviz_available() -> bool:
    return shutil.which("dot") is not None
