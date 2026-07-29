"""Planar phase portraits over the geometrized complex."""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from ..alignment import check_alignment  # noqa: E402
from ..geometrization import RectangularGeometrization  # noqa: E402
from ..pipeline import conley_morse_graph  # noqa: E402
from ..ramp import RampSystem  # noqa: E402
from ..spec import PAPER  # noqa: E402

#: Same palette as DSGRN_utils.PlotMorseSets, so figures stay comparable.
MORSE_COLOURS = [
    "#1f77b4", "#e6550d", "#31a354", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#80b1d3", "#ffffb3", "#fccde5",
    "#b3de69", "#fdae6b", "#6a3d9a", "#c49c94",
]


def plot_phase_portrait_2d(
    system: RampSystem,
    *,
    spec=PAPER,
    level: int = 3,
    ax=None,
    show_field: bool = True,
    show_nullclines: bool = True,
    show_cells: bool = True,
    show_go: bool = True,
    trajectories: list | None = None,
    field_grid: int = 26,
    title: str | None = None,
    zoom: str = "thresholds",
    xlim=None,
    ylim=None,
):
    """Draw the geometrized complex, Morse sets, field, nullclines and GO curves.

    Everything is in real ramp coordinates, so the combinatorial objects and the
    ODE objects live on the same axes and can be compared directly.
    """
    assert system.dim == 2, "planar portrait needs a two-dimensional system"
    geo = RectangularGeometrization(system)
    labelling, K = system.wall_labelling()
    result = conley_morse_graph(
        labelling=labelling, num_thresholds=K, spec=spec, level=level
    )
    stg = result.stg

    if ax is None:
        _fig, ax = plt.subplots(figsize=(7.5, 7.0))

    gb = system.global_bound

    # The outermost boxes run all the way to GB_n, which for DSGRN-sampled
    # parameters is far beyond the thresholds; drawing the whole box would leave
    # every cell of interest in one corner.  Default to the threshold range.
    if xlim is None or ylim is None:
        if zoom == "thresholds":
            span = []
            for n in range(2):
                top = system.sorted_thresholds[n][-1] + system.sorted_widths[n][-1]
                span.append((0.0, min(float(gb[n]), top * 1.25)))
        else:
            span = [(0.0, float(gb[n])) for n in range(2)]
        xlim = xlim or span[0]
        ylim = ylim or span[1]

    # -- Morse sets, drawn as the embedded blowup cells -------------------
    if show_cells:
        vertex_index = {
            v: result.morse_graph.vertex_label(v)[0]
            for v in result.morse_graph.vertices()
        }
        for cell in stg.blowup_complex(stg.dim):
            if stg.blowup_complex.rightfringe(cell):
                continue
            grade = result.graded_complex.value(cell)
            if grade not in vertex_index:
                continue
            node = vertex_index[grade]
            coords = stg.blowup_complex.coordinates(cell)
            (x0, x1), (y0, y1) = geo.rectangle(list(coords), [1, 1])
            ax.add_patch(
                Rectangle(
                    (x0, y0),
                    x1 - x0,
                    y1 - y0,
                    facecolor=MORSE_COLOURS[node % len(MORSE_COLOURS)],
                    edgecolor="none",
                    alpha=0.55,
                    zorder=1,
                )
            )

    # -- the cell grid ----------------------------------------------------
    for n, axis_draw in ((0, ax.axvline), (1, ax.axhline)):
        top = 2 * (system.num_thresholds[n] + 1) + 1
        window = xlim if n == 0 else ylim
        for c in range(top + 1):
            value = geo.embed_coordinate(n, c)
            if window[0] <= value <= window[1]:
                axis_draw(value, color="0.7", lw=0.5, zorder=0)

    # -- vector field -----------------------------------------------------
    if show_field:
        xs = np.linspace(xlim[0] + 0.02 * (xlim[1] - xlim[0]), 0.98 * xlim[1], field_grid)
        ys = np.linspace(ylim[0] + 0.02 * (ylim[1] - ylim[0]), 0.98 * ylim[1], field_grid)
        X, Y = np.meshgrid(xs, ys)
        U = np.zeros_like(X)
        V = np.zeros_like(Y)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                f = system.vector_field([X[i, j], Y[i, j]])
                U[i, j], V[i, j] = f[0], f[1]
        norm = np.hypot(U, V)
        norm[norm == 0] = 1.0
        ax.quiver(X, Y, U / norm, V / norm, color="0.35", alpha=0.5,
                  width=0.0022, scale=42, zorder=2)

    # -- nullclines -------------------------------------------------------
    if show_nullclines:
        xs = np.linspace(xlim[0], xlim[1], 420)
        ys = np.linspace(ylim[0], ylim[1], 420)
        X, Y = np.meshgrid(xs, ys)
        F0 = np.zeros_like(X)
        F1 = np.zeros_like(X)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                f = system.vector_field([X[i, j], Y[i, j]])
                F0[i, j], F1[i, j] = f[0], f[1]
        ax.contour(X, Y, F0, levels=[0.0], colors="#0b6fa4", linewidths=1.8, zorder=4)
        ax.contour(X, Y, F1, levels=[0.0], colors="#b03060", linewidths=1.8, zorder=4)

    # -- GO-manifolds and orange connectors --------------------------------
    if show_go:
        from .. import r2_manifolds as R2

        pairs, stg2 = R2.go_pairs(system, spec=spec)
        for pair in pairs:
            manifold = R2.build_go_manifold(system, pair, base_samples=3)
            colour = "#111111" if pair.external else "#5b5b5b"
            for traj in manifold.trajectories:
                ax.plot(traj[:, 0], traj[:, 1], color=colour, lw=2.4,
                        zorder=6, solid_capstyle="round",
                        ls="-" if pair.external else (0, (4, 2)))
            if manifold.base.size:
                ax.plot(manifold.base[:, 0], manifold.base[:, 1], ".",
                        color=colour, ms=6, zorder=7)

        # defn:orange-manifold: the affine connector for a D_3 triple
        for triple in R2.d3_triples(stg2):
            orange = R2.build_orange_manifold(system, triple)
            for a, b in zip(orange.anchor_g, orange.anchor_o):
                ax.plot([a[0], b[0]], [a[1], b[1]],
                        color="#ff8c00", lw=3.0, zorder=8, solid_capstyle="round")
            ax.plot(orange.sigma[:, 0], orange.sigma[:, 1], "s",
                    color="#ff8c00", ms=4, zorder=9)

    # -- trajectories ------------------------------------------------------
    for x0 in trajectories or []:
        sol = system.flow(x0, t_span=(0.0, 60.0))
        ts = np.linspace(0.0, sol.t[-1], 4000)
        path = sol.sol(ts)
        ax.plot(path[0], path[1], color="#ff7f0e", lw=1.4, alpha=0.95, zorder=5)
        ax.plot([x0[0]], [x0[1]], "o", color="#ff7f0e", ms=4, zorder=5)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("$x_1$", fontsize=13)
    ax.set_ylabel("$x_2$", fontsize=13)
    if title:
        ax.set_title(title, fontsize=12)
    return ax


def plot_alignment_margins(system: RampSystem, *, spec=PAPER, levels=(1, 2, 3), ax=None):
    """Histogram of the inward normal component over every oriented wall.

    Bars left of zero are walls the rectangular geometrization fails to align --
    the ones the R2/R3 construction must replace.
    """
    if ax is None:
        _fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for level in levels:
        report = check_alignment(system, spec=spec, level=level)
        values = [w.min_signed for w in report.walls]
        ax.hist(values, bins=45, histtype="step", lw=1.8, label=f"$\\mathcal{{F}}_{level}$")
    ax.axvline(0.0, color="k", lw=1.0, ls="--")
    ax.set_xlabel("inward normal component on the embedded wall")
    ax.set_ylabel("walls")
    ax.legend()
    return ax
