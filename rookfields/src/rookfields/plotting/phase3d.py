"""Three-dimensional views: Morse sets and GO-manifold sheets in phase space.

This fills TODO 4 of ``CubicalBlowupGraph.py:7`` ("Add 3D Morse sets plotting
capabilities"), and draws in real ramp coordinates so the swept sheets, the
trajectories, and the cells can be compared directly.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from ..geometrization import RectangularGeometrization  # noqa: E402
from ..pipeline import conley_morse_graph  # noqa: E402
from ..ramp import RampSystem  # noqa: E402
from ..spec import PAPER  # noqa: E402
from .phase2d import MORSE_COLOURS  # noqa: E402


def _box_faces(intervals):
    """The six faces of an axis-aligned box, as vertex lists."""
    (x0, x1), (y0, y1), (z0, z1) = intervals
    return [
        [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
        [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
        [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],
    ]


def plot_morse_sets_3d(
    system: RampSystem,
    *,
    spec=PAPER,
    level: int = 3,
    morse_nodes=None,
    ax=None,
    alpha: float = 0.35,
    trajectories: list | None = None,
    title: str | None = None,
):
    """Draw the Morse sets as embedded blowup cells in ``R^3``."""
    assert system.dim == 3, "this view needs a three-dimensional system"
    geo = RectangularGeometrization(system)
    labelling, K = system.wall_labelling()
    result = conley_morse_graph(
        labelling=labelling, num_thresholds=K, spec=spec, level=level
    )
    stg = result.stg

    if ax is None:
        fig = plt.figure(figsize=(8.5, 7.5))
        ax = fig.add_subplot(111, projection="3d")

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
        if morse_nodes is not None and node not in morse_nodes:
            continue
        coords = stg.blowup_complex.coordinates(cell)
        intervals = geo.rectangle(list(coords), [1, 1, 1])
        colour = MORSE_COLOURS[node % len(MORSE_COLOURS)]
        ax.add_collection3d(
            Poly3DCollection(
                _box_faces(intervals),
                facecolors=colour,
                edgecolors="0.25",
                linewidths=0.3,
                alpha=alpha,
            )
        )

    for x0 in trajectories or []:
        sol = system.flow(x0, t_span=(0.0, 80.0))
        ts = np.linspace(0.0, sol.t[-1], 6000)
        path = sol.sol(ts)
        ax.plot(path[0], path[1], path[2], color="#ff7f0e", lw=1.1, alpha=0.95)

    gb = system.global_bound
    ax.set_xlim(0, gb[0])
    ax.set_ylim(0, gb[1])
    ax.set_zlim(0, gb[2])
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_zlabel("$x_3$")
    if title:
        ax.set_title(title, fontsize=12)
    return ax


def plot_go_manifold_in_context(
    system: RampSystem,
    manifold,
    *,
    ax=None,
    pad: float = 0.35,
    title: str | None = None,
):
    """One GO sheet together with the two blowup cells whose face it replaces.

    ``b(xi)`` is the ramp-window cell and ``b(xi')`` the adjacent plateau cell;
    their shared rectangular face in the ``n_o`` direction is what
    ``defn:R2-face-coverage`` swaps for this sheet.  Drawn with the sweep
    direction stretched, since the sheet lives inside a ramp window of
    half-width ``h`` and would otherwise be a hairline.
    """
    from ..geometrization import RectangularGeometrization

    geo = RectangularGeometrization(system)
    pair = manifold.pair

    if ax is None:
        fig = plt.figure(figsize=(8.5, 7.0))
        ax = fig.add_subplot(111, projection="3d")

    for coords, shape, colour, label in (
        (pair.face_coords, pair.face_shape, "#9467bd", r"$b(\xi)$"),
        (pair.coface_coords, pair.coface_shape, "#2ca02c", r"$b(\xi')$"),
    ):
        blowup = [2 * v + w for v, w in zip(coords, shape)]
        intervals = geo.rectangle(blowup, [1] * system.dim)
        ax.add_collection3d(
            Poly3DCollection(
                _box_faces(intervals),
                facecolors=colour,
                edgecolors="0.3",
                linewidths=0.6,
                alpha=0.12,
            )
        )
        centre = [0.5 * (a + b) for a, b in intervals]
        ax.text(*centre, label, color=colour, fontsize=11)

    trajs = [t for t in manifold.trajectories if len(t) > 1]
    steps = min(len(t) for t in trajs)
    quads = []
    for a, b in zip(trajs, trajs[1:]):
        for i in range(steps - 1):
            quads.append([a[i], a[i + 1], b[i + 1], b[i]])
    ax.add_collection3d(
        Poly3DCollection(quads, facecolors="#1f77b4", edgecolors="none", alpha=0.8)
    )
    ax.plot(
        manifold.base[:, 0], manifold.base[:, 1], manifold.base[:, 2],
        color="#111111", lw=2.6, label=r"$\mathrm{Nul}^\delta_{n_o}$ (base)",
    )
    ax.plot(
        manifold.terminal[:, 0], manifold.terminal[:, 1], manifold.terminal[:, 2],
        color="#d62728", lw=2.2, ls="--", label="terminal slice on the directed wall",
    )

    lows = manifold.base.min(axis=0)
    highs = manifold.base.max(axis=0)
    for arr in (manifold.terminal,):
        lows = np.minimum(lows, arr.min(axis=0))
        highs = np.maximum(highs, arr.max(axis=0))
    for n, setter in enumerate((ax.set_xlim, ax.set_ylim, ax.set_zlim)):
        span = max(highs[n] - lows[n], 1e-3)
        setter(lows[n] - pad * span, highs[n] + pad * span)

    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_zlabel("$x_3$")
    ax.legend(loc="upper left", fontsize=9)
    if title:
        ax.set_title(title, fontsize=11)
    return ax


def plot_go_manifold_3d(
    system: RampSystem,
    manifolds,
    *,
    ax=None,
    show_base: bool = True,
    show_terminal: bool = True,
    sheet_colour: str = "#1f77b4",
    alpha: float = 0.55,
    title: str | None = None,
):
    """Render GO-manifolds as swept surfaces in ``R^3``.

    Each sheet is the union of the trajectories of ``eq:rampSysPerturbed``
    starting on the ``delta``-nullcline base, stopped at the directed
    ``n_g``-wall: the object ``defn:GO-manifold`` constructs to replace a
    rectangular ``(N-1)``-face.
    """
    assert system.dim == 3

    if ax is None:
        fig = plt.figure(figsize=(8.5, 7.5))
        ax = fig.add_subplot(111, projection="3d")

    for manifold in manifolds:
        trajs = [t for t in manifold.trajectories if len(t) > 1]
        if len(trajs) < 2:
            for t in trajs:
                ax.plot(t[:, 0], t[:, 1], t[:, 2], color=sheet_colour, lw=2.0)
            continue
        # quads between consecutive trajectories form the swept sheet
        quads = []
        steps = min(len(t) for t in trajs)
        for a, b in zip(trajs, trajs[1:]):
            for i in range(steps - 1):
                quads.append([a[i], a[i + 1], b[i + 1], b[i]])
        ax.add_collection3d(
            Poly3DCollection(
                quads,
                facecolors=sheet_colour,
                edgecolors="none",
                alpha=alpha,
            )
        )
        if show_base and manifold.base.size:
            ax.plot(
                manifold.base[:, 0], manifold.base[:, 1], manifold.base[:, 2],
                color="#111111", lw=2.4,
            )
        if show_terminal and manifold.terminal.size:
            ax.plot(
                manifold.terminal[:, 0], manifold.terminal[:, 1], manifold.terminal[:, 2],
                color="#d62728", lw=2.0, ls="--",
            )

    gb = system.global_bound
    ax.set_xlim(0, gb[0])
    ax.set_ylim(0, gb[1])
    ax.set_zlim(0, gb[2])
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_zlabel("$x_3$")
    if title:
        ax.set_title(title, fontsize=12)
    return ax
