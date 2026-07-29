"""Generate the figure set into ``rookfields/figures/``."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .. import r2_manifolds as R2  # noqa: E402
from ..networks import INTRO_PERIODIC, VAN_DER_POL  # noqa: E402
from ..ramp import RampSystem  # noqa: E402
from ..spec import PAPER  # noqa: E402
from .phase2d import plot_alignment_margins, plot_phase_portrait_2d  # noqa: E402
from .phase3d import (
    plot_go_manifold_3d,
    plot_go_manifold_in_context,
    plot_morse_sets_3d,
)  # noqa: E402

FIGURES = Path(__file__).resolve().parents[3] / "figures"

#: Widths satisfying H_0, H_1, H_2, H_3; the published ones do not.
ADMISSIBLE_WIDTH = 0.005


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=FIGURES)
    parser.add_argument("--dpi", type=int, default=170)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []

    def save(fig, name):
        path = args.out / name
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
        print(f"  wrote {path.name}")

    # ---- planar: van der Pol ------------------------------------------
    print("van der Pol (planar)")
    for tag, spec in (
        ("published", VAN_DER_POL),
        ("admissible", VAN_DER_POL.with_uniform_width(ADMISSIBLE_WIDTH)),
    ):
        system = RampSystem(spec)
        fig, ax = plt.subplots(figsize=(7.5, 7.0))
        plot_phase_portrait_2d(
            system,
            ax=ax,
            level=3,
            trajectories=[[1.0, 1.0], [7.5, 8.0], [4.0, 7.0]],
            title=(
                f"van der Pol ramp system, $\\mathcal{{F}}_3$ Morse sets "
                f"({'published widths' if tag == 'published' else 'admissible widths'})\n"
                "blue: $\\dot x_1=0$, red: $\\dot x_2=0$, black: GO-manifolds, "
                "orange: trajectories"
            ),
        )
        save(fig, f"van_der_pol_phase_{tag}.png")

    system = RampSystem(VAN_DER_POL.with_uniform_width(ADMISSIBLE_WIDTH))
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    plot_alignment_margins(system, ax=ax)
    ax.set_title(
        "Inward normal component per oriented wall, rectangular geometrization\n"
        "$\\mathcal{F}_1$ is aligned; $\\mathcal{F}_2,\\mathcal{F}_3$ need the "
        "GO / cycle-surface replacement",
        fontsize=10,
    )
    save(fig, "van_der_pol_alignment_margins.png")

    # ---- three-dimensional: the periodic-orbit example -----------------
    print("intro periodic (3D)")
    system3 = RampSystem(INTRO_PERIODIC.with_uniform_width(ADMISSIBLE_WIDTH))

    fig = plt.figure(figsize=(8.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    plot_morse_sets_3d(
        system3,
        ax=ax,
        level=3,
        trajectories=[[1.0, 1.0, 1.0], [8.0, 8.0, 8.0]],
        title="ex:periodicOrbit3 - $\\mathcal{F}_3$ Morse sets in phase space,\n"
        "with two genuine ramp trajectories",
    )
    ax.view_init(elev=22, azim=-58)
    save(fig, "intro_periodic_morse_sets_3d.png")

    pairs, _ = R2.go_pairs(system3, spec=PAPER)
    external = [p for p in pairs if p.external]
    manifolds = [
        R2.build_go_manifold(system3, p, base_samples=9, steps=60)
        for p in external[:8]
    ]
    fig = plt.figure(figsize=(8.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    plot_go_manifold_3d(
        system3,
        manifolds,
        ax=ax,
        title="GO-manifolds of defn:GO-manifold for the 3D ramp system\n"
        "black: $\\delta$-nullcline base; sheet: flow-out under "
        "eq:rampSysPerturbed; dashed: terminal slice on the directed wall",
    )
    ax.view_init(elev=20, azim=-62)
    save(fig, "intro_periodic_go_manifolds_3d.png")

    # a single GO sheet, in the context of the cells whose face it replaces
    fig = plt.figure(figsize=(8.5, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    one = R2.build_go_manifold(system3, external[0], base_samples=15, steps=90)
    plot_go_manifold_in_context(
        system3, one, ax=ax,
        title=(
            f"GO-manifold for pair ({one.pair.face}, {one.pair.coface}), "
            f"GO-pair $(n_g,n_o)=({one.pair.n_grad},{one.pair.n_opaque})$, "
            f"$r_{{n_g}}={one.pair.r_grad:+d}$\n"
            f"$\\delta={one.delta:.4f}$, $\\varepsilon={one.epsilon:.4f}$; "
            "the sheet replaces the shared rectangular face of "
            r"$b(\xi)$ and $b(\xi')$"
        ),
    )
    ax.view_init(elev=16, azim=-72)
    save(fig, "go_manifold_single.png")

    print(f"\n{len(written)} figures in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
