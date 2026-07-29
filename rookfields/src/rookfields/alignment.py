"""Numerical alignment of a geometrization with the ramp vector field.

``defn:aligned`` (geometrizationv2.tex:382-392) asks the flow to point inward at
every ``(N-1)``-cell of ``bdy(N)`` for every ``N`` in the lattice, and
``thm:inward``/``cor:inward`` (:393-433) reduce that to the sign of
``<F(x), z_zeta(x)>``.

For a *rectangular* geometrization the walls of ``X_b`` are axis-parallel, so
``z_zeta`` is a coordinate direction and the test is a sign condition on one
component of ``F``.  ``thm:R1ABlattice`` (GlobalDynR1v3.tex:15-21) asserts that
every rectangular geometrization is aligned over ``N(F_1)``; this module checks
that assertion pointwise.

Checking all lattice elements at once.  ``N(F_i)`` consists of the down-sets of
``SCC(F_i)`` pulled back through ``pi_b``.  A wall between adjacent blowup top
cells is on the boundary of *some* element exactly when the two cells carry
different SCC values, and the required crossing direction is then the one the
multivalued map prescribes.  So it suffices to test, for every wall the map
orients, that the field crosses it that way -- which covers every ``N``
simultaneously.
"""

from __future__ import annotations

import dataclasses
import itertools

import numpy as np

from .blowup import SpecCubicalBlowupGraph
from .geometrization import RectangularGeometrization
from .ramp import RampSystem
from .spec import PAPER, Spec


@dataclasses.dataclass
class WallCheck:
    """One oriented wall of ``X_b`` and the field's behaviour on it."""

    source: int
    target: int
    direction: int
    #: +1 when the map sends the cell with the larger n-coordinate to the smaller.
    required_sign: int
    min_signed: float
    max_signed: float
    samples: int
    location: tuple[float, ...]

    @property
    def ok(self) -> bool:
        return self.min_signed > 0.0


@dataclasses.dataclass
class AlignmentReport:
    level: int
    spec: str
    walls: list[WallCheck]
    skipped_double_edges: int
    samples_per_wall: int

    @property
    def ok(self) -> bool:
        return all(w.ok for w in self.walls)

    @property
    def failures(self) -> list[WallCheck]:
        return [w for w in self.walls if not w.ok]

    @property
    def margin(self) -> float:
        """The worst signed normal component over every oriented wall."""
        return min((w.min_signed for w in self.walls), default=float("inf"))

    def __str__(self) -> str:  # pragma: no cover - display only
        status = "aligned" if self.ok else f"NOT aligned ({len(self.failures)} walls)"
        return (
            f"F_{self.level} [{self.spec}]: {status}; "
            f"{len(self.walls)} oriented walls checked, "
            f"{self.skipped_double_edges} double-edge walls skipped, "
            f"worst inward component {self.margin:.6g}"
        )


def _sample_face(intervals, fixed_axis: int, fixed_value: float, per_axis: int):
    """Interior sample points of an axis-parallel face."""
    axes = []
    for n, (lo, hi) in enumerate(intervals):
        if n == fixed_axis:
            axes.append([fixed_value])
        else:
            # strictly interior, to stay off the ramp-window corners where the
            # field is only piecewise smooth
            axes.append(
                [lo + (hi - lo) * (i + 1) / (per_axis + 1) for i in range(per_axis)]
            )
    return [np.array(p, dtype=float) for p in itertools.product(*axes)]


def check_alignment(
    system: RampSystem,
    *,
    spec: Spec = PAPER,
    level: int = 1,
    samples_per_axis: int = 3,
    outer: str = "global_bound",
) -> AlignmentReport:
    """Test ``recG`` against the ramp field over every element of ``N(F_level)``."""
    geo = RectangularGeometrization(system)
    labelling, num_thresholds = system.wall_labelling(outer=outer)
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=spec, level=level
    )

    bc = stg.blowup_complex
    dim = stg.dim
    edges = set(stg.digraph.edges())
    walls: list[WallCheck] = []
    skipped = 0

    for cell in bc(dim):
        if bc.rightfringe(cell):
            continue
        coords = bc.coordinates(cell)
        for n in range(dim):
            neighbour = cell + stg.blowup_jump[n]
            if bc.rightfringe(neighbour):
                continue
            forward = (cell, neighbour) in edges
            backward = (neighbour, cell) in edges
            if forward and backward:
                skipped += 1
                continue
            if not forward and not backward:
                # The definition can leave a pair with no edge; nothing to align.
                continue

            # required_sign = +1 means F_n > 0 on the shared face (flow moves
            # from `cell` into `neighbour`, i.e. in the +n direction)
            required_sign = 1 if forward else -1

            intervals = geo.rectangle(list(coords), [1] * dim)
            face_value = geo.embed_coordinate(n, coords[n] + 1)
            points = _sample_face(intervals, n, face_value, samples_per_axis)
            signed = [required_sign * float(system.vector_field(p)[n]) for p in points]
            walls.append(
                WallCheck(
                    source=cell,
                    target=neighbour,
                    direction=n,
                    required_sign=required_sign,
                    min_signed=min(signed),
                    max_signed=max(signed),
                    samples=len(points),
                    location=tuple(float(v) for v in points[0]),
                )
            )

    return AlignmentReport(
        level=level,
        spec=spec.name,
        walls=walls,
        skipped_double_edges=skipped,
        samples_per_wall=samples_per_axis ** (dim - 1),
    )


def main(argv=None) -> int:
    import argparse
    import json
    from pathlib import Path

    from .networks import RAMP_SYSTEMS
    from .spec import LEGACY

    reports = Path(__file__).resolve().parents[2] / "reports"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--out", type=Path, default=reports / "alignment.md")
    parser.add_argument("--json", type=Path, default=reports / "alignment.json")
    args = parser.parse_args(argv)

    widths = [None, 0.02, 0.005, 0.002]
    out: list[str] = []
    data: list[dict] = []

    out.append("# Numerical alignment of the rectangular geometrization\n")
    out.append(
        "Generated by `python -m rookfields.alignment`. For each oriented wall "
        "of `X_b` the inward component of the ramp field is sampled on the "
        "embedded face; `thm:R1ABlattice` asserts every such component is "
        "positive. Walls carrying a double edge lie inside a recurrent set and "
        "carry no alignment requirement, so they are skipped.\n"
    )
    out.append(
        "The `worst inward component` column is the minimum of "
        "`required_sign * F_n(x)` over every sampled point of every oriented "
        "wall: positive means aligned, with that much margin.\n"
    )
    out.append("| system | uniform half-width | level | spec | oriented walls | skipped | worst inward component | aligned |")
    out.append("|---|---|---:|---|---:|---:|---:|:-:|")

    for name, spec in RAMP_SYSTEMS.items():
        for h in widths:
            sp = spec if h is None else spec.with_uniform_width(h)
            system = RampSystem(sp)
            for level in (1, 2, 3):
                for algo in (LEGACY, PAPER):
                    report = check_alignment(
                        system, spec=algo, level=level, samples_per_axis=args.samples
                    )
                    tag = "published" if h is None else f"{h:g}"
                    out.append(
                        f"| `{name}` | {tag} | {level} | {algo.name} | "
                        f"{len(report.walls)} | {report.skipped_double_edges} | "
                        f"{report.margin:.5g} | {'ok' if report.ok else '**no**'} |"
                    )
                    data.append(
                        {
                            "system": name,
                            "width": tag,
                            "level": level,
                            "spec": algo.name,
                            "walls": len(report.walls),
                            "skipped": report.skipped_double_edges,
                            "margin": report.margin,
                            "aligned": report.ok,
                            "failures": [
                                {
                                    "source": w.source,
                                    "target": w.target,
                                    "direction": w.direction,
                                    "min_signed": w.min_signed,
                                    "location": w.location,
                                }
                                for w in report.failures[:10]
                            ],
                        }
                    )

    out.append(
        "\n## Reading\n\n"
        "A positive worst-case component means the rectangular geometrization is "
        "aligned with the ramp field over every element of `N(F_level)`, which is "
        "the numerical content of `thm:R1ABlattice`. A negative one identifies a "
        "wall the combinatorial map orients one way and the field crosses the "
        "other -- reported per wall in the JSON.\n"
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out) + "\n")
    args.json.write_text(json.dumps(data, indent=2, default=str))
    print(f"wrote {args.out}\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
