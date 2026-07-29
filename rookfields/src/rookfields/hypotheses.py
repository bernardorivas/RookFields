"""The admissibility and width conditions Lambda(R), H_0, H_1, H_2, H_3.

Nothing in the project verified these before.  The manuscript itself notes that
``ex:ramp_van_der_pol`` "gives numerical ramp parameters and computes F_3, but
it does not verify h in H_3" -- this module closes that gap.

Sources (``RampSystemsv4.tex``):

* ``defn:admissiblehtheta``      :281-312   ordered, non-overlapping windows
* ``defn:admissible-parameters`` :347-362   Lambda(R): gamma_n theta != E_n(D_v)
* ``defn:H1``                    :493-522
* ``defn:H2``                    :559-590
* ``defn:H3``                    :617-622
* ``prop:smallh``                :624-664   constructive nonemptiness bounds

H1 note.  The literal definition uses the half-open windows
``(theta, theta+h]`` and ``[theta-h, theta)``, which permits a focal value
exactly at an endpoint.  ``issues/source-defects.md`` records that every later
use depends on the intended *strict* convention, so ``H1`` is checked with the
endpoint excluded and ``strict=False`` reproduces the literal reading.
"""

from __future__ import annotations

import dataclasses
import itertools
import math

from .blowup import SpecCubicalBlowupGraph
from .ramp import RampSystem
from .spec import PAPER, Spec


@dataclasses.dataclass
class Verdict:
    """Outcome of one hypothesis check."""

    name: str
    holds: bool
    #: Human-readable violations, empty when the hypothesis holds.
    violations: list[str] = dataclasses.field(default_factory=list)
    #: Quantitative margin: positive means satisfied with room to spare.
    margin: float | None = None
    #: How many instances the quantifier ranged over.
    checked: int = 0

    def __bool__(self) -> bool:
        return self.holds

    def __str__(self) -> str:
        status = "holds" if self.holds else f"FAILS ({len(self.violations)})"
        margin = "" if self.margin is None else f", margin {self.margin:.4g}"
        return f"{self.name}: {status} over {self.checked} instances{margin}"


def _boxes(system: RampSystem):
    return itertools.product(*[range(k + 1) for k in system.num_thresholds])


def _interior_boxes(system: RampSystem):
    """``v in prod {1,...,K(n)}``, the range H1 and Lambda(R) quantify over."""
    return itertools.product(*[range(1, k + 1) for k in system.num_thresholds])


def check_admissible_theta_h(system: RampSystem) -> Verdict:
    """``defn:admissiblehtheta``: distinct thresholds, non-overlapping windows."""
    violations = []
    margin = math.inf
    checked = 0
    for n in range(system.dim):
        thetas = system.sorted_thresholds[n]
        widths = system.sorted_widths[n]
        for k in range(len(thetas) - 1):
            checked += 1
            gap = (thetas[k + 1] - widths[k + 1]) - (thetas[k] + widths[k])
            margin = min(margin, gap)
            if gap <= 0:
                violations.append(
                    f"coordinate {n}: ramp windows {k+1} and {k+2} overlap "
                    f"({thetas[k]}+{widths[k]} >= {thetas[k+1]}-{widths[k+1]})"
                )
        if len(set(thetas)) != len(thetas):
            violations.append(f"coordinate {n}: repeated threshold in {thetas}")
    return Verdict("H_0 / admissible (theta, h)", not violations, violations,
                   None if margin is math.inf else margin, checked)


def check_lambda_R(system: RampSystem) -> Verdict:
    """``defn:admissible-parameters``: no focal value sits on a bounding threshold.

    Quantified over ``v in prod {0, ..., K(n)}`` and every ``n``, against the
    box's own two bounding thresholds ``theta_{m_{k_n}}`` and
    ``theta_{m_{k_n+1}}`` (with the sentinels ``theta_{m_0} = 0`` and
    ``theta_{m_{K+1}} = GB_n``).  This is what makes
    ``sgn(-gamma_n theta + E_n(D_v))`` -- and hence the wall labelling -- well
    defined.
    """
    violations = []
    margin = math.inf
    checked = 0
    for box in _boxes(system):
        e = system.production_on_box(box)
        for n in range(system.dim):
            k = box[n]
            for offset in (0, 1):
                checked += 1
                gap = abs(-system.gamma[n] * system.theta(n, k + offset) + float(e[n]))
                margin = min(margin, gap)
                if gap == 0.0:
                    violations.append(
                        f"box {box}, coordinate {n}: "
                        f"gamma_{n} * theta_{{m_{k + offset}}} = E_{n}(D_v) = {float(e[n]):.6g}"
                    )
    return Verdict("Lambda(R)", not violations, violations,
                   None if margin is math.inf else margin, checked)


def check_H1(system: RampSystem, *, strict: bool = True) -> Verdict:
    """``defn:H1``: focal values stay clear of the ramp half-windows.

    With ``strict=True`` the endpoints are excluded too, which is the intended
    convention the later chapters rely on.
    """
    violations = []
    margin = math.inf
    checked = 0
    for box in _interior_boxes(system):
        e = system.production_on_box(box)
        for n in range(system.dim):
            focal = float(e[n]) / system.gamma[n]
            k = box[n]
            lo_t, lo_h = system.theta(n, k), system.width(n, k)
            hi_t, hi_h = system.theta(n, k + 1), system.width(n, k + 1)
            checked += 1
            # forbidden: (theta_k, theta_k + h_k]  and  [theta_{k+1} - h_{k+1}, theta_{k+1})
            right_gap = min(abs(focal - lo_t), abs(focal - (lo_t + lo_h)))
            left_gap = min(abs(focal - (hi_t - hi_h)), abs(focal - hi_t))
            in_right = lo_t < focal <= lo_t + lo_h
            in_left = hi_t - hi_h <= focal < hi_t
            if strict:
                in_right = in_right or focal == lo_t
                in_left = in_left or focal == hi_t
            if in_right or in_left:
                violations.append(
                    f"box {box}, coordinate {n}: focal value {focal:.6g} lies in a "
                    f"ramp half-window of [{lo_t}, {hi_t}]"
                )
            else:
                margin = min(margin, right_gap, left_gap)
    return Verdict("H_1", not violations, violations,
                   None if margin is math.inf else margin, checked)


def check_H2(system: RampSystem, *, spec: Spec = PAPER) -> Verdict:
    """``defn:H2``: the external GO crossing estimate.

    Quantified over ``(xi, xi') in D(Phi) cap (X^(N-2) x X^(N-1))`` with GO-pair
    ``(n_g, n_o)`` and ``xi' not in F_2(xi)`` -- the *externally pruned* pairs
    only.  Internally pruned pairs use the short internal flow-out instead and
    carry no width bound.
    """
    labelling, num_thresholds = system.wall_labelling()
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=spec, level=2
    )
    cc = stg.cubical_complex
    dim = stg.dim

    violations = []
    margin = math.inf
    checked = 0
    external = 0

    for cc_face, cc_coface, n_grad, n_opaque in stg.indecisive_drift_pairs(
        codim=(dim - 2, dim - 1)
    ):
        # externally pruned: xi' not in F_2(xi)
        if stg.has_edge(cc_face, cc_coface):
            continue
        external += 1
        checked += 1

        coords = list(cc.coordinates(cc_coface))
        shape = cc.cell_shape(cc_coface)
        bits = [1 if shape & (1 << n) else 0 for n in range(dim)]

        L = system.L(coords, bits, n_grad)
        U = system.U(coords, bits, n_opaque)
        Xi = system.interval_length(coords, bits, n_opaque)
        if U == 0.0:
            violations.append(f"pair ({cc_face},{cc_coface}): U_{n_opaque} = 0")
            continue
        lhs = 2.0 * system.width(n_grad, coords[n_grad])
        rhs = (L / U) * (Xi / 2.0)
        margin = min(margin, rhs - lhs)
        if not lhs < rhs:
            violations.append(
                f"pair ({cc_face},{cc_coface}) GO=({n_grad},{n_opaque}): "
                f"2h = {lhs:.6g} but (L/U)(Xi/2) = {rhs:.6g}"
            )

    verdict = Verdict("H_2", not violations, violations,
                      None if margin is math.inf else margin, checked)
    verdict.violations = violations
    if external == 0:
        verdict.violations = []
        verdict.holds = True
    return verdict


def check_H3(system: RampSystem, *, spec: Spec = PAPER) -> Verdict:
    """``defn:H3``: ``H_3 = H_1`` in dimension two; the 3-cycle product bound in three.

    The denominator is
    ``-g1 g2 g3 + (g1+g2+g3)(g1 g2 + g1 g3 + g2 g3)``.
    """
    if system.dim == 2:
        v = check_H1(system)
        return Verdict("H_3 (= H_1 in dimension 2)", v.holds, v.violations, v.margin, v.checked)
    if system.dim != 3:
        return Verdict(
            "H_3",
            False,
            ["defn:H3 is only defined for N = 2 and N = 3"],
            None,
            0,
        )

    base = check_H1(system)
    if not base.holds:
        return Verdict("H_3", False, ["H_1 fails: " + v for v in base.violations], base.margin, base.checked)

    labelling, num_thresholds = system.wall_labelling()
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=spec, level=3
    )
    cc = stg.cubical_complex
    g = system.gamma
    denominator = -g[0] * g[1] * g[2] + (g[0] + g[1] + g[2]) * (
        g[0] * g[1] + g[0] * g[2] + g[1] * g[2]
    )

    violations = []
    margin = math.inf
    checked = 0

    for cell, cycle in stg.three_cycle_vertices():
        coords = list(cc.coordinates(cell))
        if any(c < 1 or c > system.num_thresholds[n] for n, c in enumerate(coords)):
            continue  # kappa^-(xi) would leave the complex
        checked += 1
        reg = stg.active_regulation_map(cell)
        kappa_plus = tuple(coords)
        kappa_minus = tuple(c - 1 for c in coords)
        e_plus = system.production_on_box(kappa_plus)
        e_minus = system.production_on_box(kappa_minus)

        lhs = 8.0
        for n in range(3):
            lhs *= system.width(n, coords[n])
        numerator = 1.0
        for n in range(3):
            numerator *= abs(float(e_plus[n]) - float(e_minus[n]))
        rhs = numerator / denominator
        margin = min(margin, rhs - lhs)
        if not lhs < rhs:
            violations.append(
                f"vertex {cell} coords {coords} cycle {cycle} (o_xi={reg}): "
                f"8*prod h = {lhs:.6g} but bound = {rhs:.6g}"
            )

    return Verdict("H_3", not violations, violations,
                   None if margin is math.inf else margin, checked)


# ---------------------------------------------------------------------------
# prop:smallh -- constructive bounds giving a width that works
# ---------------------------------------------------------------------------


def uniform_width_bound(system: RampSystem, level: int) -> float:
    """``h_0(gamma, nu, theta, i)`` from the proof of ``prop:smallh``.

    Returns a ``h_bar`` such that the uniform width ``(h_bar, ..., h_bar)`` lies
    in ``H_i``.  Implemented for ``i in {0, 1, 3}``, which is what the
    proposition states; the proof also covers ``i = 2`` but its constant ``c``
    is built from ``L/U``, which themselves depend on ``h`` -- a circularity
    the manuscript records as an open obligation, so it is not implemented here.
    """
    if level not in (0, 1, 3):
        raise ValueError(
            "prop:smallh states nonemptiness for i in {0,1,3}; the i=2 bound in "
            "its proof is circular in h and is tracked as an open obligation"
        )

    # h_0(., 0): half the smallest gap between consecutive thresholds
    h0 = math.inf
    for n in range(system.dim):
        thetas = system.sorted_thresholds[n]
        for k in range(len(thetas) - 1):
            h0 = min(h0, (thetas[k + 1] - thetas[k]) / 2.0)
    if h0 is math.inf:
        h0 = 1.0
    if level == 0:
        return h0

    # h~: distance from every focal value to the enclosing thresholds
    h_tilde = math.inf
    for box in _interior_boxes(system):
        e = system.production_on_box(box)
        for n in range(system.dim):
            focal = float(e[n]) / system.gamma[n]
            k = box[n]
            lo, hi = system.theta(n, k), system.theta(n, k + 1)
            if lo < focal < hi:
                h_tilde = min(h_tilde, abs(lo - focal), abs(hi - focal))
    h1 = min(h0, h_tilde if h_tilde is not math.inf else h0)
    if level == 1:
        return h1

    # H_3 in dimension three: the cube-root bound
    if system.dim == 2:
        return h1
    g = system.gamma
    denominator = -g[0] * g[1] * g[2] + (g[0] + g[1] + g[2]) * (
        g[0] * g[1] + g[0] * g[2] + g[1] * g[2]
    )
    labelling, num_thresholds = system.wall_labelling()
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=PAPER, level=3
    )
    cc = stg.cubical_complex
    best = math.inf
    for cell, _cycle in stg.three_cycle_vertices():
        coords = list(cc.coordinates(cell))
        if any(c < 1 or c > system.num_thresholds[n] for n, c in enumerate(coords)):
            continue
        e_plus = system.production_on_box(tuple(coords))
        e_minus = system.production_on_box(tuple(c - 1 for c in coords))
        numerator = 1.0
        for n in range(3):
            numerator *= abs(float(e_plus[n]) - float(e_minus[n]))
        best = min(best, 0.5 * (numerator / denominator) ** (1.0 / 3.0))
    return min(h1, best) if best is not math.inf else h1


def report(system: RampSystem, *, spec: Spec = PAPER) -> list[Verdict]:
    """Run every applicable check for one ramp system."""
    verdicts = [
        check_admissible_theta_h(system),
        check_lambda_R(system),
        check_H1(system),
    ]
    verdicts.append(check_H2(system, spec=spec))
    if system.dim in (2, 3):
        verdicts.append(check_H3(system, spec=spec))
    return verdicts


def largest_admissible_width(
    system_spec, *, target: int = 2, lo: float = 1e-5, hi: float = 1.0, steps: int = 40
) -> float | None:
    """Bisect for the largest uniform width satisfying ``H_0, H_1, H_{target}``.

    The width classes are nested and shrink as constraints are added, and every
    condition is of the form "``h`` small enough", so bisection is sound.
    Returns ``None`` when even ``lo`` fails.
    """
    from .ramp import RampSystem as _RS

    def ok(h: float) -> bool:
        s = _RS(system_spec.with_uniform_width(h))
        if not check_admissible_theta_h(s).holds:
            return False
        if not check_H1(s).holds:
            return False
        if target == 2:
            return check_H2(s).holds
        if target == 3:
            return check_H3(s).holds
        return True

    if not ok(lo):
        return None
    if ok(hi):
        return hi
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def main(argv=None) -> int:
    import argparse
    import json
    from pathlib import Path

    from .networks import RAMP_SYSTEMS
    from .ramp import RampSystem as _RS
    from .wall_labeling import WallLabeling

    reports = Path(__file__).resolve().parents[2] / "reports"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=reports / "hypotheses.md")
    parser.add_argument("--json", type=Path, default=reports / "hypotheses.json")
    args = parser.parse_args(argv)

    widths = [None, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    out: list[str] = []
    data: dict = {}

    out.append("# Analytic admissibility of the monograph's ramp systems\n")
    out.append(
        "Generated by `python -m rookfields.hypotheses`. Two systems in the "
        "monograph are specified directly by ramp parameters rather than through "
        "DSGRN, and the text never checks that those parameters satisfy the width "
        "conditions the theorems require.\n"
    )
    out.append(
        "The wall labelling used here follows `defn:ramp-wall-labeling` "
        "(RampSystemsv4.tex:375), whose outermost sentinel threshold is `GB_n`. "
        "`DSGRN_utils.WallLabelling.ramp_system_wall_labelling` (:44) instead uses "
        "`theta_{m_K} + 10 h_{m_K}`; see the section on that below.\n"
    )

    for name, spec in RAMP_SYSTEMS.items():
        base = _RS(spec)
        out.append(f"\n## `{name}`\n")
        out.append(f"Source: `{spec.source}`. Dimension {base.dim}, K(n) = {base.num_thresholds}.\n")
        out.append(f"Global bound GB = `{[round(float(g), 4) for g in base.global_bound]}`.\n")
        out.append("| uniform half-width | H_0 | Lambda(R) | H_1 | H_2 | H_3 |")
        out.append("|---|:-:|:-:|:-:|:-:|:-:|")
        rows = []
        for h in widths:
            sp = spec if h is None else spec.with_uniform_width(h)
            s = _RS(sp)
            v = {
                "H_0": check_admissible_theta_h(s),
                "Lambda(R)": check_lambda_R(s),
                "H_1": check_H1(s),
                "H_2": check_H2(s),
                "H_3": check_H3(s),
            }
            tag = "published" if h is None else f"{h:g}"
            mark = lambda k: "ok" if v[k].holds else "**FAILS**"
            out.append(
                f"| {tag} | {mark('H_0')} | {mark('Lambda(R)')} | {mark('H_1')} | "
                f"{mark('H_2')} | {mark('H_3')} |"
            )
            rows.append(
                {
                    "width": tag,
                    **{k: {"holds": x.holds, "checked": x.checked, "margin": x.margin,
                           "violations": x.violations[:5]} for k, x in v.items()},
                }
            )
        data[name] = {"rows": rows}

        h2 = largest_admissible_width(spec, target=2)
        h3 = largest_admissible_width(spec, target=3)
        data[name]["largest_width_H2"] = h2
        data[name]["largest_width_H3"] = h3
        out.append(
            f"\nLargest uniform half-width satisfying H_0+H_1+H_2: "
            f"`{h2:.5g}`; with H_3 instead of H_2: `{h3:.5g}`. "
            f"`prop:smallh`'s constructive H_1 bound is `{uniform_width_bound(base, 1):.5g}`.\n"
        )

        published = check_H1(base)
        if not published.holds:
            out.append("\nViolations at the published widths:\n")
            for m in published.violations[:6]:
                out.append(f"- {m}")
            if len(published.violations) > 6:
                out.append(f"- ... and {len(published.violations) - 6} more")

    # -- the WallLabelling.py sentinel ---------------------------------
    out.append("\n## The outermost sentinel threshold\n")
    out.append(
        "`defn:ramp-wall-labeling` sets `theta_{m_{K(n)+1}} = GB_n`. "
        "`WallLabelling.py:44` uses `theta_{m_K} + 10 h_{m_K}` instead. That "
        "surrogate depends on `h`, so it contradicts `prop:wall-labeling-const` "
        "(RampSystemsv4.tex:490: the labelling is constant in `h` on `H_0`), and "
        "it can break the standing strong-dissipativity assumption once `h` is "
        "small -- which is exactly the regime `H_1` forces.\n"
    )
    out.append("| system | width | sentinel | `def:wall_labeling` valid | strongly dissipative |")
    out.append("|---|---|---|:-:|:-:|")
    sentinel_rows = []
    for name, spec in RAMP_SYSTEMS.items():
        for h in [None, 0.1, 0.05, 0.02, 0.005]:
            sp = spec if h is None else spec.with_uniform_width(h)
            s = _RS(sp)
            for outer in ("legacy", "global_bound"):
                try:
                    lab, K = s.wall_labelling(outer=outer)
                except ValueError as exc:
                    out.append(
                        f"| `{name}` | {'published' if h is None else h} | "
                        f"`{outer}` | error | {exc} |"
                    )
                    continue
                wl = WallLabeling(lab, K)
                valid, dissipative = wl.is_valid(), wl.is_strongly_dissipative()
                out.append(
                    f"| `{name}` | {'published' if h is None else h} | `{outer}` | "
                    f"{'ok' if valid else '**no**'} | {'ok' if dissipative else '**no**'} |"
                )
                sentinel_rows.append(
                    {"system": name, "width": h, "sentinel": outer,
                     "valid": valid, "strongly_dissipative": dissipative}
                )
    data["sentinel"] = sentinel_rows

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out) + "\n")
    args.json.write_text(json.dumps(data, indent=2, default=str))
    print(f"wrote {args.out}\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
