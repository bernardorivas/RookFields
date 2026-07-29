"""R3: cycle geometry, and the spectral claims the chapter rests on.

``GlobalDynR3v5.tex`` proposes, for a semi-opaque cell whose regulation map has
a nontrivial cycle:

* ``N = 2``: an internal quadratic/elliptical surface ``M^i_eps(xi)`` for
  negative feedback, an external swept surface ``M^e_eps(xi)`` for positive
  (``defn:compatible-recG-F3N2``, :462-478);
* ``N = 3``: stable/unstable level sets ``L_s^{-1}(eps)``, ``L_u^{-1}(eps)`` for
  3-cycles, and projected families for 2-cycles (``defn:F3N3-Geo``, :499-520).

Two of its statements are recorded in ``issues/source-defects.md`` as false.
Both are cheap to settle numerically and are settled here:

1. "The positive 3-cycle H3 bound does not force a complex stable pair.  The
   explicit parameters ``(gamma_1,gamma_2,gamma_3) = (1,2,100)`` and signed
   slope product ``P = 31000`` satisfy the bound but yield three real
   eigenvalues."
2. "The vector used as a purported transverse vector of the swept surface is one
   of its generating tangent directions."
"""

from __future__ import annotations

import dataclasses
import itertools

import numpy as np

from .blowup import SpecCubicalBlowupGraph
from .ramp import RampSystem
from .spec import PAPER, Spec


# ---------------------------------------------------------------------------
# the cyclic-feedback linearisation
# ---------------------------------------------------------------------------


def cyclic_feedback_matrix(gamma, slopes, cycle) -> np.ndarray:
    """The Jacobian of a cyclic feedback system at an equilibrium.

    For a cycle ``sigma`` on ``{0..N-1}`` the linearisation is

        J = -diag(gamma) + S,     S[sigma(n), n] = slopes[n],

    so its characteristic equation is ``prod_n (lam + gamma_n) = P`` with
    ``P = prod_n slopes[n]`` the signed slope product (``eq:F3-bounds`` bounds
    exactly this product).
    """
    n = len(gamma)
    J = -np.diag(np.asarray(gamma, dtype=float))
    for k, source in enumerate(cycle):
        target = cycle[(k + 1) % len(cycle)]
        J[target, source] += slopes[k]
    return J


def eigenvalues_from_product(gamma, P: float) -> np.ndarray:
    """Roots of ``prod_n (lam + gamma_n) = P``."""
    poly = np.array([1.0])
    for g in gamma:
        poly = np.convolve(poly, [1.0, float(g)])
    poly[-1] -= float(P)
    return np.roots(poly)


def has_complex_pair(eigenvalues, tol: float = 1e-9) -> bool:
    return bool(np.any(np.abs(np.imag(eigenvalues)) > tol))


@dataclasses.dataclass
class SpectralVerdict:
    gamma: tuple[float, ...]
    slope_product: float
    h3_bound: float
    satisfies_h3: bool
    eigenvalues: np.ndarray
    complex_pair: bool

    def __str__(self) -> str:  # pragma: no cover - display only
        ev = ", ".join(f"{v:.4g}" for v in np.sort_complex(self.eigenvalues))
        return (
            f"gamma={self.gamma}, P={self.slope_product:g}: "
            f"H3 bound {self.h3_bound:.6g}, satisfied={self.satisfies_h3}, "
            f"complex pair={self.complex_pair}; eigenvalues {ev}"
        )


def h3_denominator(gamma) -> float:
    """``-g1 g2 g3 + (g1+g2+g3)(g1 g2 + g1 g3 + g2 g3)`` from ``eq:F3-bounds``."""
    g1, g2, g3 = gamma
    return -g1 * g2 * g3 + (g1 + g2 + g3) * (g1 * g2 + g1 * g3 + g2 * g3)


def check_h3_forces_complex_pair(gamma, slope_product: float) -> SpectralVerdict:
    """Does the ``H_3`` bound imply a complex stable pair?

    ``eq:F3-bounds`` is ``8 prod h < prod |Delta E_n| / D(gamma)``.  Writing the
    slope of the ``n``-th ramp across its window as ``Delta E_n / (2 h_n)``, the
    signed slope product is ``P = prod Delta E_n / (8 prod h)``, so the bound is
    exactly ``|P| > D(gamma)``.  The R3 chapter uses that to claim the cyclic
    equilibrium has a complex conjugate stable pair.
    """
    bound = h3_denominator(gamma)
    eigenvalues = eigenvalues_from_product(gamma, slope_product)
    return SpectralVerdict(
        gamma=tuple(float(g) for g in gamma),
        slope_product=float(slope_product),
        h3_bound=bound,
        satisfies_h3=abs(slope_product) > bound,
        eigenvalues=eigenvalues,
        complex_pair=has_complex_pair(eigenvalues),
    )


def search_real_spectrum_counterexamples(
    gamma_grid=None, product_grid=None, limit: int = 20
) -> list[SpectralVerdict]:
    """Parameters satisfying ``H_3`` whose cyclic equilibrium has real spectrum.

    Each hit refutes the claim that the ``H_3`` bound forces a complex pair.
    """
    if gamma_grid is None:
        values = [0.5, 1.0, 2.0, 5.0, 20.0, 100.0]
        gamma_grid = [g for g in itertools.product(values, repeat=3)]
    out: list[SpectralVerdict] = []
    for gamma in gamma_grid:
        bound = h3_denominator(gamma)
        products = product_grid or [
            bound * f for f in (1.001, 1.01, 1.1, 1.5, 2.0, 5.0)
        ]
        for P in products:
            verdict = check_h3_forces_complex_pair(gamma, P)
            if verdict.satisfies_h3 and not verdict.complex_pair:
                out.append(verdict)
                if len(out) >= limit:
                    return out
    return out


# ---------------------------------------------------------------------------
# swept surfaces: tangent vs normal
# ---------------------------------------------------------------------------


def swept_surface(generator, flow_field, *, samples: int = 25, steps: int = 60,
                  t_max: float = 1.0):
    """Sweep a curve under a flow, returning the surface as a grid of points.

    ``generator(s)`` parameterises the generating curve for ``s`` in ``[0,1]``;
    ``flow_field(x)`` is the field that sweeps it.
    """
    from scipy.integrate import solve_ivp

    rows = []
    for s in np.linspace(0.0, 1.0, samples):
        x0 = np.asarray(generator(float(s)), dtype=float)
        sol = solve_ivp(
            lambda _t, x: flow_field(x), (0.0, t_max), x0,
            dense_output=True, rtol=1e-9, atol=1e-12,
        )
        rows.append(sol.sol(np.linspace(0.0, t_max, steps)).T)
    return np.array(rows)  # (samples, steps, N)


def tangent_frame(surface: np.ndarray, i: int, j: int) -> tuple[np.ndarray, np.ndarray]:
    """Two independent tangents of the swept surface at grid point ``(i, j)``."""
    samples, steps, _ = surface.shape
    i0, i1 = max(i - 1, 0), min(i + 1, samples - 1)
    j0, j1 = max(j - 1, 0), min(j + 1, steps - 1)
    along_generator = surface[i1, j] - surface[i0, j]
    along_flow = surface[i, j1] - surface[i, j0]
    return along_generator, along_flow


def is_transverse(surface: np.ndarray, candidate, *, tol: float = 1e-6) -> dict:
    """Is ``candidate`` really transverse to the swept surface, or tangent to it?

    Reports the largest component of ``candidate`` inside the tangent plane; a
    genuinely transverse vector has a substantial normal component everywhere.
    """
    candidate = np.asarray(candidate, dtype=float)
    candidate = candidate / np.linalg.norm(candidate)
    worst_normal = np.inf
    samples, steps, _ = surface.shape
    for i in range(samples):
        for j in range(steps):
            t1, t2 = tangent_frame(surface, i, j)
            basis = np.array([t1, t2])
            q, _ = np.linalg.qr(basis.T)
            tangential = q @ (q.T @ candidate)
            normal_part = float(np.linalg.norm(candidate - tangential))
            worst_normal = min(worst_normal, normal_part)
    return {
        "min_normal_component": worst_normal,
        "transverse": worst_normal > tol,
    }


# ---------------------------------------------------------------------------
# locating the cycle cells
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class CycleCell:
    cell: int
    coords: tuple[int, ...]
    cycle: tuple[int, ...]
    signature: int
    regulation_map: dict

    @property
    def feedback(self) -> str:
        """``defn:feedback``: the sign of the product of rook components."""
        return "positive" if self.signature > 0 else "negative"


def cycle_cells_from_stg(stg: SpecCubicalBlowupGraph) -> list[CycleCell]:
    """Semi-opaque cells whose regulation map has a nontrivial cycle.

    The signature ``delta_sigma = prod_{i in S_sigma} Phi_i(xi, [v, 1])``
    (``defn:feedback``, PGradingR3v4.tex:28-46) distinguishes positive from
    negative feedback, which is what selects the R3 surface type.

    These are the only cells the P3 grading and the R3 geometry act on: for a
    wall labelling without them, ``F_3 = F_2`` and the R3 chapter has nothing
    to construct.
    """
    cc = stg.cubical_complex
    out: list[CycleCell] = []
    for cell in cc:
        if cc.rightfringe(cell):
            continue
        cycles = stg.nontrivial_cycles(cell)
        if not cycles:
            continue
        coords = tuple(cc.coordinates(cell))
        top = cc.cell_index(list(coords), stg.top_shape)
        for cycle in cycles:
            signature = 1
            for n in cycle:
                signature *= stg.rook_field_component(cell, top, n)
            out.append(
                CycleCell(
                    cell=cell,
                    coords=coords,
                    cycle=tuple(cycle),
                    signature=signature,
                    regulation_map=stg.active_regulation_map(cell),
                )
            )
    return out


def cycle_cells(
    system: RampSystem, *, spec: Spec = PAPER, outer: str = "global_bound"
) -> list[CycleCell]:
    """:func:`cycle_cells_from_stg` for a directly specified ramp system."""
    labelling, num_thresholds = system.wall_labelling(outer=outer)
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=spec, level=3
    )
    return cycle_cells_from_stg(stg)


def cycle_cells_for_parameter(name: str, index: int, *, spec: Spec = PAPER):
    """:func:`cycle_cells_from_stg` for a DSGRN parameter node."""
    from .networks import parameter

    stg = SpecCubicalBlowupGraph(parameter(name, index), spec=spec, level=3)
    return cycle_cells_from_stg(stg), stg


def report(system: RampSystem, *, spec: Spec = PAPER) -> dict:
    """Summarise the R3 cycle data for one ramp system."""
    cells = cycle_cells(system, spec=spec)
    by_length: dict[int, int] = {}
    by_feedback: dict[str, int] = {}
    for c in cells:
        by_length[len(c.cycle)] = by_length.get(len(c.cycle), 0) + 1
        by_feedback[c.feedback] = by_feedback.get(c.feedback, 0) + 1
    return {
        "cycle_cells": len(cells),
        "by_cycle_length": by_length,
        "by_feedback": by_feedback,
    }
