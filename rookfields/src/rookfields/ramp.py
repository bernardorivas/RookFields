"""Ramp systems: vector field, flow, and the cell-complex bookkeeping.

``defn:ramp_system`` (RampSystemsv4.tex:66-112) is

    dot x_n = F_n(x) = -gamma_n x_n + E_n(x),   gamma_n > 0,

with each ``E_n`` an interaction function of one-variable ramps.

Indexing convention.  For coordinate ``n`` all output thresholds are pooled over
the targets of ``n`` and sorted (``rem:K(n)``, RampSystemsv4.tex:254), giving
``K(n)`` thresholds.  The manuscript indexes them ``theta_{m_k, n, j_k}`` for
``k = 1..K(n)``, with sentinels ``theta_{m_0} = 0`` and
``theta_{m_{K+1}} = GB_n`` (``defn:ramp-wall-labeling``, :375-390).  A cell
coordinate ``v_n = k`` therefore means:

* ``w_n = 0`` (inessential): the cell sits on threshold ``theta_{m_k}``;
* ``w_n = 1`` (essential): the cell spans the plateau box between
  ``theta_{m_k}`` and ``theta_{m_{k+1}}``.

This matches ``DSGRN_utils.WallLabelling.ramp_system_wall_labelling``, whose
box ``k`` has left wall ``theta_s[k-1]`` and right wall ``theta_s[k]``.
"""

from __future__ import annotations

import dataclasses
import itertools
from functools import cached_property

import numpy as np
from DSGRN_utils.RampFunction import ramp_function

from .networks import RampSpec


@dataclasses.dataclass
class RampSystem:
    """A ramp system with the bookkeeping the analytic hypotheses need."""

    spec: RampSpec

    # -- basic structure -------------------------------------------------

    @property
    def dim(self) -> int:
        return self.spec.dim

    @property
    def gamma(self) -> tuple[float, ...]:
        return self.spec.gamma

    @cached_property
    def _tables(self):
        blank = lambda t: [[list(e) if e else [] for e in row] for row in t]
        return blank(self.spec.nu), blank(self.spec.theta), blank(self.spec.h)

    @cached_property
    def sorted_thresholds(self) -> list[list[float]]:
        """``theta_{m_1,n} < ... < theta_{m_{K(n)},n}``, the output thresholds of ``n``."""
        _nu, theta, h = self._tables
        out = []
        for n in range(self.dim):
            pooled = [(t, hh) for k in range(self.dim)
                      for t, hh in zip(theta[k][n], h[k][n])]
            pooled.sort()
            out.append([t for t, _ in pooled])
        return out

    @cached_property
    def sorted_widths(self) -> list[list[float]]:
        """The half-widths ``h_{m_k,n}``, in the same order."""
        _nu, theta, h = self._tables
        out = []
        for n in range(self.dim):
            pooled = [(t, hh) for k in range(self.dim)
                      for t, hh in zip(theta[k][n], h[k][n])]
            pooled.sort()
            out.append([hh for _, hh in pooled])
        return out

    @property
    def num_thresholds(self) -> list[int]:
        """``K(n)``."""
        return [len(t) for t in self.sorted_thresholds]

    # -- the vector field -------------------------------------------------

    def production(self, x) -> np.ndarray:
        """``E(x)``, the interaction term (type-I, additive)."""
        return np.asarray(self.spec.production(list(np.asarray(x, dtype=float))), dtype=float)

    def vector_field(self, x) -> np.ndarray:
        """``F(x) = -Gamma x + E(x)``."""
        x = np.asarray(x, dtype=float)
        return -np.asarray(self.gamma, dtype=float) * x + self.production(x)

    def jacobian(self, x, eps: float = 1e-7) -> np.ndarray:
        """Numerical Jacobian.

        The field is only piecewise smooth -- it is not differentiable at ramp
        endpoints (a defect the manuscript records for the R2 derivative
        claims) -- so a central difference is the honest tool here.
        """
        x = np.asarray(x, dtype=float)
        out = np.zeros((self.dim, self.dim))
        for j in range(self.dim):
            e = np.zeros(self.dim)
            e[j] = eps
            out[:, j] = (self.vector_field(x + e) - self.vector_field(x - e)) / (2 * eps)
        return out

    # -- global bound and box (eq:GAB, lem:global-bound) ------------------

    @cached_property
    def max_production(self) -> np.ndarray:
        """``M^E_n = max E_n``, respecting the interaction structure."""
        return np.asarray(self.spec.max_production, dtype=float)

    @cached_property
    def global_bound(self) -> np.ndarray:
        """``GB_n = max{ M^E_n / gamma_n, max ramp endpoint } + 1`` (eq:GAB)."""
        return np.asarray(self.spec.global_bound, dtype=float)

    def absorption_bound(self, x0, t) -> np.ndarray:
        """``x_n(t) <= G^b_n + (x_n(0) - G^b_n)^+ e^{-gamma_n t}``.

        The repaired scalar-comparison estimate.  The manuscript's product-box
        argument is recorded as defective: being outside the box does not make
        every component of the field negative, and ``E`` is not constant there.
        """
        x0 = np.asarray(x0, dtype=float)
        gb = self.max_production / np.asarray(self.gamma, dtype=float)
        excess = np.maximum(x0 - gb, 0.0)
        return gb + excess * np.exp(-np.asarray(self.gamma, dtype=float) * t)

    # -- thresholds with sentinels ---------------------------------------

    def theta(self, n: int, k: int) -> float:
        """``theta_{m_k, n, j_k}``, 1-indexed, with the wall-labeling sentinels."""
        K = self.num_thresholds[n]
        if k <= 0:
            return 0.0
        if k > K:
            return float(self.global_bound[n])
        return self.sorted_thresholds[n][k - 1]

    def width(self, n: int, k: int) -> float:
        """``h_{m_k, n, j_k}``, 1-indexed; the sentinels carry zero width."""
        K = self.num_thresholds[n]
        if k <= 0 or k > K:
            return 0.0
        return self.sorted_widths[n][k - 1]

    # -- I_n, Xi_n, L_n, U_n ----------------------------------------------

    def interval(self, coords, shape_bits, n: int) -> tuple[float, float]:
        """``I_n(xi)`` (RampSystemsv4.tex:534-540).

        ``shape_bits[n] == 0`` means ``n in J_i(xi)``.
        """
        k = coords[n]
        if shape_bits[n] == 0:
            return (self.theta(n, k) - self.width(n, k), self.theta(n, k) + self.width(n, k))
        return (self.theta(n, k) + self.width(n, k), self.theta(n, k + 1) - self.width(n, k + 1))

    def interval_length(self, coords, shape_bits, n: int) -> float:
        """``Xi_n(xi)`` (eq:IntervalLength)."""
        lo, hi = self.interval(coords, shape_bits, n)
        return hi - lo

    def plateau_point(self, box_coords) -> np.ndarray:
        """An interior point of the plateau box ``D_v`` where every ramp is flat."""
        out = np.zeros(self.dim)
        for n, k in enumerate(box_coords):
            lo = self.theta(n, k) + self.width(n, k)
            hi = self.theta(n, k + 1) - self.width(n, k + 1)
            out[n] = 0.5 * (lo + hi)
        return out

    def production_on_box(self, box_coords) -> np.ndarray:
        """``E(D_v)``: the constant production value on a plateau box."""
        return self.production(self.plateau_point(box_coords))

    def top_star_boxes(self, coords, shape_bits) -> list[tuple[int, ...]]:
        """``Top_X(xi)`` as plateau-box coordinate tuples.

        A cell extends to the boxes obtained by dropping each inessential
        coordinate to either side of its threshold.
        """
        options = []
        for n in range(self.dim):
            if shape_bits[n] == 1:
                options.append([coords[n]])
            else:
                options.append([coords[n] - 1, coords[n]])
        return [
            tuple(c)
            for c in itertools.product(*options)
            if all(0 <= c[n] <= self.num_thresholds[n] for n in range(self.dim))
        ]

    def _component_extremes(self, coords, shape_bits, n: int) -> tuple[float, float]:
        """``(min, max)`` of ``| -gamma_n x_n + E_n(mu) |`` over ``I_n(xi)`` and ``Top_X(xi)``."""
        lo, hi = self.interval(coords, shape_bits, n)
        lows, highs = [], []
        for box in self.top_star_boxes(coords, shape_bits):
            e_n = float(self.production_on_box(box)[n])
            # -gamma_n x_n + e_n is affine and decreasing in x_n
            a = -self.gamma[n] * lo + e_n
            b = -self.gamma[n] * hi + e_n
            values = [abs(a), abs(b)]
            if a * b <= 0:  # the affine function vanishes inside the interval
                lows.append(0.0)
            else:
                lows.append(min(values))
            highs.append(max(values))
        return (min(lows), max(highs))

    def L(self, coords, shape_bits, n: int) -> float:
        """``L_n(xi)`` (eq:Ln)."""
        return self._component_extremes(coords, shape_bits, n)[0]

    def U(self, coords, shape_bits, n: int) -> float:
        """``U_n(xi)`` (eq:Un)."""
        return self._component_extremes(coords, shape_bits, n)[1]

    # -- flow --------------------------------------------------------------

    def flow(self, x0, t_span=(0.0, 50.0), *, max_step=None, rtol=1e-8, atol=1e-10, **kw):
        """Integrate the ramp field.

        The field is Lipschitz but only piecewise smooth, so a tight tolerance
        and a bounded step keep the solver from stepping across a ramp window.
        """
        from scipy.integrate import solve_ivp

        if max_step is None:
            smallest = min(
                (2 * w for row in self._tables[2] for e in row for w in e), default=0.1
            )
            max_step = max(smallest / 4.0, 1e-3)
        return solve_ivp(
            lambda _t, x: self.vector_field(x),
            t_span,
            np.asarray(x0, dtype=float),
            max_step=max_step,
            rtol=rtol,
            atol=atol,
            dense_output=True,
            **kw,
        )

    # -- the ramp induced wall labelling -----------------------------------

    def wall_labelling(self, *, outer: str = "global_bound"):
        """``defn:ramp-wall-labeling`` (RampSystemsv4.tex:375-390).

        ``omega(xi, mu) = sgn(-gamma_n theta_{m_{k_n},n,j_{k_n}} + E_n(mu))``
        with sentinels ``theta_{m_0} = 0`` and ``theta_{m_{K(n)+1}} = GB_n``.

        ``outer`` selects the outermost sentinel:

        ``"global_bound"``
            ``GB_n``, as the definition states.
        ``"legacy"``
            ``theta_{m_K} + 10 h_{m_K}``, which is what
            ``DSGRN_utils.WallLabelling.ramp_system_wall_labelling`` (:44) uses.
            That surrogate depends on ``h``, so it contradicts
            ``prop:wall-labeling-const`` (the labelling is constant in ``h``
            on ``H_0``) and can produce a labelling that is not strongly
            dissipative once ``h`` is small enough to satisfy ``H_1``.

        Returns ``(labelling, num_thresholds)`` in the DSGRN_utils bitmask
        encoding: bit ``n`` set when the left ``n``-wall is absorbing, bit
        ``n + dim`` set when the right ``n``-wall is.
        """
        if outer == "legacy":
            return self.spec.wall_labelling(legacy=True)
        if outer != "global_bound":
            raise ValueError(f"unknown outer sentinel {outer!r}")

        K = self.num_thresholds
        num_boxes = [k + 1 for k in K]
        pv = [1]
        for b in num_boxes:
            pv.append(pv[-1] * b)

        labelling = [0] * pv[-1]
        for box in itertools.product(*[range(b) for b in num_boxes]):
            index = sum(c * pv[n] for n, c in enumerate(box))
            e = self.production_on_box(box)
            label = 0
            for n in range(self.dim):
                left = -self.gamma[n] * self.theta(n, box[n]) + float(e[n])
                right = -self.gamma[n] * self.theta(n, box[n] + 1) + float(e[n])
                if left == 0.0 or right == 0.0:
                    raise ValueError(
                        f"ramp system vanishes on a wall of box {box}, coordinate {n}: "
                        "the parameter is outside Lambda(R)"
                    )
                if left < 0:  # omega = -1: flow moves left, out of the box
                    label |= 1 << n
                if right > 0:  # omega = +1: flow moves right, out of the box
                    label |= 1 << (n + self.dim)
            labelling[index] = label
        return labelling, list(K)


def scalar_ramp(x: float, nu, theta, h) -> float:
    """One ramp factor ``r(x; nu, theta, h)`` (``defn:rampfunction``)."""
    return ramp_function(x, list(nu), list(theta), list(h))
