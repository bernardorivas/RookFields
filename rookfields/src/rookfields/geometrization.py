"""Rectangular geometrization of the blowup complex, and the Janus subdivision.

``RectangularGeo.tex`` builds a geometrization of ``X_b`` in four steps:

1. sentinel thresholds at both ends (eq:m0nj0, eq:mKnjK, :21-29) chosen so the
   embedded complex exactly fills ``prod [0, GB_n]``;
2. the vertex embedding ``e`` (``defn:vertexEmb``, :56-67);
3. rectangles ``Rec([v,w]) = prod I_n`` and model faces ``K(w)`` (:79-106);
4. inductive affine charts ``a_zeta : K(w) -> Rec(zeta)`` agreeing on faces
   (:108-122), composed with any ``Geo`` preserving each rectangle to give
   ``recG`` (``defn:rectGeoXb``, :135-149).

``Janus3v5.tex`` then subdivides each blowup interval into eight
(``defn:JanusComplex``, :9-20), with the corrected closed cubical carrier
``eq:barS`` and the minimal-coarse-carrier chart ``a_zeta(y) = (k + y)/8``
(:109-121).

Nothing in the project implemented any of this before.
"""

from __future__ import annotations

import dataclasses
import itertools
from functools import lru_cache

import numpy as np

from .ramp import RampSystem


@dataclasses.dataclass
class RectangularGeometrization:
    """``recG(gamma, nu, theta, h)`` for a ramp system.

    The identity is used for ``Geo``; ``defn:rectGeoXb`` allows any homeomorphism
    preserving each rectangle, and every statement checked here is invariant
    under that choice because it only uses the rectangles themselves.
    """

    system: RampSystem

    # -- sentinels (eq:m0nj0, eq:mKnjK) ----------------------------------

    def theta(self, n: int, k: int) -> float:
        """``theta_{m_k, n, j_k}`` with the *geometrization* sentinels.

        These differ from the wall-labelling sentinels of
        ``defn:ramp-wall-labeling`` (which are ``0`` and ``GB_n``): here they
        are chosen so that ``theta_0 - h_0 = 0`` and
        ``theta_{K+1} + h_{K+1} = GB_n``.
        """
        s = self.system
        K = s.num_thresholds[n]
        if k <= 0:
            return (s.sorted_thresholds[n][0] - s.sorted_widths[n][0]) / 4.0
        if k > K:
            return float(s.global_bound[n]) - 0.25
        return s.sorted_thresholds[n][k - 1]

    def width(self, n: int, k: int) -> float:
        s = self.system
        K = s.num_thresholds[n]
        if k <= 0:
            return (s.sorted_thresholds[n][0] - s.sorted_widths[n][0]) / 4.0
        if k > K:
            return 0.25
        return s.sorted_widths[n][k - 1]

    # -- defn:vertexEmb ---------------------------------------------------

    def embed_coordinate(self, n: int, blowup_coord: int) -> float:
        """``e_n(v)`` for a blowup vertex coordinate ``v_n = 2 k_n + j_n``."""
        k, j = divmod(int(blowup_coord), 2)
        return self.theta(n, k) + (self.width(n, k) if j == 1 else -self.width(n, k))

    def embed_vertex(self, blowup_coords) -> np.ndarray:
        """``e(v)`` for a vertex of ``X_b``."""
        return np.array(
            [self.embed_coordinate(n, c) for n, c in enumerate(blowup_coords)],
            dtype=float,
        )

    # -- rectangles (eq:Rxi) ----------------------------------------------

    def rectangle(self, blowup_coords, shape_bits) -> list[tuple[float, float]]:
        """``Rec([v, w])`` as a list of intervals; degenerate where ``w_n = 0``."""
        out = []
        for n, (c, w) in enumerate(zip(blowup_coords, shape_bits)):
            lo = self.embed_coordinate(n, c)
            hi = self.embed_coordinate(n, c + 1) if w else lo
            out.append((lo, hi))
        return out

    def chart(self, blowup_coords, shape_bits):
        """``a_zeta : K(w) -> Rec(zeta)``, the inductive affine map.

        The induction of ``RectangularGeo.tex:108-122`` has a closed form: the
        affine map is coordinatewise, sending the unit interval in each
        essential direction onto that direction's embedded interval and every
        inessential direction to its single vertex value.  That formula
        restricts correctly to every face, which is exactly what the induction
        demands.
        """
        rect = self.rectangle(blowup_coords, shape_bits)

        def a(y) -> np.ndarray:
            y = np.atleast_1d(np.asarray(y, dtype=float))
            out = np.empty(len(rect))
            cursor = 0
            for n, (lo, hi) in enumerate(rect):
                if shape_bits[n]:
                    out[n] = lo + (hi - lo) * y[cursor]
                    cursor += 1
                else:
                    out[n] = lo
            return out

        return a

    def cell_image(self, blowup_coords, shape_bits) -> list[tuple[float, float]]:
        """``g(zeta)``, the embedded image of a blowup cell."""
        return self.rectangle(blowup_coords, shape_bits)

    def blowup_cell_image_of(self, cc_coords, cc_shape_bits) -> list[tuple[float, float]]:
        """``g(b(xi))`` for a cell ``xi`` of ``X`` (``rem:g(b(xi))``).

        Agrees with ``I_n(xi)`` of ``RampSystemsv4.tex:534`` on every cell that
        does not reference a sentinel index -- see :meth:`uses_sentinel`.  The
        two chapters give the outer sentinel different values, so the two
        formulas disagree at the outermost cells; the test suite pins that down
        rather than papering over it.
        """
        blowup_coords = [2 * v + w for v, w in zip(cc_coords, cc_shape_bits)]
        return self.rectangle(blowup_coords, [1] * len(blowup_coords))

    def uses_sentinel(self, cc_coords, cc_shape_bits) -> bool:
        """Whether ``I_n`` for this cell references a sentinel threshold index.

        ``defn:ramp-wall-labeling`` (RampSystemsv4.tex:388) sets
        ``theta_{m_{K+1}} = GB_n``; ``eq:mKnjK`` (RectangularGeo.tex:27) sets
        ``theta_{m_{K+1}} = GB_n - 1/4`` with ``h_{m_{K+1}} = 1/4`` so that
        ``theta + h = GB_n``.  Likewise at the lower end.  Both are internally
        consistent, but a formula quantified over ``v in prod {1..K(n)}`` and
        involving ``theta_{m_{k+1}}`` -- such as ``I_n``, ``Xi_n`` and hence the
        ``H_2`` inequality -- picks up whichever convention its chapter uses.
        """
        s = self.system
        for n, (v, w) in enumerate(zip(cc_coords, cc_shape_bits)):
            top = v + (1 if w else 0)
            if v <= 0 or top > s.num_thresholds[n]:
                return True
        return False

    # -- sanity ------------------------------------------------------------

    def covers_bounding_box(self, tol: float = 1e-12) -> bool:
        """``eq:rectangular_geo_boundary``: the complex fills ``prod [0, GB_n]``."""
        s = self.system
        for n in range(s.dim):
            if abs(self.embed_coordinate(n, 0) - 0.0) > tol:
                return False
            top = 2 * (s.num_thresholds[n] + 1) + 1
            if abs(self.embed_coordinate(n, top) - float(s.global_bound[n])) > tol:
                return False
        return True

    def intervals_are_ordered(self) -> bool:
        """Consecutive embedded vertices strictly increase, so no cell degenerates."""
        s = self.system
        for n in range(s.dim):
            top = 2 * (s.num_thresholds[n] + 1) + 1
            values = [self.embed_coordinate(n, c) for c in range(top + 1)]
            if any(b <= a for a, b in zip(values, values[1:])):
                return False
        return True


# ---------------------------------------------------------------------------
# Janus complex
# ---------------------------------------------------------------------------

JANUS_FACTOR = 8


@dataclasses.dataclass
class JanusComplex:
    """The factor-eight subdivision ``X_J`` of ``X_b`` (``defn:JanusComplex``).

    Coordinates run over ``0 .. 8 (Kbar(n) + 1)`` where
    ``Kbar(n) = 2 (K(n) + 1)`` is the blowup range; the endpoint
    ``8 (Kbar(n) + 1)`` is the corrected one recorded in the R2 dossier
    (the earlier draft stopped one step short, which double-assigned the fine
    top cell based at ``8(v+1)`` to both adjacent coarse cells).
    """

    system: RampSystem

    @property
    def dim(self) -> int:
        return self.system.dim

    @property
    def blowup_range(self) -> list[int]:
        """``Kbar(n)``: the largest blowup coordinate."""
        return [2 * (k + 1) + 1 for k in self.system.num_thresholds]

    @property
    def vertex_range(self) -> list[int]:
        """The largest Janus vertex coordinate, ``8 (Kbar(n) + 1)``."""
        return [JANUS_FACTOR * (r + 1) for r in self.blowup_range]

    def carrier(self, blowup_coords, shape_bits) -> list[tuple[range, ...]]:
        """``bar S([v, w])`` (eq:barS), the closed cubical carrier.

        ``{ [u, q] : q <= w, 8 v_n <= u_n <= 8 (v_n + w_n) - q_n }``.
        Returned as, for each ``q <= w``, the tuple of allowed ``u_n`` ranges.
        """
        dim = self.dim
        out = []
        for q in itertools.product(*[[0, 1] if w else [0] for w in shape_bits]):
            ranges = []
            for n in range(dim):
                lo = JANUS_FACTOR * blowup_coords[n]
                hi = JANUS_FACTOR * (blowup_coords[n] + shape_bits[n]) - q[n]
                ranges.append(range(lo, hi + 1))
            out.append((q, tuple(ranges)))
        return out

    def top_slice(self, blowup_coords, shape_bits) -> set[tuple[int, ...]]:
        """``bar S^{(N)}``: the top-dimensional fine cells of the carrier.

        These are the ``q = w`` members, i.e. the fine cells of full dimension
        inside the coarse cell.
        """
        for q, ranges in self.carrier(blowup_coords, shape_bits):
            if list(q) == list(shape_bits):
                return {tuple(u) for u in itertools.product(*ranges)}
        return set()

    def chart_offset(self, fine_coords, coarse_coords) -> np.ndarray:
        """``bar a_zeta(y) = (k + y) / 8`` relative to the minimal coarse carrier."""
        k = np.asarray(fine_coords, dtype=float) - JANUS_FACTOR * np.asarray(
            coarse_coords, dtype=float
        )
        return k / JANUS_FACTOR

    def _blowup_top_cells(self) -> tuple:
        return tuple(
            itertools.product(*[range(r + 1) for r in self.blowup_range])
        )

    def adjacent_top_slices_are_disjoint(self, limit: int | None = 400) -> bool:
        """Regression test: adjacent coarse top cells share no fine top cell.

        The pre-correction ``eq:barS`` assigned the fine top cell based at
        ``8(v+1)`` to both neighbours, which is what broke well-definedness of
        ``pi_J`` and of the top-cell branch of ``pi_2``.
        """
        cells = self._blowup_top_cells()
        if limit is not None:
            cells = cells[:limit]
        ones = [1] * self.dim
        for coords in cells:
            mine = self.top_slice(list(coords), ones)
            for n in range(self.dim):
                neighbour = list(coords)
                neighbour[n] += 1
                if neighbour[n] > self.blowup_range[n]:
                    continue
                if mine & self.top_slice(neighbour, ones):
                    return False
        return True

    def subdivision_chain_endpoint_ok(self) -> bool:
        """The one-dimensional subdivision chain reaches the right endpoint.

        A coarse edge from ``v_n`` to ``v_n + 1`` must contain exactly the eight
        fine edges based at ``8 v_n, ..., 8 v_n + 7``.
        """
        ones = [1] * self.dim
        for n in range(self.dim):
            coords = [0] * self.dim
            for v in range(self.blowup_range[n] + 1):
                coords[n] = v
                edge_shape = [0] * self.dim
                edge_shape[n] = 1
                bases = sorted(
                    u[n] for u in self.top_slice(coords, edge_shape)
                )
                if bases != list(range(JANUS_FACTOR * v, JANUS_FACTOR * v + JANUS_FACTOR)):
                    return False
            coords[n] = 0
        return True
