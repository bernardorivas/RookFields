"""Abstract wall labelings: validity, dissipativity, and enumeration.

``def:wall_labeling`` (RookFields6.tex:514-536) constrains a function
``omega : W(X) -> {+-1}`` by requiring, at every vertex ``sigma``, a local
inducement map ``o~_sigma : {1..N} -> {1..N}`` with

  (i)  for ``n``-adjacent ``mu, mu' in Top_X(sigma)`` and ``k != n``,
       ``k != o~_sigma(n)``:
       ``omega(mu^-_k, mu) = omega(mu'^-_k, mu')`` and
       ``omega(mu^+_k, mu) = omega(mu'^+_k, mu')``;
  (ii) for ``n``-walls ``(xi, mu), (xi, mu') in W(sigma)`` with
       ``n != o~_sigma(n)``: ``omega(xi, mu) = omega(xi, mu')``.

``defn:dissipativewall`` (RookFields6.tex:618-628) additionally forces
``omega(xi, mu) = -p(xi, mu)`` on every wall with ``xi in bdy(X)``.

Theorems such as ``thm:2dDone`` (CombinatorialDynamics.tex:1135) are stated for
*every* wall labeling on an abstract cubical complex, not only for the ones
DSGRN produces.  Testing them therefore needs a generator of abstract labelings,
which is what this module provides.

The encoding matches ``DSGRN_utils``: one integer per top cell of ``X``, with
bit ``n`` set when the left ``n``-wall is absorbing and bit ``n + dim`` set when
the right ``n``-wall is.
"""

from __future__ import annotations

import itertools
import random
from functools import lru_cache

import pychomp


class WallLabeling:
    """A candidate wall labeling on ``X(I)`` with ``I = prod {0..K(n)+1}``."""

    def __init__(self, labelling: list[int], num_thresholds: list[int]):
        self.labelling = list(labelling)
        self.num_thresholds = list(num_thresholds)
        self.dim = len(num_thresholds)
        self.num_boxes = [k + 1 for k in num_thresholds]
        self.limits = [k + 1 for k in num_thresholds]
        # One padding layer, exactly as CubicalBlowupGraph.__init__ (:66) does,
        # so the honest complex is a subcomplex.
        self.complex = pychomp.CubicalComplex([k + 1 for k in self.num_boxes])
        self.pv = [1]
        for k in self.num_boxes:
            self.pv.append(self.pv[-1] * k)
        self.top_shape = 2**self.dim - 1

    # -- basic accessors ------------------------------------------------

    def is_fringe(self, cell) -> bool:
        return self.complex.rightfringe(cell)

    def label(self, top_cell, n: int, side: int) -> int:
        """``omega(mu^side_n, mu)``, matching ``CubicalBlowupGraph.wall_label``."""
        coords = self.complex.coordinates(top_cell)
        if any(coords[k] == self.limits[k] for k in range(self.dim)):
            return side  # fringe walls are absorbing
        index = sum(c * self.pv[k] for k, c in enumerate(coords))
        mask = 1 << (n + (self.dim if side == 1 else 0))
        return side if self.labelling[index] & mask else -side

    @lru_cache(maxsize=None)
    def _real_top_cells(self) -> tuple:
        return tuple(c for c in self.complex(self.dim) if not self.is_fringe(c))

    def vertices(self) -> list:
        return [c for c in self.complex(0) if not self.is_fringe(c)]

    def top_star(self, cell) -> list:
        return [c for c in self.complex.topstar(cell) if not self.is_fringe(c)]

    # -- def:wall_labeling ----------------------------------------------

    def _adjacent_top_pairs(self, cell, n):
        """``n``-adjacent pairs of top cells in ``Top_X(cell)``."""
        star = set(self.top_star(cell))
        jump = [1]
        for k in self.complex.boxes():
            jump.append(jump[-1] * k)
        return [(m, m + jump[n]) for m in star if m + jump[n] in star]

    def _n_walls_at(self, sigma, n):
        """``n``-walls ``xi`` with ``sigma <= xi``, with their two top cells."""
        walls = []
        for cell in self.complex.star({sigma}):
            if self.is_fringe(cell):
                continue
            if self.complex.cell_dim(cell) != self.dim - 1:
                continue
            shape = self.complex.cell_shape(cell)
            # an n-wall has n inessential and every other direction essential
            if shape & (1 << n):
                continue
            tops = self.top_star(cell)
            if len(tops) == 2:
                walls.append((cell, tops))
        return walls

    def _omega(self, wall, top) -> int:
        """``omega(xi, mu)`` for a codimension-one ``xi`` and ``mu`` above it."""
        shape = self.complex.cell_shape(wall)
        n = next(k for k in range(self.dim) if not (shape & (1 << k)))
        wc = self.complex.coordinates(wall)
        tc = self.complex.coordinates(top)
        side = -1 if wc[n] == tc[n] else 1
        return self.label(top, n, side)

    def local_inducement_maps(self, sigma) -> list[dict[int, int]]:
        """Every ``o~_sigma`` compatible with this labeling at ``sigma``.

        Empty iff ``def:wall_labeling`` fails at ``sigma``.
        """
        # Condition (i) and (ii) constrain o~_sigma(n) independently for each n,
        # so the admissible maps form a product and we test one n at a time.
        allowed: list[list[int]] = []
        for n in range(self.dim):
            targets = []
            for candidate in range(self.dim):
                if self._conditions_hold(sigma, n, candidate):
                    targets.append(candidate)
            if not targets:
                return []
            allowed.append(targets)
        return [dict(enumerate(choice)) for choice in itertools.product(*allowed)]

    def _conditions_hold(self, sigma, n: int, candidate: int) -> bool:
        # (i)
        for mu, mu_prime in self._adjacent_top_pairs(sigma, n):
            for k in range(self.dim):
                if k == n or k == candidate:
                    continue
                if self.label(mu, k, -1) != self.label(mu_prime, k, -1):
                    return False
                if self.label(mu, k, 1) != self.label(mu_prime, k, 1):
                    return False
        # (ii)
        if n != candidate:
            for wall, tops in self._n_walls_at(sigma, n):
                values = {self._omega(wall, t) for t in tops}
                if len(values) > 1:
                    return False
        return True

    def is_valid(self) -> bool:
        """``def:wall_labeling``: a local inducement map exists at every vertex."""
        return all(self.local_inducement_maps(v) for v in self.vertices())

    # -- defn:dissipativewall -------------------------------------------

    def is_strongly_dissipative(self) -> bool:
        """``omega(xi, mu) = -p(xi, mu)`` for every wall with ``xi in bdy(X)``.

        Concretely: the outermost walls must point inward, i.e. their absorbing
        bit must be clear.
        """
        for top in self._real_top_cells():
            coords = self.complex.coordinates(top)
            for n in range(self.dim):
                if coords[n] == 0 and self.label(top, n, -1) == -1:
                    return False
                if coords[n] == self.num_thresholds[n] and self.label(top, n, 1) == 1:
                    return False
        return True


def _dissipativity_mask(num_thresholds: list[int]) -> tuple[list[int], list[int]]:
    """Bits forced clear by strong dissipativity, per top-cell label slot."""
    dim = len(num_thresholds)
    num_boxes = [k + 1 for k in num_thresholds]
    pv = [1]
    for k in num_boxes:
        pv.append(pv[-1] * k)
    forced_clear = [0] * (pv[-1])
    for coords in itertools.product(*[range(b) for b in num_boxes]):
        index = sum(c * pv[k] for k, c in enumerate(coords))
        mask = 0
        for n in range(dim):
            if coords[n] == 0:
                mask |= 1 << n
            if coords[n] == num_thresholds[n]:
                mask |= 1 << (n + dim)
        forced_clear[index] = mask
    return forced_clear, num_boxes


def enumerate_wall_labelings(
    num_thresholds: list[int],
    *,
    dissipative: bool = True,
    valid_only: bool = True,
    limit: int | None = None,
    seed: int | None = None,
):
    """Yield ``WallLabeling`` objects on ``X(I)``.

    With ``seed`` the free bits are sampled uniformly at random (``limit``
    samples); otherwise the free bits are enumerated exhaustively in order.
    """
    dim = len(num_thresholds)
    forced_clear, num_boxes = _dissipativity_mask(num_thresholds)
    n_slots = len(forced_clear)
    width = 2 * dim

    free_bits: list[tuple[int, int]] = []
    for slot in range(n_slots):
        for bit in range(width):
            if dissipative and (forced_clear[slot] & (1 << bit)):
                continue
            free_bits.append((slot, bit))

    def build(assignment: int) -> WallLabeling:
        labelling = [0] * n_slots
        for position, (slot, bit) in enumerate(free_bits):
            if assignment & (1 << position):
                labelling[slot] |= 1 << bit
        return WallLabeling(labelling, num_thresholds)

    total = 1 << len(free_bits)
    if seed is not None:
        rng = random.Random(seed)
        produced = 0
        target = limit if limit is not None else 1000
        while produced < target:
            candidate = build(rng.randrange(total))
            if valid_only and not candidate.is_valid():
                continue
            if dissipative and not candidate.is_strongly_dissipative():
                continue
            produced += 1
            yield candidate
        return

    produced = 0
    for assignment in range(total):
        candidate = build(assignment)
        if valid_only and not candidate.is_valid():
            continue
        if dissipative and not candidate.is_strongly_dissipative():
            continue
        produced += 1
        yield candidate
        if limit is not None and produced >= limit:
            return
