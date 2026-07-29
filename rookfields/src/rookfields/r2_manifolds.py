"""GO-manifolds: the geometric objects the R2 chapter constructs.

``GlobalDynR2v9.tex`` replaces selected rectangular ``(N-1)``-faces by sheets
swept out from a nullcline under a perturbed field.  For a codimension-two
indecisive-drift pair ``(xi, xi')`` with GO-pair ``(n_g, n_o)``:

* ``defn:admissible-delta`` (:147-172) shrinks or enlarges every non-immutable
  interval by ``delta``, giving the control region

      Q^delta(xi, I) = prod_{n in I} I_n(b(xi)) x prod_{n not in I} I_n^delta(xi),
      Q^delta(xi, xi') = Q^delta(xi, I) cap Q^delta(xi', I);

* ``defn:Ixi`` (:173-182) fixes the immutable indices ``I = {n_g, n_o}`` in
  codimension two;

* the base is the ``delta``-nullcline
  ``Nul_{n_o}^delta(xi,xi') = { x in Q^delta : -gamma_{n_o} x_{n_o} + E_{n_o}(x) = 0 }``
  (:203), which ``cor:nullcline-graph-delta`` presents as a graph in the
  ``n_g`` direction;

* ``defn:GO-manifold`` (:723-780) flows that base under ``eq:rampSysPerturbed``

      dot x_n = 0                                              n != n_o, n_g
      dot x_{n_o} = -gamma_{n_o} x_{n_o} + E_{n_o}(x)
      dot x_{n_g} = (-gamma_{n_g} - r_{n_g} eps) x_{n_g} + E_{n_g}(x)

  until the directed ``n_g``-wall at
  ``theta_{n_o,n_g,j} + r_{n_g} h_{n_o,n_g,j}`` is hit, and takes the union of
  the trajectories.

Since ``n_o in J_i(xi) cap J_e(xi')``, the two ``n_o``-intervals meet in a
single point, so ``Q^delta(xi, xi')`` is degenerate in ``n_o`` and the base has
dimension ``N-2``.  Sweeping it adds one dimension, so the GO-manifold is the
codimension-one sheet that replaces the rectangular face.
"""

from __future__ import annotations

import dataclasses
import itertools

import numpy as np
from scipy.optimize import brentq

from .blowup import SpecCubicalBlowupGraph
from .geometrization import RectangularGeometrization
from .ramp import RampSystem
from .spec import PAPER, Spec


@dataclasses.dataclass
class GOPair:
    """A codimension-two indecisive-drift pair, with everything it determines."""

    face: int
    coface: int
    n_grad: int
    n_opaque: int
    face_coords: tuple[int, ...]
    face_shape: tuple[int, ...]
    coface_coords: tuple[int, ...]
    coface_shape: tuple[int, ...]
    r_grad: int
    external: bool

    @property
    def kind(self) -> str:
        return "external" if self.external else "internal"


@dataclasses.dataclass
class GOManifold:
    """A computed GO-manifold: its base, its trajectories, and its terminal slice."""

    pair: GOPair
    delta: float
    epsilon: float
    base: np.ndarray            # (m, N) sample of Nul_{n_o}^delta
    trajectories: list          # list of (k, N) arrays, one per base point
    terminal: np.ndarray        # (m, N) terminal slice
    hit_times: np.ndarray       # (m,) first hitting times
    directed_wall: float        # the n_g value defining the directed wall
    control_region: list        # Q^delta(xi, xi') as intervals

    @property
    def reached_wall(self) -> bool:
        """External case: every base point must reach the directed wall."""
        return bool(np.all(np.isfinite(self.hit_times)))

    @property
    def dimension(self) -> int:
        return self.base.shape[1]


# ---------------------------------------------------------------------------
# admissible spatial perturbation
# ---------------------------------------------------------------------------


def max_admissible_delta(system: RampSystem) -> float:
    """The supremum of ``delta`` allowed by ``eq:delta-restriction``.

    ``delta < Xi_n(mu)/2`` for every top cell ``mu`` and every ``n``, where
    ``Xi_n(mu)`` is the length of ``I_n(b(mu))`` -- for a top cell, the plateau
    between consecutive ramp windows.
    """
    best = float("inf")
    for box in itertools.product(*[range(k + 1) for k in system.num_thresholds]):
        ones = [1] * system.dim
        for n in range(system.dim):
            best = min(best, system.interval_length(list(box), ones, n) / 2.0)
    return best


def delta_intervals(system: RampSystem, coords, shape_bits, delta: float):
    """``I_n^delta(xi)`` (eq:Idelta)."""
    out = []
    for n in range(system.dim):
        k = coords[n]
        if shape_bits[n] == 0:  # n in J_i
            out.append(
                (
                    system.theta(n, k) - system.width(n, k) - delta,
                    system.theta(n, k) + system.width(n, k) + delta,
                )
            )
        else:  # n in J_e
            out.append(
                (
                    system.theta(n, k) + system.width(n, k) + delta,
                    system.theta(n, k + 1) - system.width(n, k + 1) - delta,
                )
            )
    return out


def control_region(system: RampSystem, pair: GOPair, delta: float):
    """``Q^delta(xi, xi')`` (eq:Qdelta-pair) as a list of intervals."""
    immutable = {pair.n_grad, pair.n_opaque} if _is_codim_two(pair, system.dim) else {pair.n_opaque}

    def q(coords, shape_bits):
        plain = [system.interval(list(coords), list(shape_bits), n) for n in range(system.dim)]
        shrunk = delta_intervals(system, coords, shape_bits, delta)
        return [plain[n] if n in immutable else shrunk[n] for n in range(system.dim)]

    a = q(pair.face_coords, pair.face_shape)
    b = q(pair.coface_coords, pair.coface_shape)
    return [(max(a[n][0], b[n][0]), min(a[n][1], b[n][1])) for n in range(system.dim)]


def _is_codim_two(pair: GOPair, dim: int) -> bool:
    return sum(pair.face_shape) == dim - 2


# ---------------------------------------------------------------------------
# locating the pairs
# ---------------------------------------------------------------------------


def go_pairs(
    system: RampSystem, *, spec: Spec = PAPER, outer: str = "global_bound"
) -> tuple[list[GOPair], SpecCubicalBlowupGraph]:
    """Every codimension-two indecisive-drift pair, with its pruning side."""
    labelling, num_thresholds = system.wall_labelling(outer=outer)
    stg = SpecCubicalBlowupGraph(
        labelling=labelling, num_thresholds=num_thresholds, spec=spec, level=2
    )
    cc = stg.cubical_complex
    dim = stg.dim
    out: list[GOPair] = []

    for face, coface, n_grad, n_opaque in stg.indecisive_drift_pairs(
        codim=(dim - 2, dim - 1)
    ):
        face_coords = tuple(cc.coordinates(face))
        coface_coords = tuple(cc.coordinates(coface))
        fs = cc.cell_shape(face)
        cs = cc.cell_shape(coface)
        face_shape = tuple(1 if fs & (1 << n) else 0 for n in range(dim))
        coface_shape = tuple(1 if cs & (1 << n) else 0 for n in range(dim))

        values = stg.rook_value_set(coface, n_grad)
        if len(values) != 1:
            # defn:GO-manifold requires R_{n_g}(xi') to be a singleton.
            continue
        r_grad = next(iter(values))
        if r_grad not in (-1, 1):
            continue

        # external: xi' not in F_2(xi); internal: xi not in F_2(xi')
        external = not stg.has_edge(face, coface)
        internal = not stg.has_edge(coface, face)
        if not (external or internal):
            continue

        out.append(
            GOPair(
                face=face,
                coface=coface,
                n_grad=n_grad,
                n_opaque=n_opaque,
                face_coords=face_coords,
                face_shape=face_shape,
                coface_coords=coface_coords,
                coface_shape=coface_shape,
                r_grad=r_grad,
                external=external,
            )
        )
    return out, stg


# ---------------------------------------------------------------------------
# the perturbed field
# ---------------------------------------------------------------------------


def perturbed_field(system: RampSystem, pair: GOPair, epsilon: float):
    """``eq:rampSysPerturbed``: a planar field in ``(n_o, n_g)``, frozen elsewhere."""
    n_o, n_g, r = pair.n_opaque, pair.n_grad, pair.r_grad

    def F_eps(x):
        x = np.asarray(x, dtype=float)
        e = system.production(x)
        out = np.zeros_like(x)
        out[n_o] = -system.gamma[n_o] * x[n_o] + e[n_o]
        out[n_g] = (-system.gamma[n_g] - r * epsilon) * x[n_g] + e[n_g]
        return out

    return F_eps


def max_admissible_epsilon(
    system: RampSystem, pair: GOPair, delta: float, *, samples: int = 7
) -> float:
    """The largest ``eps`` with ``r_{n_g} F_{eps,n_g} > 0`` on the control region.

    ``defn:admissible-GO-perturbation`` (GlobalDynR2v9.tex:1161-1180) requires
    that sign condition on ``Q^delta(xi', I)`` in the external case (on
    ``Q^delta(xi, I)`` in the internal one).  Since

        r F_{eps,n_g} = r F_{n_g} - eps x_{n_g},

    the condition is ``eps < r F_{n_g}(x) / x_{n_g}`` throughout the region, so
    the supremum is the minimum of that ratio.  A larger ``eps`` reverses the
    ``n_g``-motion and the sweep never reaches the directed wall.
    """
    immutable = {pair.n_grad, pair.n_opaque} if _is_codim_two(pair, system.dim) else {pair.n_opaque}
    coords = pair.coface_coords if pair.external else pair.face_coords
    shape = pair.coface_shape if pair.external else pair.face_shape

    plain = [system.interval(list(coords), list(shape), n) for n in range(system.dim)]
    shrunk = delta_intervals(system, coords, shape, delta)
    region = [plain[n] if n in immutable else shrunk[n] for n in range(system.dim)]

    grids = []
    for lo, hi in region:
        if hi <= lo:
            grids.append([0.5 * (lo + hi)])
        else:
            grids.append(list(np.linspace(lo, hi, samples)))

    n_g, r = pair.n_grad, pair.r_grad
    best = float("inf")
    for point in itertools.product(*grids):
        x = np.array(point, dtype=float)
        if x[n_g] <= 0:
            continue
        best = min(best, r * float(system.vector_field(x)[n_g]) / x[n_g])
    return best


def directed_wall_value(system: RampSystem, pair: GOPair) -> float:
    """``theta_{n_o,n_g,j} + r_{n_g} h_{n_o,n_g,j}``.

    The threshold of coordinate ``n_g`` that the GO-pair activates is the one at
    the face's ``n_g``-coordinate; the directed wall is the endpoint of its ramp
    window on the side ``r_{n_g}`` points to.
    """
    k = pair.face_coords[pair.n_grad]
    return system.theta(pair.n_grad, k) + pair.r_grad * system.width(pair.n_grad, k)


# ---------------------------------------------------------------------------
# the delta-nullcline base
# ---------------------------------------------------------------------------


def nullcline_base(
    system: RampSystem, pair: GOPair, delta: float, *, samples: int = 9
) -> np.ndarray:
    """Sample ``Nul_{n_o}^delta(xi, xi')`` as a graph in the ``n_g`` direction.

    ``Q^delta(xi, xi')`` is degenerate in ``n_o`` (the ramp window and the
    adjacent plateau meet in one point), so the base is cut out of the remaining
    ``N-1`` coordinates by one equation and is a graph over the ``N-2``
    coordinates other than ``n_g`` (``cor:nullcline-graph-delta``).
    """
    Q = control_region(system, pair, delta)
    n_o, n_g = pair.n_opaque, pair.n_grad
    x_no = 0.5 * (Q[n_o][0] + Q[n_o][1])

    free = [n for n in range(system.dim) if n not in (n_o, n_g)]
    grids = []
    for n in free:
        lo, hi = Q[n]
        if hi <= lo:
            grids.append([0.5 * (lo + hi)])
        else:
            grids.append(
                [lo + (hi - lo) * (i + 1) / (samples + 1) for i in range(samples)]
            )

    lo_g, hi_g = Q[n_g]
    points = []
    for combo in itertools.product(*grids) if free else [()]:
        base = np.zeros(system.dim)
        base[n_o] = x_no
        for n, v in zip(free, combo):
            base[n] = v

        def f(t: float) -> float:
            y = base.copy()
            y[n_g] = t
            return float(-system.gamma[n_o] * y[n_o] + system.production(y)[n_o])

        a, b = f(lo_g), f(hi_g)
        if a == 0.0:
            root = lo_g
        elif b == 0.0:
            root = hi_g
        elif a * b > 0:
            continue  # no crossing on this fibre
        else:
            root = brentq(f, lo_g, hi_g, xtol=1e-13, rtol=1e-13)
        point = base.copy()
        point[n_g] = root
        points.append(point)

    return np.array(points) if points else np.zeros((0, system.dim))


# ---------------------------------------------------------------------------
# the GO-manifold
# ---------------------------------------------------------------------------


def build_go_manifold(
    system: RampSystem,
    pair: GOPair,
    *,
    delta: float | None = None,
    epsilon: float | None = None,
    epsilon_fraction: float = 0.5,
    base_samples: int = 9,
    steps: int = 120,
    max_time: float = 200.0,
) -> GOManifold:
    """Construct ``M_{delta,eps}(xi, xi')`` of ``eq:GO-manifold``.

    With ``epsilon=None`` an admissible perturbation is chosen automatically as
    ``epsilon_fraction`` times the supremum allowed by
    ``defn:admissible-GO-perturbation``.
    """
    from scipy.integrate import solve_ivp

    if delta is None:
        delta = 0.25 * max_admissible_delta(system)
    if epsilon is None:
        epsilon = epsilon_fraction * max_admissible_epsilon(system, pair, delta)

    Q = control_region(system, pair, delta)
    base = nullcline_base(system, pair, delta, samples=base_samples)
    F_eps = perturbed_field(system, pair, epsilon)
    wall = directed_wall_value(system, pair)
    n_g = pair.n_grad

    def event(_t, x):
        return x[n_g] - wall

    event.terminal = True
    event.direction = float(pair.r_grad)

    trajectories = []
    terminal = []
    hit_times = []
    for x0 in base:
        sol = solve_ivp(
            lambda t, x: F_eps(x),
            (0.0, max_time),
            x0,
            events=event,
            dense_output=True,
            rtol=1e-9,
            atol=1e-12,
            max_step=max(2 * system.sorted_widths[n_g][0], 1e-3),
        )
        if sol.t_events[0].size:
            T = float(sol.t_events[0][0])
        else:
            T = float("inf")
        stop = T if np.isfinite(T) else float(sol.t[-1])
        ts = np.linspace(0.0, stop, steps)
        trajectories.append(sol.sol(ts).T)
        terminal.append(sol.sol(stop))
        hit_times.append(T)

    return GOManifold(
        pair=pair,
        delta=delta,
        epsilon=epsilon,
        base=base,
        trajectories=trajectories,
        terminal=np.array(terminal) if terminal else np.zeros((0, system.dim)),
        hit_times=np.array(hit_times),
        directed_wall=wall,
        control_region=Q,
    )


# ---------------------------------------------------------------------------
# checks against the R2 statements
# ---------------------------------------------------------------------------


def crossing_sign(system: RampSystem, manifold: GOManifold) -> dict:
    """``prop:GO-crossing-sign``: ``<F, z> = eps x_{n_g} |z_{n_g}| > 0`` off the base.

    ``lem:GO-product`` (GlobalDynR2v9.tex:780-805) says the transverse
    coordinate directions ``W_og`` are tangent to the sheet, so a normal ``z``
    has no transverse components; it also annihilates ``F_eps``, which sweeps
    the sheet.  Hence only the ``(n_o, n_g)`` block of

        F - F_eps = r_{n_g} eps x_{n_g} e^{(n_g)}

    contributes, giving ``<F, z> = r_{n_g} eps x_{n_g} z_{n_g}``.  Since
    ``r_{n_g} z_{n_g} = |z_{n_g}|``, that is positive.

    Checked numerically: within the ``(n_o, n_g)`` block the deviation must be
    exactly ``r eps x_{n_g}`` in the ``n_g`` slot and zero in the ``n_o`` slot.
    The transverse components of ``F`` are irrelevant here, and are reported
    separately as evidence for ``lem:GO-product``.
    """
    n_g, n_o, r = manifold.pair.n_grad, manifold.pair.n_opaque, manifold.pair.r_grad
    F_eps = perturbed_field(system, manifold.pair, manifold.epsilon)
    worst = float("inf")
    residual_no = 0.0
    residual_ng = 0.0
    checked = 0
    for traj in manifold.trajectories:
        for x in traj[1:]:  # skip the base point itself
            deviation = system.vector_field(x) - F_eps(x)
            expected = r * manifold.epsilon * float(x[n_g])
            residual_no = max(residual_no, abs(float(deviation[n_o])))
            residual_ng = max(residual_ng, abs(float(deviation[n_g]) - expected))
            worst = min(worst, r * float(deviation[n_g]))
            checked += 1
    return {
        "checked": checked,
        "min_signed_deviation": worst,
        "residual_n_o": residual_no,
        "residual_n_g": residual_ng,
        "holds": (
            checked > 0
            and worst > 0
            and residual_no < 1e-9
            and residual_ng < 1e-9
        ),
    }


def product_structure(system: RampSystem, manifold: GOManifold) -> dict:
    """``lem:GO-product``: the transverse directions are tangent to the sheet.

    In codimension two every ``k`` outside ``{n_o, n_g}`` lies in a gap on both
    sides of the shared wall, so every ramp with source ``k`` is constant on the
    control rectangle and neither the base nor its flow-out depends on ``x_k``.
    Checked by confirming the trajectories never move a transverse coordinate,
    and that the base really is a product in those coordinates.
    """
    n_g, n_o = manifold.pair.n_grad, manifold.pair.n_opaque
    transverse = [n for n in range(system.dim) if n not in (n_o, n_g)]
    drift = 0.0
    for traj in manifold.trajectories:
        for n in transverse:
            drift = max(drift, float(np.max(np.abs(traj[:, n] - traj[0, n]))))

    # The base should be a graph whose n_g value is independent of the
    # transverse coordinates.
    spread = 0.0
    if manifold.base.size and transverse:
        spread = float(np.max(manifold.base[:, n_g]) - np.min(manifold.base[:, n_g]))
    return {
        "transverse_drift": drift,
        "base_n_g_spread": spread,
        "holds": drift < 1e-9 and spread < 1e-9,
    }


def transversality(system: RampSystem, manifold: GOManifold) -> dict:
    """The true field crosses the sheet away from its base."""
    result = crossing_sign(system, manifold)
    return {
        "transverse_away_from_base": result["holds"],
        "min_signed_deviation": result["min_signed_deviation"],
        "samples": result["checked"],
    }


def base_is_on_nullcline(system: RampSystem, manifold: GOManifold) -> float:
    """Largest ``|F_{n_o}|`` over the sampled base; should be numerically zero."""
    n_o = manifold.pair.n_opaque
    if manifold.base.size == 0:
        return 0.0
    return float(
        max(abs(float(system.vector_field(x)[n_o])) for x in manifold.base)
    )


def summarise(system: RampSystem, *, spec: Spec = PAPER, **kw) -> list[dict]:
    """Build every codimension-two GO-manifold and report the R2 checks."""
    pairs, _stg = go_pairs(system, spec=spec)
    rows = []
    for pair in pairs:
        manifold = build_go_manifold(system, pair, **kw)
        rows.append(
            {
                "face": pair.face,
                "coface": pair.coface,
                "n_grad": pair.n_grad,
                "n_opaque": pair.n_opaque,
                "r_grad": pair.r_grad,
                "kind": pair.kind,
                "base_points": int(manifold.base.shape[0]),
                "base_residual": base_is_on_nullcline(system, manifold),
                "reached_directed_wall": manifold.reached_wall,
                "max_hit_time": float(np.max(manifold.hit_times))
                if manifold.hit_times.size
                else None,
                "directed_wall": manifold.directed_wall,
                **transversality(system, manifold),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# orange manifolds (defn:orange-manifold, GlobalDynR2v9.tex:1181-1313)
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class OrangeTriple:
    """A triple of ``D_3(rook)`` (``defn:triple-indecisive``).

    ``(xi, xi')`` is externally pruned and ``(xi, xi'')`` internally pruned, both
    with the same lower cell ``xi``; the orange connector mediates between the
    two resulting GO faces.
    """

    face: int
    external_coface: int
    internal_coface: int
    n_grad: int
    n_opaque: int
    r_grad: int
    face_coords: tuple
    face_shape: tuple
    p_opaque: int


@dataclasses.dataclass
class OrangeManifold:
    triple: OrangeTriple
    eta: float
    alpha: float
    sigma: "np.ndarray"
    anchor_g: "np.ndarray"
    anchor_o: "np.ndarray"

    def points(self, samples: int = 12) -> "np.ndarray":
        """The sheet, as an ``(m, samples, N)`` grid."""
        ts = np.linspace(0.0, 1.0, samples)
        return np.array(
            [np.array([(1 - t) * a + t * b for t in ts])
             for a, b in zip(self.anchor_g, self.anchor_o)]
        )


def band_endpoint(system: RampSystem, coords, shape_bits, n: int, side: int) -> float:
    """``Theta_n^{side}(xi)``: the ``side`` endpoint of ``I_n(b(xi))``."""
    lo, hi = system.interval(list(coords), list(shape_bits), n)
    return hi if side == 1 else lo


def d3_triples(stg: SpecCubicalBlowupGraph) -> list[OrangeTriple]:
    """``D_3(rook)``: pairs sharing a lower cell that F_2 prunes in opposite ways.

    ``defn:triple-indecisive`` (CombinatorialDynamics.tex:827-833):

        D_3 = { (xi, xi', xi'') : (xi,xi'), (xi,xi'') in D(Phi),
                                  xi' not in F_2(xi), xi not in F_2(xi'') }.
    """
    cc = stg.cubical_complex
    dim = stg.dim
    by_face: dict = {}
    for face, coface, n_grad, n_opaque in stg.indecisive_drift_pairs():
        by_face.setdefault(face, []).append((coface, n_grad, n_opaque))

    out = []
    for face, entries in by_face.items():
        external = [e for e in entries if not stg.has_edge(face, e[0])]
        internal = [e for e in entries if not stg.has_edge(e[0], face)]
        if not external or not internal:
            continue
        face_coords = tuple(cc.coordinates(face))
        fs = cc.cell_shape(face)
        face_shape = tuple(1 if fs & (1 << n) else 0 for n in range(dim))
        for coface_e, n_grad, n_opaque in external:
            values = stg.rook_value_set(coface_e, n_grad)
            if len(values) != 1:
                continue
            r_grad = next(iter(values))
            for coface_i, _ng2, _no2 in internal:
                if coface_i == coface_e:
                    continue
                internal_coords = cc.coordinates(coface_i)
                # p_{n_o}(xi, xi''): -1 when xi is the left n_o-wall of xi''
                p_opaque = -1 if face_coords[n_opaque] == internal_coords[n_opaque] else 1
                out.append(
                    OrangeTriple(
                        face=face,
                        external_coface=coface_e,
                        internal_coface=coface_i,
                        n_grad=n_grad,
                        n_opaque=n_opaque,
                        r_grad=r_grad,
                        face_coords=face_coords,
                        face_shape=face_shape,
                        p_opaque=p_opaque,
                    )
                )
    return out


def build_orange_manifold(
    system: RampSystem,
    triple: OrangeTriple,
    *,
    delta=None,
    eta=None,
    alpha=None,
    samples: int = 7,
) -> OrangeManifold:
    """``O_{delta,eps,eta,alpha}(Sigma; xi, xi'')`` of ``eq:orange-manifold``.

    The codimension-two flat is

        Pi(xi, xi'') = H_{n_g}^{r_{n_g}}(xi) cap H_{n_o}^{p_{n_o}(xi,xi'')}(xi),

    ``Sigma`` is a compact disk inside ``Pi cap closure(Q^delta(xi, xi''))``, and
    the manifold is the straight-line interpolation between the two anchors
    ``Sigma + r_{n_g} eta e^{(n_g)}`` and ``Sigma + p_{n_o} alpha e^{(n_o)}``.
    """
    if delta is None:
        delta = 0.25 * max_admissible_delta(system)
    n_g, n_o = triple.n_grad, triple.n_opaque

    coords, shape = list(triple.face_coords), list(triple.face_shape)
    pi_g = band_endpoint(system, coords, shape, n_g, triple.r_grad)
    pi_o = band_endpoint(system, coords, shape, n_o, triple.p_opaque)

    region = delta_intervals(system, coords, shape, delta)
    free = [n for n in range(system.dim) if n not in (n_g, n_o)]
    grids = []
    for n in free:
        lo, hi = region[n]
        grids.append([0.5 * (lo + hi)] if hi <= lo else list(np.linspace(lo, hi, samples)))

    sigma = []
    for combo in itertools.product(*grids) if free else [()]:
        point = np.zeros(system.dim)
        point[n_g] = pi_g
        point[n_o] = pi_o
        for n, v in zip(free, combo):
            point[n] = v
        sigma.append(point)
    sigma = np.array(sigma)

    if eta is None:
        eta = 0.5 * system.width(n_g, coords[n_g])
    if alpha is None:
        alpha = 0.5 * system.width(n_o, coords[n_o])

    e_g = np.zeros(system.dim)
    e_g[n_g] = 1.0
    e_o = np.zeros(system.dim)
    e_o[n_o] = 1.0
    return OrangeManifold(
        triple=triple,
        eta=eta,
        alpha=alpha,
        sigma=sigma,
        anchor_g=sigma + triple.r_grad * eta * e_g,
        anchor_o=sigma + triple.p_opaque * alpha * e_o,
    )
