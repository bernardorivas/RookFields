"""Realising a DSGRN parameter node by an explicit ramp system.

``defn:DSGRN_ramp`` (RampSystemsv4.tex:955-970) replaces each switching function
by a type-1 ramp with the same plateaus.  ``prop:DSGRN_wall_constant`` (:975)
then says the induced wall labelling is constant on ``R(k) x H_0``, so it must
agree with the one DSGRN computes combinatorially.

That is what makes the geometrization pictures possible: it turns a purely
combinatorial parameter node into a genuine vector field on the geometrized
complex.
"""

from __future__ import annotations

import pytest

from rookfields import PAPER, conley_morse_graph
from rookfields.networks import parameter_graph, realises_dsgrn_labelling
from rookfields.ramp import RampSystem

CASES = [
    ("toggle", [0, 3, 4, 8]),
    ("N2_a", [47, 752, 900]),
    ("N2_b", [974, 1109]),
    ("cycle3", [0, 10, 100]),
    ("N3_B", [2_472_286, 2_472_287]),
]


@pytest.mark.parametrize("name,indices", CASES, ids=[c[0] for c in CASES])
def test_ramp_system_realises_the_dsgrn_labelling(name, indices):
    """`prop:DSGRN_wall_constant`, checked on the sampled point."""
    for index in indices:
        _spec, realised = realises_dsgrn_labelling(name, index, 0.05)
        assert realised, f"{name}[{index}] not realised by a ramp system"


@pytest.mark.parametrize("name,indices", CASES, ids=[c[0] for c in CASES])
def test_realised_system_gives_the_same_morse_graph(name, indices):
    """The ramp realisation reproduces the node's Morse graph, not just its labelling."""
    for index in indices:
        spec, realised = realises_dsgrn_labelling(name, index, 0.05)
        assert realised
        system = RampSystem(spec)
        labelling, K = system.wall_labelling()

        from_ramp = conley_morse_graph(
            labelling=labelling, num_thresholds=K, spec=PAPER, level=3
        )
        from_dsgrn = conley_morse_graph(
            parameter_graph(name).parameter(index), spec=PAPER, level=3
        )
        assert from_ramp.conley_indices == from_dsgrn.conley_indices
        assert from_ramp.edges == from_dsgrn.edges


def test_product_interaction_is_respected():
    """`E_n` must use the network's factor structure, not a bare sum.

    `2 : (~1) 2` is a product of two factors; summing the plateaus instead
    gives the wrong `M^E_n`, hence the wrong global bound `GB_n`, hence the
    wrong label on the outermost wall.
    """
    from rookfields.networks import RampSpec, network, parameter

    net = network("N2_a")
    spec = RampSpec.from_dsgrn(net, parameter("N2_a", 752), 0.02)
    assert spec.logic == ((( 0, 1),), ((0,), (1,)))

    nu, theta, _h = spec._tables()
    # node 2 is a product, so its maximum production is a product of factor maxima
    expected = max(nu[1][0]) * max(nu[1][1])
    assert spec.max_production[1] == pytest.approx(expected)
    # node 1 is a single additive factor
    assert spec.max_production[0] == pytest.approx(max(nu[0][0]) + max(nu[0][1]))
    del theta
