"""The two specs must agree with the two modes of `DSGRN_utils`.

`code/DSGRN_utils` now implements the manuscript's definitions by default, with
`legacy=True` restoring the superseded behaviour.  Two guarantees are pinned
here:

* ``LEGACY`` reproduces ``CubicalBlowupGraph(..., legacy=True)`` bit-for-bit, so
  every figure and number produced before the alignment still comes out
  identical;
* ``PAPER`` reproduces the new default, so ``rookfields`` and plain
  ``DSGRN_utils`` cannot drift apart.
"""

from __future__ import annotations

import pytest
from DSGRN_utils.CubicalBlowupGraph import CubicalBlowupGraph

from rookfields import LEGACY, PAPER, SpecCubicalBlowupGraph
from rookfields import networks

FAMILIES = ["toggle", "repressilator", "cycle3", "cycle4"]
LEVELS = [0, 1, 2, 3, 4]


def sample(name, limit=25):
    total = networks.parameter_graph(name).size()
    if total <= limit:
        return list(range(total))
    step = max(1, total // limit)
    return list(range(0, total, step))[:limit]


@pytest.mark.parametrize("name", FAMILIES)
@pytest.mark.parametrize("level", LEVELS)
def test_legacy_matches_dsgrn_utils_legacy_mode(name, level):
    for index in sample(name):
        p = networks.parameter(name, index)
        expected = set(
            CubicalBlowupGraph(p, level=level, legacy=True).digraph.edges()
        )
        actual = set(
            SpecCubicalBlowupGraph(p, spec=LEGACY, level=level).digraph.edges()
        )
        assert actual == expected, f"{name}[{index}] level={level}"


@pytest.mark.parametrize("name", FAMILIES)
@pytest.mark.parametrize("level", LEVELS)
def test_paper_matches_dsgrn_utils_default(name, level):
    """`DSGRN_utils` now implements the manuscript's definitions by default."""
    for index in sample(name):
        p = networks.parameter(name, index)
        expected = set(CubicalBlowupGraph(p, level=level).digraph.edges())
        actual = set(
            SpecCubicalBlowupGraph(p, spec=PAPER, level=level).digraph.edges()
        )
        assert actual == expected, f"{name}[{index}] level={level}"


@pytest.mark.parametrize("name", FAMILIES)
def test_legacy_matches_on_vertices_too(name):
    """Vertex sets must match as well, not only edges."""
    for index in sample(name, limit=10):
        p = networks.parameter(name, index)
        expected = set(
            CubicalBlowupGraph(p, level=4, legacy=True).digraph.vertices()
        )
        actual = set(SpecCubicalBlowupGraph(p, spec=LEGACY, level=4).digraph.vertices())
        assert actual == expected
