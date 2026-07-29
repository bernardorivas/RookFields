"""Rook Fields: paper-faithful reference implementation and simulation tools.

Companion code for ``paper/Rook_Field_Paper_v2``.  Three jobs:

1. compute the combinatorial multivalued maps ``F_0 .. F_4`` under both the
   ``DSGRN_utils`` implementation and the live manuscript definitions, and diff
   them (:mod:`rookfields.spec`, :mod:`rookfields.blowup`, :mod:`rookfields.audit`);
2. reproduce every computational example in the monograph and compare against
   the claimed output (:mod:`rookfields.examples`);
3. construct, render, and numerically validate the geometrizations of Part III
   (:mod:`rookfields.ramp`, :mod:`rookfields.geometrization`,
   :mod:`rookfields.alignment`, :mod:`rookfields.r2_manifolds`,
   :mod:`rookfields.r3_manifolds`).
"""

from .blowup import SpecCubicalBlowupGraph
from .pipeline import MorseResult, conley_morse_graph
from .spec import ALL_SPECS, ISOLATED, LEGACY, PAPER, Spec, by_name

__all__ = [
    "ALL_SPECS",
    "ISOLATED",
    "LEGACY",
    "PAPER",
    "MorseResult",
    "Spec",
    "SpecCubicalBlowupGraph",
    "by_name",
    "conley_morse_graph",
]

__version__ = "0.1.0"
