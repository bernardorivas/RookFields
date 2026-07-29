"""Rendering the geometrized complex, the Morse sets, and the R2/R3 manifolds.

Unlike ``DSGRN_utils.PlotMorseSets``, which draws the cell complex in index
coordinates with artificial alternating widths (``real_coord``,
PlotMorseSets.py:132), everything here is drawn in the *real* phase-space
coordinates supplied by the rectangular geometrization.  That is what makes it
possible to overlay the vector field, the nullclines, genuine trajectories and
the swept GO-manifolds on the same axes.
"""

from .phase2d import plot_phase_portrait_2d
from .phase3d import (
    plot_go_manifold_3d,
    plot_go_manifold_in_context,
    plot_morse_sets_3d,
)

__all__ = [
    "plot_go_manifold_in_context",
    "plot_go_manifold_3d",
    "plot_morse_sets_3d",
    "plot_phase_portrait_2d",
]
