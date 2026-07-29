### IsomorphismQuery.py
### MIT LICENSE 2024 Marcio Gameiro

import DSGRN
import DSGRN_utils
from collections import defaultdict

def morse_graph_signature(morse_graph):
    """Return a hashable invariant of a Morse graph.

    DSGRN.isomorphic_morse_graphs expects the string vertex labels this package
    used before the labels became lists, so it cannot be applied directly.  The
    signature below compares the same data it did: the multiset of Conley
    indices together with the reachability order on the labelled nodes.
    """
    index_of = {v: morse_graph.vertex_label(v)[0] for v in morse_graph.vertices()}
    conley = {index_of[v]: tuple(morse_graph.vertex_label(v)[2]) for v in morse_graph.vertices()}
    edges = frozenset((index_of[v], index_of[u])
                      for v in morse_graph.vertices()
                      for u in morse_graph.adjacencies(v))
    return (tuple(sorted(conley.items())), edges)


def IsomorphismQuery(network, param_indices=None, level=4):
    """Return a list of sets of parameters with isomorphic Morse graphs"""
    parameter_graph = DSGRN.ParameterGraph(network)
    if param_indices == None:
        param_indices = range(parameter_graph.size())
    # Isomorphism classes, keyed by Morse graph signature
    isomorphism_classes = defaultdict(set)
    for par_index in param_indices:
        parameter = parameter_graph.parameter(par_index)
        morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter, level=level)
        # A parameter belongs to exactly one class
        isomorphism_classes[morse_graph_signature(morse_graph)].add(par_index)
    return list(isomorphism_classes.values())
