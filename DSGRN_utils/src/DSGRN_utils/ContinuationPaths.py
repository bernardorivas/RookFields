### ContinuationPaths.py
### MIT LICENSE 2025 Marcio Gameiro

import DSGRN_utils

def is_morse_set(par_node, morse_set, parameter_graph, level):
    """Check if morse_set is a Morse set for the Morse graph at parameter par_node"""
    # Compute Morse graph
    parameter = parameter_graph.parameter(par_node)
    morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter, level=level)
    # Get list of stable Morse nodes
    stable_nodes = [v for v in morse_graph.vertices() if not morse_graph.adjacencies(v)]
    for vert in stable_nodes:
        # Get Morse node (vertex index)
        morse_node = morse_graph.vertex_label(vert)[0]
        # Get set of cells in the Morse set corresponding to the Morse node morse_node
        vertices, edges = DSGRN_utils.MorseComponent(morse_graph, stg, graded_complex, morse_node=morse_node)
        # Check if Morse sets are equal
        if {v for v in vertices} == morse_set:
            return True
    return False

def continuation_path(path, parameter_graph, level=3):
    """Check if a path is a valid continuation path, i. e., check
       if the Morse set in the initial node of the path is a Morse
       set for all nodes in the path.
    """
    # Compute Morse graph for initial parameter in the path
    parameter = parameter_graph.parameter(path[0])
    morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter, level=level)
    # Check if Morse graph contains a single node
    if len(morse_graph.vertices()) != 1:
        return False
    # Get set of cells in the Morse set corresponding to the Morse node
    vertices, edges = DSGRN_utils.MorseComponent(morse_graph, stg, graded_complex, morse_node=0)
    # # Check if the Morse set consists of a single cell
    # if len(vertices) != 1:
    #     return False
    # Get set of cells in Morse set at initial parameter
    initial_morse_set = {v for v in vertices}
    for par_node in path[1:]:
        # Check if initial_morse_set is a Morse set at par_node
        if not is_morse_set(par_node, initial_morse_set, parameter_graph, level):
            return False
    return True

