### ActiveInputEdges.py
### MIT LICENSE 2025 Marcio Gameiro

import DSGRN

def active_input(partial_order, input_index):
    """Check if input with index input_index is active in partial order"""
    # Get number of outputs from partial order
    n_outputs = len([p for p in partial_order if p < 0])
    # Get threshold indices in partial order
    thres_indices = [partial_order.index(-(k + 1)) for k in range(n_outputs)]
    # Get input polynomials with lower value at input_index
    polys_lower = [p for p in partial_order if p >= 0 and (p & (1 << input_index)) == 0]
    for p_lower in polys_lower:
        # Get corresponding p upper
        p_upper = p_lower ^ (1 << input_index)
        # Get sorted indices of p_lower and p_upper in partial order
        i1, i2 = sorted([partial_order.index(p_lower), partial_order.index(p_upper)])
        # Check if there is a threshold in between them
        if any(i1 < t_index < i2 for t_index in thres_indices):
            return True
    return False

def active_input_indices(partial_order, n_inputs):
    """Get indices of active inputs in partial order"""
    active_indices = [index for index in range(n_inputs) if active_input(partial_order, index)]
    return active_indices

def number_active_edges(par_index, parameter_graph):
    """Get number of active edges for parameter index"""
    parameter = parameter_graph.parameter(par_index)
    network = parameter_graph.network()
    num_active_edges = 0
    for node in range(network.size()):
        # Get number of inputs and outputs
        n_inputs = len(network.inputs(node))
        n_outputs = len(network.outputs(node))
        # Treat the no out edge case as one out edge
        n_outputs += 1 if n_outputs == 0 else 0
        # Get factor graph hex code
        hex_code = parameter.logic()[node].hex()
        # Compute partial order corresponding to hex code
        partial_order = DSGRN.hex2partial(hex_code, n_inputs, n_outputs)
        num_active_edges += len(active_input_indices(partial_order, n_inputs))
    return num_active_edges

def paths_from_node(start_node, parameter_graph, max_weight=None):
    """Get all paths with increasing node weights starting at start_node,
       where the weight of a node is the number of active edges.
    """
    dim = parameter_graph.network().size()
    # Set max weight if not given
    if max_weight == None:
        max_weight = sum(max(len(parameter_graph.network().inputs(k)), 1) for k in range(dim))
    # Check if target node, i. e., if weight >= max weight
    target = lambda node: number_active_edges(node, parameter_graph) >= max_weight
    # target = lambda node: number_active_edges(node, parameter_graph) == max_weight
    # Check if allowed next node, i. e., if weight > current weight
    allowed = lambda node, current_weight: number_active_edges(node, parameter_graph) > current_weight
    # Stack with tuples (current_node, current_path)
    paths_stack = [(start_node, [start_node])]
    paths = []
    while paths_stack:
        current_node, current_path = paths_stack.pop()
        if target(current_node):
            paths.append(current_path)
        current_weight = number_active_edges(current_node, parameter_graph)
        adjacencies = parameter_graph.adjacencies(current_node)
        for node in [node for node in adjacencies if allowed(node, current_weight)]:
            paths_stack.append((node, current_path + [node]))
    return paths

def get_all_paths(starting_nodes, parameter_graph, max_weight=None):
    """Get all paths starting at a node in starting_nodes"""
    paths = {}
    for start_node in starting_nodes:
        paths[start_node] = paths_from_node(start_node, parameter_graph, max_weight=max_weight)
    return paths
