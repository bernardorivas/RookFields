#!/usr/bin/env python3
"""
DSGRN-GINsim Integration Example

Network example:
- X1 : (X1+X2)    # X1 is activated by X1 OR X2
- X2 : (~X1)(X2)  # X2 is activated by X2 AND NOT X1
"""

import ginsim
import re
from itertools import product
from dataclasses import dataclass


@dataclass
class DSGRNNetwork:
    """Parsed DSGRN network representation."""
    nodes: list[str]
    edges: list[tuple[str, str, int]]  # (source, target, sign)
    rules: dict[str, str]


class DSGRNGINsimConverter:
    """Convert between DSGRN and GINsim formats."""
    
    def parse_dsgrn(self, spec: str) -> DSGRNNetwork:
        """Parse DSGRN specification."""
        nodes = []
        edges = []
        rules = {}
        
        for line in spec.strip().split('\n'):
            if ':' not in line:
                continue
                
            target, regulation = line.split(':', 1)
            target = target.strip()
            regulation = regulation.strip()
            
            nodes.append(target)
            rules[target] = regulation
            
            # Extract regulators
            for match in re.finditer(r'(~?)(\w+)', regulation):
                negated, source = match.groups()
                if source != target and source[0].isupper():  # Assume uppercase = node name
                    sign = -1 if negated else 1
                    edges.append((source, target, sign))
        
        return DSGRNNetwork(
            nodes=list(set(nodes)),
            edges=list(set(edges)),
            rules=rules
        )
    
    def dsgrn_to_ginsim(self, spec: str) -> ginsim.RegulatoryGraph:
        """Convert DSGRN to GINsim model."""
        network = self.parse_dsgrn(spec)
        
        # Create GINsim graph
        lrg = ginsim.RegulatoryGraph()
        
        # Add nodes
        for node in network.nodes:
            lrg.add_node(node, 1)  # Boolean by default
        
        # Add edges
        for source, target, sign in network.edges:
            lrg.add_edge(source, target, sign=sign)
        
        return lrg
    
    def ginsim_to_dsgrn(self, lrg: ginsim.RegulatoryGraph) -> str:
        """Convert GINsim to DSGRN format."""
        lines = []
        
        for node in lrg.nodes():
            # Get regulators
            activators = []
            inhibitors = []
            
            for edge in lrg.edges():
                if edge[1] == node:
                    if lrg.edges[edge]['sign'] > 0:
                        activators.append(edge[0])
                    else:
                        inhibitors.append(edge[0])
            
            # Build expression
            terms = []
            if activators:
                terms.append(f"({'+'.join(activators)})")
            for inh in inhibitors:
                terms.append(f"(~{inh})")
            
            expr = ''.join(terms) or '1'
            lines.append(f"{node} : {expr}")
        
        return '\n'.join(lines)


class BooleanAnalyzer:
    """Analyze Boolean network dynamics."""
    
    def __init__(self, network: DSGRNNetwork):
        self.network = network
        self.node_order = sorted(network.nodes)
        self.state_space_size = 2 ** len(self.node_order)
    
    def evaluate_rule(self, rule: str, state: dict[str, int]) -> int:
        """Evaluate a DSGRN rule given current state."""
        # Simple evaluation for our example rules
        # X1 : (X1+X2) means X1 OR X2
        # X2 : (~X1)(X2) means (NOT X1) AND X2
        
        # Replace variables with values
        expr = rule
        for node, value in state.items():
            expr = expr.replace(f'~{node}', str(1 - value))
            expr = expr.replace(node, str(value))
        
        # Evaluate OR operations
        if '+' in expr:
            # (X1+X2) pattern
            match = re.search(r'\(([01+]+)\)', expr)
            if match:
                or_expr = match.group(1)
                result = 1 if '1' in or_expr else 0
                expr = expr.replace(match.group(0), str(result))
        
        # Evaluate AND operations (multiplication)
        while '(' in expr and ')' in expr:
            expr = re.sub(r'\(([01])\)', r'\1', expr)
        
        # Multiple terms mean AND
        if re.match(r'^[01]+$', expr):
            return 0 if '0' in expr else 1
        
        return int(expr) if expr in '01' else 0
    
    def compute_dynamics(self) -> dict:
        """Compute full dynamics of the Boolean network."""
        results = {
            'states': [],
            'transitions': [],
            'stable_states': [],
            'state_graph': {}
        }
        
        # Generate all possible states
        for state_tuple in product([0, 1], repeat=len(self.node_order)):
            state = dict(zip(self.node_order, state_tuple))
            
            # Compute next state
            next_state = {}
            for node in self.node_order:
                if node in self.network.rules:
                    next_state[node] = self.evaluate_rule(self.network.rules[node], state)
                else:
                    next_state[node] = state[node]
            
            # Store results
            state_info = {
                'state': state,
                'next_state': next_state,
                'is_stable': state == next_state
            }
            results['states'].append(state_info)
            
            if state == next_state:
                results['stable_states'].append(state)
            
            # Build state graph (for asynchronous dynamics)
            state_key = tuple(state[n] for n in self.node_order)
            results['state_graph'][state_key] = []
            
            # Add asynchronous transitions (one variable changes at a time)
            for node in self.node_order:
                if next_state[node] != state[node]:
                    async_next = state.copy()
                    async_next[node] = next_state[node]
                    async_key = tuple(async_next[n] for n in self.node_order)
                    results['state_graph'][state_key].append(async_key)
        
        return results
    
    def find_attractors(self, dynamics: dict) -> list[list[tuple]]:
        """Find all attractors (stable states and cycles)."""
        attractors = []
        visited = set()
        
        # Stable states are single-state attractors
        for stable in dynamics['stable_states']:
            state_tuple = tuple(stable[n] for n in self.node_order)
            attractors.append([state_tuple])
            visited.add(state_tuple)
        
        # Find cycles using DFS
        state_graph = dynamics['state_graph']
        
        def find_cycle_from(start: tuple, path: list[tuple]) -> list[tuple] | None:
            if start in path:
                # Found a cycle
                cycle_start = path.index(start)
                return path[cycle_start:]
            
            if start in visited or not state_graph[start]:
                return None
            
            for next_state in state_graph[start]:
                if cycle := find_cycle_from(next_state, path + [start]):
                    return cycle
            
            return None
        
        # Check each unvisited state for cycles
        for state in state_graph:
            if state not in visited:
                if cycle := find_cycle_from(state, []):
                    attractors.append(cycle)
                    visited.update(cycle)
        
        return attractors


def main():
    """Run the DSGRN-GINsim integration demo."""
    
    # Define DSGRN network
    dsgrn_spec = """
    X1 : (X1+X2)
    X2 : (~X1)(X2)
    """
    
    print("DSGRN Network Specification:")
    print(dsgrn_spec)
    print("\n" + "="*50 + "\n")
    
    # Convert to GINsim
    converter = DSGRNGINsimConverter()
    network = converter.parse_dsgrn(dsgrn_spec)
    lrg = converter.dsgrn_to_ginsim(dsgrn_spec)
    
    print("Parsed Network Structure:")
    print(f"Nodes: {network.nodes}")
    print(f"Edges: {network.edges}")
    print(f"Rules: {network.rules}")
    print("\n" + "="*50 + "\n")
    
    # Analyze dynamics
    analyzer = BooleanAnalyzer(network)
    dynamics = analyzer.compute_dynamics()
    
    print("Boolean Network Analysis:")
    print(f"State space size: {analyzer.state_space_size}")
    print(f"\nTruth Table:")
    print("X1 X2 | X1' X2'")
    print("------|--------")
    
    for state_info in dynamics['states']:
        s = state_info['state']
        ns = state_info['next_state']
        stable = "*" if state_info['is_stable'] else " "
        print(f"{s['X1']}  {s['X2']}  | {ns['X1']}   {ns['X2']}  {stable}")
    
    print(f"\nStable states: {len(dynamics['stable_states'])}")
    for i, stable in enumerate(dynamics['stable_states']):
        print(f"  {i+1}. X1={stable['X1']}, X2={stable['X2']}")
    
    # Find attractors
    attractors = analyzer.find_attractors(dynamics)
    print(f"\nAttractors: {len(attractors)}")
    for i, attractor in enumerate(attractors):
        if len(attractor) == 1:
            print(f"  {i+1}. Fixed point: {attractor[0]}")
        else:
            print(f"  {i+1}. Cycle of length {len(attractor)}: {attractor}")
    
    # Convert back to DSGRN
    reconstructed = converter.ginsim_to_dsgrn(lrg)
    print(f"\nReconstructed DSGRN specification:")
    print(reconstructed)
    
    # Show asynchronous state graph
    print("\n" + "="*50 + "\n")
    print("Asynchronous State Transition Graph:")
    for state, successors in dynamics['state_graph'].items():
        if successors:
            print(f"{state} -> {successors}")
    
    return dynamics


if __name__ == "__main__":
    # Note: This uses basic ginsim functionality
    # Some features may require the full COLOMOTO docker environment
    try:
        results = main()
    except ImportError as e:
        print(f"Error: {e}")
        print("Install ginsim with: pip install ginsim")