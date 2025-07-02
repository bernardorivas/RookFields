#!/usr/bin/env python3
"""
Repressilator network analysis using DSGRN_utils
Generates JSON database with Morse graph analysis
"""

import DSGRN
import DSGRN_utils
from DSGRN_utils.SaveDatabaseJSON_CubicalBlowup import save_morse_graph_database_json
import json
import time
import sys

def main():
    # Define the repressilator network specification
    net_spec = """1 : ~3 :E
                  2 : ~1 :E 
                  3 : ~2 :E """
    network = DSGRN.Network(net_spec)
    
    # Create parameter graph
    parameter_graph = DSGRN.ParameterGraph(network)
    
    # Set parameter index (repressilator has only 1 parameter)
    par_index = 0
    
    # Output filename for JSON database
    output_filename = "repressilator_morse_database.json"
    
    print(f"\nGenerating Morse graph database for parameter {par_index}...")
    start_time = time.time()
    
    try:
        # Generate and save the database
        save_morse_graph_database_json(
            network=network,
            database_fname=output_filename,
            param_indices=[par_index],
            level=3  # Default level as used in ConleyMorseGraph
        )
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        print(f"Successfully generated database: {output_filename}")
        print(f"Time elapsed: {elapsed_time:.2f} seconds")
        
        # Optionally print some statistics about the generated file
        with open(output_filename, 'r') as f:
            data = json.load(f)
            if 'dynamics_database' in data and len(data['dynamics_database']) > 0:
                dynamics = data['dynamics_database'][0]
                num_morse_nodes = len(dynamics['morse_graph'])
                num_morse_sets = len(dynamics['morse_sets'])
                print(f"\nMorse graph has {num_morse_nodes} nodes")
                print(f"Found {num_morse_sets} Morse sets")
        
    except Exception as e:
        import traceback
        print(f"Error generating database: {e}")
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()