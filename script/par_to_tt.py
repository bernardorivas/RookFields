## Convert a DSGRN parameter to a truth table

import DSGRN
import DSGRN_utils
import os
import random
import itertools

# Load the boolean network
script_dir = os.path.dirname(__file__)
# network_file = os.path.join(script_dir, "network.txt")
network_file = os.path.join(script_dir, "network2.txt")

# Define the network
network = DSGRN.Network(network_file)

# Draw the network
G_network = DSGRN.DrawGraph(network)
G_network.render(os.path.join(script_dir, "network"), format="png", view=False, cleanup=True)

# Define the parameter graph
parameter_graph = DSGRN.ParameterGraph(network)
print('Number of parameter nodes:', parameter_graph.size())

# Pick a parameter index
par_index = int(parameter_graph.size() * random.random())
# par_index = 42
print('Parameter index:', par_index)

# Define the dsgrn parameter
parameter = parameter_graph.parameter(par_index)

# From here on we should be able to compute everything we need for the truth table

# Get the labelling, dimension and number of thresholds
labelling = parameter.labelling()
dim = network.size()
num_thresholds = [len(network.outputs(i)) for i in range(dim)]

# Calculate place values (pv) for indexing into the labelling list.
# This helps convert multi-dimensional coordinates to a single list index.
pv = [1]
for k in num_thresholds:
    pv.append(pv[-1] * (k + 1))
print('pv:', pv)

# Build the truth table string
# Start with the header line (e.g., "variable_1 variable_2 ... variable_dim")
header = " ".join(network.name(i) for i in range(dim))
tt_lines = [header]

# Generate all possible binary inputs (corners of the hypercube)
# e.g., for dim=2, this will be [(0,0), (0,1), (1,0), (1,1)]
# e.g., for dim=3, this will be [(0,0,0), (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1), (1,1,0), (1,1,1)]
inputs = list(itertools.product([0, 1], repeat=dim))

# Iterate through each input combination to create a row in the truth table
for input_tuple in inputs:
    input_str = "".join(map(str, input_tuple))
    output_str = ""

    # For each node/dimension, determine its output value
    for i, input_val in enumerate(input_tuple):
        # a. Determine the coordinates of the extreme cell in the parameter space
        coords = [num_thresholds[j] if input_tuple[j] == 1 else 0 for j in range(dim)]
        
        # b. Decide which wall to check based on the input value
        # If input is 0, check the right wall (+1); if 1, check the left wall (-1).
        side_to_check = 1 if input_val == 0 else -1
        
        # c. Calculate the wall label directly
        label_index = sum(c * pv[k] for k, c in enumerate(coords))
        integer_label = labelling[label_index]
        
        bit_to_check = i if side_to_check == -1 else i + dim
        bitmask = 1 << bit_to_check
        is_absorbing = (integer_label & bitmask) != 0
        
        label = side_to_check if is_absorbing else -side_to_check
        
        # d. Map the label to the truth table output ('1' for ON, '0' for OFF)
        output_char = '1' if label == 1 else '0'
        output_str += output_char
        
    tt_lines.append(f"{input_str} {output_str}")
    
truth_table_content = "\n".join(tt_lines)

# Write the final string to the output file
# Generate output filename based on input network file
network_base = os.path.splitext(os.path.basename(network_file))[0]
output_filename = os.path.join(script_dir, f"{network_base}.tt")
with open(output_filename, "w") as f:
    f.write(truth_table_content)

print(f"\nTruth table successfully generated and saved to '{output_filename}'")
print("\n--- File Content ---")
print(truth_table_content)