import DSGRN
import os
import itertools
import random

# --- Main execution ---

# 1. Define the network from the provided string
network_spec = "1 : 1 + 2 : B\n2 : (~1)2 : B"
network = DSGRN.Network(network_spec)

# 2. Define the parameter graph for the network
parameter_graph = DSGRN.ParameterGraph(network)
print(f'Number of parameter nodes: {parameter_graph.size()}')

# 3. Pick the specified parameter index
par_index = 42
print(f'Parameter index: {par_index}')

# 4. Get the specific DSGRN parameter from the graph
parameter = parameter_graph.parameter(par_index)

# --- Generate the truth table directly from the parameter ---

# 5. Extract essential information for calculation
labelling = parameter.labelling()
dim = network.size()
num_thresholds = [len(network.outputs(i)) for i in range(dim)]

# 6. Calculate place values (pv) for indexing into the labelling list.
# This helps convert multi-dimensional coordinates to a single list index.
pv = [1]
for k in num_thresholds:
    pv.append(pv[-1] * (k + 1))

# 7. Build the truth table string
# Start with the header line (e.g., "1 2")
header = " ".join(network.name(i) for i in range(dim))
tt_lines = [header]

# Generate all possible binary inputs (corners of the hypercube)
# e.g., for dim=2, this will be [(0,0), (0,1), (1,0), (1,1)]
inputs = list(itertools.product([0, 1], repeat=dim))

# 8. Iterate through each input combination to create a row in the truth table
for input_tuple in inputs:
    input_str = "".join(map(str, input_tuple))
    output_str = ""

    # Add detailed debugging for the (0,0) case
    if input_tuple == (0, 0):
        print("\n--- Detailed Calculation for Input (0,0) ---")

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
        
        # Print debugging information for the (0,0) case
        if input_tuple == (0, 0):
            print(f"\n  Node {i+1} (input bit {input_val}):")
            print(f"    Coordinates of extreme cell: {coords}")
            print(f"    Wall to check: {'Right' if side_to_check == 1 else 'Left'} (side = {side_to_check})")
            print(f"    Label Index for these coordinates: {label_index}")
            print(f"    Integer Label from labelling list: {integer_label} (binary: {bin(integer_label)})")
            print(f"    Bit to check for this wall: {bit_to_check}")
            print(f"    Bitmask: {bitmask} (binary: {bin(bitmask)})")
            print(f"    Is wall absorbing? (integer_label & bitmask != 0): {is_absorbing}")
            print(f"    Resulting wall label: {label}")
            print(f"    Final output character: '{output_char}'")

    if input_tuple == (0, 0):
        print("\n--- End of Detailed Calculation ---")

    tt_lines.append(f"{input_str} {output_str}")
    
truth_table_content = "\n".join(tt_lines)

# 9. Write the final string to the output file
output_filename = "output.tt"
with open(output_filename, "w") as f:
    f.write(truth_table_content)

print(f"\nTruth table successfully generated and saved to '{output_filename}'")
print("\n--- File Content ---")
print(truth_table_content)
