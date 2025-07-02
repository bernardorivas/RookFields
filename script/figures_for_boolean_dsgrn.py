#!/usr/bin/env python3
"""
Figures for Boolean/DSGRN Analysis

This script generates figures for analyzing boolean networks using DSGRN.
All figures are saved to the figures/ directory.
"""

import os
import time
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from collections import defaultdict

import DSGRN
import DSGRN_utils
import graphviz

# Global figure settings
FIGURE_SETTINGS = {
    'dpi': 300,
    'format': 'png',
    'bbox_inches': 'tight',
    'facecolor': 'white',
    'edgecolor': 'none'
}

# Default color list for Morse graphs and sets
DEFAULT_COLOR_LIST = ['#1f77b4', '#e6550d', '#31a354', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                      '#bcbd22', '#80b1d3', '#ffffb3', '#fccde5', '#b3de69', '#fdae6b', '#6a3d9a', '#c49c94',
                      '#fb8072', '#dbdb8d', '#bc80bd', '#ffed6f', '#637939', '#c5b0d5', '#636363', '#c7c7c7',
                      '#8dd3c7', '#b15928', '#e8cb32', '#9e9ac8', '#74c476', '#ff7f0e', '#9edae5', '#90d743',
                      '#e7969c', '#17becf', '#7b4173', '#8ca252', '#ad494a', '#8c6d31', '#a55194', '#00cc49']

# Color settings
MORSE_PLOT_SETTINGS = {
    'clist': DEFAULT_COLOR_LIST, 
    'self_arrow_clr': 'black',
    'arrow_clr': 'black',
    'plot_bdry_cells': True
}

# Output directory
OUTPUT_DIR = 'figures'

def save_current_figure(filename):
    """Save the current matplotlib figure with standard settings."""
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(filepath, **FIGURE_SETTINGS)
    print(f"Saved figure: {filepath}")
    plt.close()


def save_graphviz_figure(gv_source, filename):
    """Save a graphviz.Source object as a PNG file."""
    # Remove extension
    base_filename = filename.rsplit('.', 1)[0]
    filepath = os.path.join(OUTPUT_DIR, base_filename)
    
    # Render the graphviz source to PNG
    gv_source.render(filepath, format='png', cleanup=True)
    
    # Graphviz adds the extension, so the final file is at filepath.png
    print(f"Saved figure: {filepath}.png")

def get_morse_set_sizes(morse_graph, stg, graded_complex):
    """Get the size (number of cells) in each morse set."""
    # Get the grading value
    fringe_node_grade = graded_complex.value(stg.complex().size() - 1)
    
    # Get sets of cells with the same grading (same SCC)
    grading_cells = defaultdict(set)
    for cell in stg.digraph.vertices():
        val = graded_complex.value(cell)
        if val == fringe_node_grade:
            continue
        grading_cells[val].add(cell)
    
    # Map morse graph vertices to their grading values and count cells
    morse_set_sizes = {}
    for v in morse_graph.vertices():
        # The vertex v is the grading value
        if v in grading_cells:
            morse_set_sizes[v] = len(grading_cells[v])
        else:
            morse_set_sizes[v] = 0
    
    return morse_set_sizes


def get_custom_labels(morse_graph, morse_set_sizes):
    """Generate custom labels based on morse set size and Conley index."""
    custom_labels = {}
    
    for v in morse_graph.vertices():
        # Get the original label
        original_label = morse_graph.vertex_label(v)
        
        # Extract Conley index (part after ":")
        conley_part = original_label.split(':')[1].strip()
        # Parse the tuple
        conley_index = eval(conley_part)
        
        # Get morse set size
        size = morse_set_sizes.get(v, 0)
        
        # Apply labeling rules
        if size > 2:
            if conley_index == (1, 1, 0):
                new_label = "FC"  # Full Cycle
            elif conley_index == (0, 0, 0):
                # new_label = "XC"
                new_label = ""
            else:
                # new_label = "FP+FC"  # Partial Cycle
                new_label = ""
        else:
            new_label = "FP"  # Fixed Point (size == 1)
        
        custom_labels[original_label] = new_label
    
    return custom_labels


def apply_custom_labels(gv_source, custom_labels):
    """Apply custom labels to a graphviz source object."""
    # Modify the source string
    modified_source = gv_source.source
    
    for old_label, new_label in custom_labels.items():
        # Replace label="old_label" with label="new_label"
        modified_source = modified_source.replace(f'label="{old_label}"', f'label="{new_label}"')
    
    # Create new graphviz source with modified labels
    return graphviz.Source(modified_source)


def generate_figure_1():
    """Figure 1: Toggle Switch"""
    net_spec = """v1 : ~v2
                  v2 : ~v1 """
    
    network = DSGRN.Network(net_spec)
    parameter_graph = DSGRN.ParameterGraph(network)
    par_index = 4
    parameter = parameter_graph.parameter(par_index)

    morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter)
    
    # Get morse set sizes and create custom labels
    morse_set_sizes = get_morse_set_sizes(morse_graph, stg, graded_complex)
    custom_labels = get_custom_labels(morse_graph, morse_set_sizes)
    
    # Create morse graph with default labels
    gv_source = DSGRN_utils.PlotMorseGraph(morse_graph, clist=DEFAULT_COLOR_LIST)
    
    # Apply custom labels
    gv_source_custom = apply_custom_labels(gv_source, custom_labels)
    
    save_graphviz_figure(gv_source_custom, 'Example1_MorseGraph.png')
    
    # Plot Morse sets
    DSGRN_utils.PlotMorseSets(morse_graph, stg, graded_complex, **MORSE_PLOT_SETTINGS)
    save_current_figure('Example1_MorseSets.png')

def generate_figure_2():
    """Figure 2: Self-activation with repression"""
    net_spec = """v1 : v1+v2
                  v2 : ~v1 """
    
    network = DSGRN.Network(net_spec)
    parameter_graph = DSGRN.ParameterGraph(network)

    partial_orders = ['(p0, p1, p2, p3, t0, t1)', '(p0, t0, p1)']
    par_index = DSGRN.index_from_partial_orders(parameter_graph, partial_orders)
    parameter = parameter_graph.parameter(par_index)
    morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter)

    # Get morse set sizes and create custom labels
    morse_set_sizes = get_morse_set_sizes(morse_graph, stg, graded_complex)
    custom_labels = get_custom_labels(morse_graph, morse_set_sizes)
    
    # Create morse graph with default labels
    gv_source = DSGRN_utils.PlotMorseGraph(morse_graph, clist=DEFAULT_COLOR_LIST)
    
    # Apply custom labels
    gv_source_custom = apply_custom_labels(gv_source, custom_labels)
    
    save_graphviz_figure(gv_source_custom, 'Example2_MorseGraph.png')
    
    DSGRN_utils.PlotMorseSets(morse_graph, stg, graded_complex, **MORSE_PLOT_SETTINGS)
    save_current_figure('Example2_MorseSets.png')


def generate_figure_3():
    """Figure 3: Self-activation with repression"""
    net_spec = """v1 : v1+v2
                  v2 : ~v1 """
    
    network = DSGRN.Network(net_spec)
    parameter_graph = DSGRN.ParameterGraph(network)

    partial_orders = ['(p0, p2, t0, p1, t1, p3)', '(p0, t0, p1)']
    par_index = DSGRN.index_from_partial_orders(parameter_graph, partial_orders)
    parameter = parameter_graph.parameter(par_index)
    morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter)

    # Get morse set sizes and create custom labels
    morse_set_sizes = get_morse_set_sizes(morse_graph, stg, graded_complex)
    custom_labels = get_custom_labels(morse_graph, morse_set_sizes)
    
    # Create morse graph with default labels
    gv_source = DSGRN_utils.PlotMorseGraph(morse_graph, clist=DEFAULT_COLOR_LIST)
    
    # Apply custom labels
    gv_source_custom = apply_custom_labels(gv_source, custom_labels)
    
    save_graphviz_figure(gv_source_custom, 'Example3_MorseGraph.png')
    
    DSGRN_utils.PlotMorseSets(morse_graph, stg, graded_complex, **MORSE_PLOT_SETTINGS)
    save_current_figure('Example3_MorseSets.png')

def generate_figure_4():
    """Figure 4: Self-activation with AND logic"""
    
    net_spec = """v1 : v1+v2
                  v2 : (~v1)v2 """
    
    network = DSGRN.Network(net_spec)
    
    parameter_graph = DSGRN.ParameterGraph(network)
    
    configs = [
        {
            'name': 'config_1',
            'partial_orders': ['(p0, p1, p2, t0, t1, p3)', '(p0, p1, p2, t0, t1, p3)']
        },
        {
            'name': 'config_2',
            'partial_orders': ['(p0, p1, p2, t1, t0, p3)', '(p0, p1, p2, t0, t1, p3)']
        },
        {
            'name': 'config_3',
            'partial_orders': ['(p0, p1, p2, t0, t1, p3)', '(p0, p1, p2, t1, t0, p3)']
        },
        {
            'name': 'config_4',
            'partial_orders': ['(p0, p1, p2, t1, t0, p3)', '(p0, p1, p2, t1, t0, p3)']
        }
    ]
    
    for config in configs:
        par_index = DSGRN.index_from_partial_orders(parameter_graph, config['partial_orders'])
        parameter = parameter_graph.parameter(par_index)
        
        morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter)
        
        # Get morse set sizes and create custom labels
        morse_set_sizes = get_morse_set_sizes(morse_graph, stg, graded_complex)
        custom_labels = get_custom_labels(morse_graph, morse_set_sizes)
        
        # Plot Morse graph
        gv_source = DSGRN_utils.PlotMorseGraph(morse_graph, clist=DEFAULT_COLOR_LIST)
        
        # Apply custom labels
        gv_source_custom = apply_custom_labels(gv_source, custom_labels)
        
        save_graphviz_figure(gv_source_custom, f'Example4_{config["name"]}_MorseGraph.png')
        
        # Plot Morse sets
        DSGRN_utils.PlotMorseSets(morse_graph, stg, graded_complex, **MORSE_PLOT_SETTINGS)
        save_current_figure(f'Example4_{config["name"]}_MorseSets.png')


def main():
    # Generate figures
    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    generate_figure_4()
    

if __name__ == "__main__":
    main()