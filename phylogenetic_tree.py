import sys
import pandas as pd
import numpy as np
import pylab
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import matplotlib
import matplotlib.pyplot as plt
import Alignment

sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/26/18.
    Email : mn7697np@gmail.com
"""
NJ_trees = []
UPGMA_trees = []


def convert_tu_lower_triangular(data_frame):  # Convert data_frame to a lower triangular matrix
    data_frame = data_frame.values.T.tolist()  # Convert data_frame to list of lists
    for i in range(5):
        data_frame[i] = data_frame[i][:i + 1]  # Remove unused data from data_frame --> Convert to lower triangular
    return data_frame


def construct_tree(gene_name, type='UPGMA'):  # Construct Tree with specific type (Default = UPGMA)
    edit_matrix = pd.read_csv("./Output/edit_matrices/" + gene_name + ".csv", header=None)  # read edit matrix file
    constructor = DistanceTreeConstructor()  # Create a tree constructor object
    edit_matrix = convert_tu_lower_triangular(edit_matrix)  # Convert Edit Distance matrix to lower triangular
    distance_matrix = DistanceMatrix(names=['Bundibugyo', 'Reston', 'Sudan', 'TaiForest', 'Zaire'], matrix=edit_matrix)
    if type == 'NJ':  # Neighbor-Joining Alogrithm
        tree = constructor.nj(distance_matrix)
    else:  # UPGMA Algorithm
        tree = constructor.upgma(distance_matrix)
    save_tree(tree, type + '_' + gene_name)  # Save Tree into a file
    return tree


def save_tree(tree, filename):
    Phylo.draw(tree)  # Draw the tree
    plt.title(filename)  # set Title for figure
    plt.savefig('./Output/images/' + filename + '.png', dpi=100)  # Save tree in an image
    plt.close()  # Close the figure


for gene_name in Alignment.gene_names:  # For all genes
    NJ_trees.append(construct_tree(gene_name, "NJ"))  # Construct NJ Tree
    UPGMA_trees.append(construct_tree(gene_name, "UPGMA"))  # Construct UPGMA Tree
