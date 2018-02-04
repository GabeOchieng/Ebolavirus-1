import sys
import pandas as pd
import numpy as np
import pylab
from Bio import Phylo
from Bio import SeqIO
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
    size = len(data_frame[0])
    for i in range(size):
        data_frame[i] = data_frame[i][:i + 1]  # Remove unused data from data_frame --> Convert to lower triangular
    return data_frame


def construct_tree(gene_name, with_marburg=1, algorithm='UPGMA'):  # Construct Tree with specific type (Default = UPGMA)
    if with_marburg == 1:
        print('Constructing Tree with All Viruses without Marburg')
        filename = algorithm + '_' + gene_name
        names = ['Bundibugyo', 'Reston', 'Sudan', 'TaiForest', 'Zaire']
    else:
        print('Constructing {0}\'s Tree with All Viruses with Marburg'.format(gene_name))
        filename = algorithm + '_' + gene_name + '_with_Marburg'
        names = ['Bundibugyo', 'Reston', 'Sudan', 'TaiForest', 'Zaire', 'Marburg']
        marburg_genome = SeqIO.read("./Data/Marburg_genome.fasta", "fasta")
        Alignment.read_data()
        print('Aligning Genes for marburg_genome')
        gene_name += '_with_marburg'
        Alignment.read_genes(marburg_genome)
    print('Reading edit matrix and construct tree')
    edit_matrix = pd.read_csv("./Output/edit_matrices/" + gene_name + ".csv", header=None)  # read edit matrix file
    constructor = DistanceTreeConstructor()  # Create a tree constructor object
    edit_matrix = convert_tu_lower_triangular(edit_matrix)  # Convert Edit Distance matrix to lower triangular
    distance_matrix = DistanceMatrix(names=names, matrix=edit_matrix)
    if algorithm == 'NJ':  # Neighbor-Joining Alogrithm
        tree = constructor.nj(distance_matrix)
    else:  # UPGMA Algorithm
        tree = constructor.upgma(distance_matrix)
    save_tree(tree, filename)  # Save Tree into a file
    return tree


def save_tree(tree, filename):
    Phylo.draw_graphviz(tree, prog='dot')  # Draw the tree
    plt.title(filename)  # set Title for figure
    plt.savefig('./Output/images/' + filename + '.png', dpi=100)  # Save tree in an image
    plt.close()  # Close the figure


def construct_trees_for_all_genes_without_marburg():
    for gene_name in Alignment.gene_names:  # For all genes
        NJ_trees.append(construct_tree(gene_name, algorithm="NJ"))  # Construct NJ Tree
        UPGMA_trees.append(construct_tree(gene_name, algorithm="UPGMA"))  # Construct UPGMA Tree

if __name__ == '__main__':
    construct_trees_for_all_genes_without_marburg() # Section 3.1 in pdf
