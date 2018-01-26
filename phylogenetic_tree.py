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


def convert_tu_lower_triangular(data_frame):
    data_frame = data_frame.values.T.tolist()
    for i in range(5):
        data_frame[i] = data_frame[i][:i + 1]
    return data_frame


def construct_tree(gene_name, type='UPGMA'):
    edit_matrix = pd.read_csv("./Output/edit_matrices/" + gene_name + ".csv", header=None)
    constructor = DistanceTreeConstructor()
    edit_matrix = convert_tu_lower_triangular(edit_matrix)
    distance_matrix = DistanceMatrix(names=['Bundibugyo', 'Reston', 'Sudan', 'Tai', 'Zaire'], matrix=edit_matrix)
    if type == 'NJ':
        tree = constructor.nj(distance_matrix)
    else:
        tree = constructor.upgma(distance_matrix)
    save_tree(tree, type + '_' + gene_name)


def save_tree(tree, filename):
    Phylo.draw_graphviz(tree, prog='dot')
    plt.savefig('./Output/images/' + filename + '.png', dpi=100)
    plt.close()

for gene_name in Alignment.gene_names:
    construct_tree(gene_name, "NJ")
    construct_tree(gene_name, "UPGMA")
