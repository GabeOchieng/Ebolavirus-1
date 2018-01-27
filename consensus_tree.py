import sys
from Bio import Phylo
from Bio.Phylo.Consensus import *
import phylogenetic_tree

sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/27/18.
    Email : mn7697np@gmail.com
"""

""" ---PERFORMS SEMI-GLOBAL ALIGNMENT TO THE GIVEN SEQUENCES """

NJ_trees = phylogenetic_tree.NJ_trees  # All NJ Trees in a list
UPGMA_trees = phylogenetic_tree.UPGMA_trees  # All UPGMA Trees in a list
NJ_trees = adam_consensus(NJ_trees)  # Merge NJ Trees using Adam Consensus Algorithm
UPGMA_trees = adam_consensus(UPGMA_trees)  # Merge UPGMA Trees using Adam Consensus Algorithm
Phylo.draw(UPGMA_trees)  # Draw merged UPGMA Tree
Phylo.draw(NJ_trees)  # Draw merged NJ Tree
