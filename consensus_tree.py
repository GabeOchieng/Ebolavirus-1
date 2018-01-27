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



NJ_trees = phylogenetic_tree.NJ_trees
UPGMA_trees = phylogenetic_tree.UPGMA_trees
NJ_trees = adam_consensus(NJ_trees)
UPGMA_trees = adam_consensus(UPGMA_trees)
Phylo.draw(UPGMA_trees)
Phylo.draw(NJ_trees)