import sys
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import SeqIO
import phylogenetic_tree, Alignment

sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/27/18.
    Email : mn7697np@gmail.com
"""

""" ---PERFORMS SEMI-GLOBAL ALIGNMENT TO THE GIVEN SEQUENCES """

NJ_trees = phylogenetic_tree.NJ_trees  # All NJ Trees in a list
UPGMA_trees = phylogenetic_tree.UPGMA_trees  # All UPGMA Trees in a list
NJ_tree = adam_consensus(NJ_trees)  # Merge NJ Trees using Adam Consensus Algorithm
UPGMA_tree = adam_consensus(UPGMA_trees)  # Merge UPGMA Trees using Adam Consensus Algorithm
Phylo.draw(UPGMA_tree)  # Draw merged UPGMA Tree
Phylo.draw(NJ_tree)  # Draw merged NJ Tree

final_tree = adam_consensus([NJ_tree, UPGMA_tree]) # Merge UPGMA && NJ Trees
Phylo.draw(final_tree) # Draw Final Tree


marburg_genome = SeqIO.read("./Data/Marburg_genome.fasta", "fasta")
all_genomes = Alignment.ebolavirus_genomes + [marburg_genome]
Alignment.read_genes(all_genomes)
Alignment.global_align()
for gene_name in Alignment.gene_names:  # For all genes
    NJ_trees.append(phylogenetic_tree.construct_tree(gene_name,names=2, type="NJ"))  # Construct NJ Tree
    UPGMA_trees.append(phylogenetic_tree.construct_tree(gene_name, names=2, type="UPGMA"))  # Construct UPGMA Tree



