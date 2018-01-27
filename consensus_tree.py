import sys
from Bio import Phylo
from Bio import pairwise2
from Bio.Phylo.Consensus import *
from Bio import SeqIO
import phylogenetic_tree, Alignment

sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/27/18.
    Email : mn7697np@gmail.com
"""

NJ_trees = []
UPGMA_trees = []


def merge_all_trees():
    phylogenetic_tree.construct_trees_for_all_genes_without_marburg()
    NJ_trees = phylogenetic_tree.NJ_trees  # All NJ Trees in a list
    UPGMA_trees = phylogenetic_tree.UPGMA_trees  # All UPGMA Trees in a list
    NJ_tree = majority_consensus(NJ_trees,
                                 0.4)  # Merge NJ Trees using Majority Consensus Algorithm (0.4 is best for practical)
    UPGMA_tree = majority_consensus(UPGMA_trees,
                                    0.4)  # Merge UPGMA Trees using Majority Consensus Algorithm (0.4 is best for practical)
    phylogenetic_tree.save_tree(UPGMA_tree, 'UPGMA_Merged')  # Draw merged UPGMA Tree --> it is not a good merge
    phylogenetic_tree.save_tree(NJ_tree, 'NJ_Merged')  # Draw merged NJ Tree

    # final_tree = majority_consensus([NJ_tree, UPGMA_tree], 0.4)  # Merge UPGMA && NJ Trees (Not Recommended!)
    # Phylo.draw_graphviz(final_tree) # Draw Final Tree


def align_all_ebola_genomes():  # All all ebola genomes to each other
    edm = [[0 for i in range(5)] for j in range(5)]  # New edit matrix
    Alignment.read_data()
    all_genomes = Alignment.ebolavirus_genomes
    g1_id = 0
    for genome1 in all_genomes:
        g2_id = 0
        for genome2 in all_genomes:
            if genome1.name != genome2.name and g2_id > g1_id:
                print('Aligning {0} with {1}'.format(genome1.name, genome2.name))
                alignments = pairwise2.align.globalms(genome1, genome2, 1, -1, -1, 0)  # Biopython package
                alignment = alignments[0]  # first alignment
                score = alignment[2]  # score of alignment
                # a, b, score = Alignment.global_alignment(genome1.seq, genome2.seq) # Global Alignment
                edm[g1_id][g2_id] = score
                edm[g2_id][g1_id] = score
            g2_id += 1
        g1_id += 1
    print("genomes aligned!")
    Alignment.save_edit_matrix("all_ebola_genomes", edm)  # Save edit matrix to file
    phylogenetic_tree.construct_tree("all_ebola_genomes", algorithm="UPGMA")  # Construct Tree
    phylogenetic_tree.construct_tree("all_ebola_genomes", algorithm="NJ")  # Construct Tree


def construct_trees_for_all_genes_with_marburg():
    for gene_name in Alignment.gene_names:  # For all genes
        NJ_trees.append(
            phylogenetic_tree.construct_tree(gene_name, with_marburg=2, algorithm="NJ"))  # Construct NJ Tree
        UPGMA_trees.append(
            phylogenetic_tree.construct_tree(gene_name, with_marburg=2, algorithm="UPGMA"))  # Construct UPGMA Tree


if __name__ == '__main__':
    merge_all_trees()
    align_all_ebola_genomes()
    construct_trees_for_all_genes_with_marburg()
