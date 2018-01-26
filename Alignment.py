import sys
from Bio import SeqIO
from Bio import pairwise2
import pandas as pd
import numpy as np


sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/26/18.
    Email : mn7697np@gmail.com
"""

marburg_genes = []
ebolavirus_genomes = []
all_genes = {}
edit_distance_matrices = [[[0 for i in range(5)] for j in range(5)] for k in range(7)]

def align_and_find_genes(genome):
    f = open('./Output/' + genome.name + '.csv', "w")
    start = 0
    for gene in marburg_genes:
        len_gene = len(gene.seq)
        end = (start + len_gene * 3) if len(genome) > (start + len_gene * 3) else len(genome)
        gene_str = str(gene.seq)
        genome_str = str(genome.seq)[start: start + len_gene * 3]
        alignments = pairwise2.align.localmd(gene_str, genome_str, 1, -1, -1, -1, 0, 0)
        final_alignment = alignments[0]  # align1, align2, score, begin, end
        begin_idx = final_alignment[3]
        end_idx = final_alignment[4]
        f.write(gene.name + "," + str(start + begin_idx) + "," + str(start + end_idx) + "\n")
        start = end - len_gene
    f.close()


def read_data():
    global ebolavirus_genomes, marburg_genes
    # read Marburg genes data
    marburg_genes = []
    for seq_record in SeqIO.parse("./Data/Marburg_Genes.fasta", "fasta"):
        marburg_genes.append(seq_record)

    # read Ebolaviruses genome data as seq_record object
    Bundibugyo_genome = SeqIO.read("./Data/Bundibugyo_genome.fasta", "fasta")
    Reston_genome = SeqIO.read("./Data/Reston_genome.fasta", "fasta")
    Sudan_genome = SeqIO.read("./Data/Sudan_genome.fasta", "fasta")
    TaiForest_genome = SeqIO.read("./Data/TaiForest_genome.fasta", "fasta")
    Zaire_genome = SeqIO.read("./Data/Zaire_genome.fasta", "fasta")
    ebolavirus_genomes = [Bundibugyo_genome, Reston_genome, Sudan_genome, TaiForest_genome, Zaire_genome]


def start_aligning():
    for genome in ebolavirus_genomes:
        align_and_find_genes(genome)
        break


if __name__ == '__main__':
    read_data()
    for gene in marburg_genes:
        i = 0
        genes = []
        for genome in ebolavirus_genomes:
            indices = pd.read_csv("./Output/" + genome.name + ".csv", header=None)
            begin_idx = int(indices.loc[i, 1])
            end_idx = int(indices.loc[i, 2])
            new_record = SeqIO.SeqRecord(genome.seq[begin_idx: end_idx])
            new_record.name = genome.name
            genes.append(new_record)
            i += 1
        all_genes[gene.name] = genes
    gene_id = 0
    for gene in all_genes.values():
        g1_id = 0
        for genome1 in gene:
            g2_id = 0
            for genome2 in gene:
                if genome1.seq != genome2.seq:
                    alignments = pairwise2.align.globalms(genome1, genome2, 0, -1, -1, -1)
                    alignment = alignments[0]
                    score = alignment[2]
                    edit_distance = -1 * score
                    edit_distance_matrices[gene_id][g1_id][g2_id] = edit_distance
                    edit_distance_matrices[gene_id][g2_id][g1_id] = edit_distance
                g2_id += 1
            g1_id += 1
        gene_id += 1
    i = 0
    for name, gene in all_genes.items():
        edit_matrix = np.array(edit_distance_matrices[i])
        np.savetxt("./Output/edit_matrices/" + name + ".csv", edit_matrix, delimiter=",", fmt='%d')
        i += 1