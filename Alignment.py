import sys
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/26/18.
    Email : mn7697np@gmail.com
"""

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

for genome in ebolavirus_genomes:
    f = open('./Output/' + genome.name + '.txt', "w")
    start = 0
    for gene in marburg_genes:
        len_gene = len(gene.seq)
        end = (start + len_gene * 3) if len(genome) > (start + len_gene * 3) else len(genome)
        gene_str = str(gene.seq)
        genome_str = str(genome.seq)[start: start + len_gene * 3]
        alignments = pairwise2.align.localmd(gene_str, genome_str, 1, -1, -1, -1, 0, 0)
        final_alignment = alignments[0]     # align1, align2, score, begin, end
        begin_idx = final_alignment[3]
        end_idx = final_alignment[4]
        f.write(gene.name + " :" + str(start + begin_idx) + " - " + str(start + end_idx) + "\n")
        start = end - len_gene
    f.close()

