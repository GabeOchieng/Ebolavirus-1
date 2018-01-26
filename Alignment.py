import sys
from Bio import SeqIO
from Bio import AlignIO
sys.setrecursionlimit(10000000)

"""
    Created by Mohsen Naghipourfar on 1/26/18.
    Email : mn7697np@gmail.com
"""

marburg_genes = []
for seq_record in SeqIO.parse("./Data/Marburg_Genes.fasta", "fasta"):
    marburg_genes.append(seq_record)



