"""

description: Script que extrai a sequência fasta do arquivo gkb utilizando a lib BioPython.

"""

import sys
from Bio import SeqIO

def extract_fasta(input_file_gbk: str, output_file_fasta: str):
    SeqIO.convert(input_file_gbk, "genbank", output_file_fasta, "fasta")
    return output_file_fasta

if __name__ == "__main__":

    input_file_gbk = sys.argv[1]
    output_file_fasta = sys.argv[2]

    extract_fasta(input_file_gbk, output_file_fasta)