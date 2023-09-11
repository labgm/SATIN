#!/usr/bin/ python

# -*- coding: utf-8 -*-
""" 

author: carlos willian

data:03/05/2023
proposito: Converter arquivos no formato genbank (.gbk, , .gbff, .gb, ...) para o formato ptt
input: genbank
output: txt file (tsv)
"""

from Bio import SeqIO
import sys
from Bio.SeqFeature import SimpleLocation, CompoundLocation

def convert_gbktoptt(file_gbk: str, output_file: str):

    output_file = output_file.rstrip('/')

    gb_record = None
    with open(file_gbk) as handle:
        for gb_record_ in SeqIO.parse(handle, "gb"):
            gb_record = gb_record_ #Para um multigbk sera selecionado somente a primeira anotacao
            break

    with open(output_file, "w") as handle:
        # handle.write("# This is a converted PTT file\n")
        handle.write("#Start    End    Strand    Length    Gene    locus_tag    Product\n")

        for feature in gb_record.features:

            if feature.type == "CDS":

                    if isinstance(feature.location, CompoundLocation): # If compound location CDS

                        # For Each parts CompoundLocation
                        for part in feature.location.parts:
                            
                            start = part.start.position + 1
                            end = part.end.position
                            strand = part.strand
                            length = end - start + 1

                            if strand == 1:
                                strand_char = "+"
                            else:
                                strand_char = "-"

                            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                            product = feature.qualifiers.get("product", [""])[0]
                            gene = feature.qualifiers.get("gene", [""])[0]
                            handle.write(f"{start}\t{end}\t{strand_char}\t{length}\t{gene}\t{locus_tag}\t{product}\n")


                    elif isinstance(feature.location, SimpleLocation): # If simple location CDS
                   
                        start = feature.location.start.position + 1
                        end = feature.location.end.position
                        strand = feature.location.strand
                        length = end - start + 1
                    
                        if strand == 1:
                            strand_char = "+"
                        else:
                            strand_char = "-"

                        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                        product = feature.qualifiers.get("product", [""])[0]
                        gene = feature.qualifiers.get("gene", [""])[0]
                        handle.write(f"{start}\t{end}\t{strand_char}\t{length}\t{gene}\t{locus_tag}\t{product}\n")


if __name__=='__main__':
    
    file_gbk = str(sys.argv[1])
    output_file = str(sys.argv[2])

    convert_gbktoptt(file_gbk, output_file)