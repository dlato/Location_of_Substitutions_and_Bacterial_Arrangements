#!/usr/local/bin/python3
#
# Program will go through the genbank file and will obtain the locus
# ID for each gene as well as the start and end positions
# run as "get_gene_pos_locus.py *.gbk > *gene_locus_location" 
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import Bio
import numpy as np
from Bio import SeqIO
# FOR DEBUGGING ONLY pdb.set_trace()

faa_filename = sys.argv[2]
output_handle = open(faa_filename, "w")
for gb_record in SeqIO.parse(open(sys.argv[1],"r"), "genbank") :
    genes = gb_record.features
#    print(genes)
    for seq_feature in gb_record.features :
        pseudo = 0
        if seq_feature.type=="CDS":
            if (seq_feature.qualifiers):
                if ('pseudo' in seq_feature.qualifiers):
                    continue
                #gene is not a pseudogene
                else:
                    end = seq_feature.location.end
                    start = seq_feature.location.start
                    start = start + 1
                    tmp = end + start
                    midpoint = int(tmp/2)
                    if seq_feature.strand == 1:
                        cod_val = 0
                    else:
                        cod_val = 1
                    output_handle.write("%s\t%s\t%s\t%s\t%s\n" % (
                           start,
                           seq_feature.location.end,
                           midpoint,
        #                   seq_feature.qualifiers['gene'][0],
                           seq_feature.qualifiers['locus_tag'][0],
                           cod_val))

output_handle.write("0\t0\t0\t0\t0\n")

output_handle.close()
