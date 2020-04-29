#!/usr/local/bin/python3
#
#program will go through the .txt file and get the
#genome position associated with each section, which is equivelent to
#the position for each gene (basically)
#run as: python3 ~/CODE/dN_dS_info_with_genome_pos.py (ref taxa number) > out
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import pandas as pd #importing pandas
import numpy as np #importing Numpy
import csv
import os #importing operating system commands
import glob #importing glob which helps read files
from numpy import nan
from Bio import SeqIO
from Bio import AlignIO
from functools import partial
from Bio.Phylo.PAML import codeml #import the biopython package that
#will read the codeml file
# FOR DEBUGGING ONLY pdb.set_trace()
######################################################################

#ref taxa form command line arg
ref_taxa = int(sys.argv[1])

#below will open all the txt/mafft files and get the info we need from
#them
block_array = []#list of blocks
gene_array = []#list of genes
sec_array = []#list of section of gene
sec_array2 = []#list of section of gene
genome_pos_array = []#list of genome position of gene
block_starts = []
block_ends = []
block_midpoint = []
taxa_names = []
rev_comp = []
#keeping track of if the gene is a complement or not
complement_array = []
for filepath in glob.glob(os.path.join('', 'Block*_good_sec_*.txt')):
    tax_count = 0
    with open(filepath) as f:
        ###########################################################
        #parse the file name to get what block and gene and section it
        #is from
        bac_name = os.path.basename(filepath)
        #store 'a' which is the '_' character so we can sub it out later
        a = re.compile("\_")
        file_name = a.sub(" ",bac_name)
        #split the ID so we separate the actual ID from the rest of
        #the file name
        split_id = re.split(" ",file_name)
        #append the block and stuff to various arrays
        block_array.append(split_id[0])
        gene_array.append(split_id[4])
        sec_name = split_id[6]
        sec_name2 = split_id[9]
        #removing the file extention in the section name
        b = re.compile("\.txt")
        sec_name_mod = b.sub("",sec_name)
        sec_array.append(sec_name_mod)
        sec_name2_mod = b.sub("",sec_name2)
        sec_array2.append(sec_name2_mod)
        ############################################################
        #now to actually read the mafft file and get the info
        for line in f:
            if re.match('^>', line):
                tax_count += 1
                a = re.compile("\:\n")
                tmp_id = a.sub("",line)
                a = re.compile("\:")
                tmp_id = a.sub(" ",line)
                b = re.compile("\-")
                tmp_id = b.sub(" ",tmp_id)
                b = re.compile(">")
                tmp_id = b.sub("",tmp_id)
                b = re.compile("\n")
                tmp_id = b.sub("",tmp_id)
                tmp_line = re.split(" ",tmp_id)
                #making sure it matches the ref taxa number
                if (tax_count == ref_taxa):
                    #append start and stop and taxa name to list so we can use it later
                    taxa_names.append(tmp_line[0])
                    block_starts.append(tmp_line[1])
                    block_ends.append(tmp_line[2])
                    #finding the midpoint of the gene
                    midpoint = int(tmp_line[2]) + int(tmp_line[1])
                    midpoint = midpoint / 2
                    block_midpoint.append(midpoint)
                    #storing info about if each seq is a reverse complement
                    if (len(tmp_line) > 3):
                        rev_comp.append(1)
                        tmp_rev_comp = tmp_line[3]
                    else:
                        rev_comp.append(0)

###############################################################
#open the gene name file with the starts and stops and gene names
gene_name = []
gene_name_start = []
gene_name_end = []
for filepath in glob.glob(os.path.join('', '*_coding_start_end.txt')):
    with open(filepath) as tsv:
        col_count = 0
        for line in csv.reader(tsv, dialect="excel-tab"): #You can also use delimiter="\t" rather than giving a dialect.
            gene_name_start.append(int(line[0]))
            gene_name_end.append(int(line[1]))
            gene_name.append(line[3])
###############################################################
#go through each midpoint and see if it is in a gene and then save it
#as a new gene name col
gene_name_blocks = []
for b in range(0,len(block_midpoint)):
    for c in range(0, len(gene_name_start)):
        if (block_midpoint[b] >= gene_name_start[c] and block_midpoint[b] <= gene_name_end[c]):
            gene_name_blocks.append(gene_name[c])
            break
###############################################################
#now to print out the data in a nice dataframe
#block  gene    sec dS  dN  omega   sec_len
print("block\tgene\tsec\tsec_two\tstart\tend\tmidpoint\tgene_name")
for s in range(0,len(block_array)):
    print(block_array[s],gene_array[s],sec_array[s],sec_array2[s],block_starts[s],block_ends[s],block_midpoint[s],gene_name_blocks[s],sep='\t')

