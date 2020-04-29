#!/usr/local/bin/python3
#
#program will go through the output of the .mlc file and get the
#necessary dN and dS info and put this, plus the gene number and
#sequence length into one big file
#run as: python3 ~/CODE/get_genes_from_mafft_to_align.py *.mafft > out
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


#below will open all the codeml files and get the info we need from
#them
block_array = []#list of blocks
gene_array = []#list of genes
sec_array = []#list of section of gene
sec_array2 = []#list of section of gene
sec_len_array = []#list of section lengths of gene
dN_array = []
dS_array = []
omega_array = []
#keeping track of if the gene is a complement or not
complement_array = []
for filepath in glob.glob(os.path.join('', 'Block*.mlc')):
    with open(filepath) as f:
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
        b = re.compile("\.mlc")
        sec_name_mod = b.sub("",sec_name)
        sec_array.append(sec_name_mod)
        sec_name2_mod = b.sub("",sec_name2)
        sec_array2.append(sec_name2_mod)
        ############################################################
        #now to actually read the codeml output file and get the info
        results = codeml.read(bac_name)
        #if the key is in the dictionary
        if ('NSsites' in results):
            try:
                dN_array.append(results['NSsites'][0]['parameters']['dN'])
            except KeyError:
                dN_array.append("NA")
            try:
                dS_array.append(results['NSsites'][0]['parameters']['dS'])
            except KeyError:
                dS_array.append("NA")
            try:
                omega_array.append(results['NSsites'][0]['parameters']['omega'])
            except KeyError:
                omega_array.append("NA")
        else:
            dN_array.append("NA")
            dS_array.append("NA")
            omega_array.append("NA")
        #############################################################
        #now getting the length of each section in bp that was used in the codeml
        #calculation (not including stop codons)
        tmp_sec_lens = []
        cwd = os.getcwd()
        for line in f:
            #dealing with the fact that bacillus has 7 seqs
            if re.findall('Bacillus', cwd):
                if re.match('^      7', line):
                    tmp_line = re.split(" ",line)
                    tmp_var = tmp_line[-1]
                    b = re.compile("\n")
                    tmp_var = b.sub("",tmp_var)
                    tmp_sec_lens.append(tmp_var)
#            if re.match('^      7', line):
#                tmp_line = re.split(" ",line)
#                tmp_var = tmp_line[-1]
#                b = re.compile("\n")
#                tmp_var = b.sub("",tmp_var)
#                tmp_sec_lens.append(tmp_var)
            if re.match('^      6', line):
                tmp_line = re.split(" ",line)
                tmp_var = tmp_line[-1]
                b = re.compile("\n")
                tmp_var = b.sub("",tmp_var)
                tmp_sec_lens.append(tmp_var)
                #dealing with the fact that strep now only has 3 seqs
            if re.match('^      3', line):
                tmp_line = re.split(" ",line)
                tmp_var = tmp_line[-1]
                b = re.compile("\n")
                tmp_var = b.sub("",tmp_var)
                tmp_sec_lens.append(tmp_var)
            if re.match('After deleting gaps. ', line):
                tmp_line = re.split(" ",line)
                tmp_sec_lens.append(tmp_line[3])
        sec_len_array.append(min(tmp_sec_lens))

###############################################################
#now to print out the data in a nice dataframe
#block  gene    sec dS  dN  omega   sec_len
print("block\tgene\tsec\tsec_two\tdS\tdN\tomega\tsec_len")
for s in range(0,len(block_array)):
    print(block_array[s],gene_array[s],sec_array[s],sec_array2[s],dS_array[s],dN_array[s],omega_array[s],sec_len_array[s],sep='\t')

