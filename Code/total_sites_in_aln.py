#!/usr/local/bin/python3
#
#program will go through *_coding_start_end.txt files and grab the gene
#start, end, and name info.
#will also go through the *.mafft alnignment file and grab the
#sequences and use the ref seq to split up the aln into multiple aln
#files, one for each gene
#run as: python3 ~/CODE/get_genes_from_mafft_to_align.py *.mafft REF_TAXA_NUMBER > out
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
# FOR DEBUGGING ONLY pdb.set_trace()
ref_taxa = int(sys.argv[2])

#below will open the mafft file and get the aln info
align = AlignIO.read(sys.argv[1], "fasta")
aln_len = len(align[ref_taxa].seq)
print(aln_len)
