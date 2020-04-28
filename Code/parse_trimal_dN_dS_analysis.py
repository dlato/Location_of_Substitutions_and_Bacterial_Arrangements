#!/usr/local/bin/python3
#
#program will take the assessed alignment from Gblocks and split it
#into chunks with proper genomic position locations
#run as: python3 ~/CODE/parse_Gblocks.py *.txt *.galn REF_TAXA_NUMBER > out
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
from collections import defaultdict
# FOR DEBUGGING ONLY pdb.set_trace()


#define function that will create a seq of numbers from a to b
def createList(r1, r2): 
    #seq needs to be reversed
    if r1 > r2:
        tmp_list = list(range(r2, r1+1))
        tmp_list = tmp_list[::-1]
        return tmp_list
    else:
        return list(range(r1, r2+1)) 

#below will open the split .txt aln file and get the aln info
align = AlignIO.read(sys.argv[1], "fasta")
aln_len = len(align[1].seq)

#getting the block name so we can print to separate
#files later
tmp_var = sys.argv[1]
sub_var = re.compile("\.txt")
tmp_var = sub_var.sub("",tmp_var)
block_name = tmp_var
sys.stderr.write(block_name)
sys.stderr.write("\n")
#ref taxa number
ref_taxa = int(sys.argv[3])

### ALIGNMENT STARTS AT POS 0 IN THE ARRAY
### TRIMAL STARTS AT POS 0 IN THE ARRAY
#edit the taxa IDs to get the start and end pos of each block
block_starts = []
block_ends = []
taxa_names = []
rev_comp = []
for i in range(0,len(align)):
    tmp_id = align[i].id
    #store 'a' which is the ':' character so we can sub it out later
    a = re.compile("\:")
    tmp_id = a.sub(" ",tmp_id)
    #store 'b' which is the '-' character so we can sub it out later
    b = re.compile("\-")
    tmp_id = b.sub(" ",tmp_id)
    #split the ID so we separate the actual ID from the start and end pos
    split_id = re.split(" ",tmp_id)
    #append start and stop and taxa name to list so we can use it later
    taxa_names.append(split_id[0])
    block_starts.append(split_id[1])
    block_ends.append(split_id[2])
    #storing info about if each seq is a reverse complement
    if (len(split_id) > 3):
        rev_comp.append(1)
        tmp_rev_comp = split_id[3]
    else:
        rev_comp.append(0)

block_starts = list(map(int,block_starts))
block_ends = list(map(int,block_ends))
num_taxa = len(block_ends)

#now make a list of lists to keep track of the codon classification
# codon 1, 2, 3
codon_class = []
for n in range(0, len(block_starts)):
    #rev comp
    if rev_comp[n] == 1:
        repetitions=[]
        torep=[3,2,1]
        for i in range(aln_len):
            repetitions=repetitions+torep
    #not rev comp
    else:
        repetitions=[]
        torep=[1,2,3]
        for i in range(aln_len):
            repetitions=repetitions+torep
    codon_class.append(repetitions)

#below will open the trimal .trimalcols file and get the aln info
#this file contains just the column numbers that trimal kept or
#threw out
with open(sys.argv[2], 'r') as f:
    reader = csv.reader(f)
    T_aln = list(reader)
T_aln[0] = [ int(x) for x in T_aln[0] ]
#make a new list that mimics the Gblocks output with a . for a
#non-conserved column, and a # for a conserved colunn
G_aln = ["."] * aln_len
for s in T_aln[0]:
    G_aln[s] = '#'

sys.stderr.write("done reading in the files! You go girl!\n")

#going through the Gblocks aln and getting the good and bad section
#starts and stops
G_starts = []
G_ends = []
G_class = []
G_codon = []
#G_aln = G_aln[0]
py_aln_len = len(G_aln) -1 
for count in range(0,len(G_aln)):
    G_cur = G_aln[count]
    G_last = G_aln[count-1]
    codon = codon_class[ref_taxa]
    #we are at the beginning of the aln
    if count == 0:
        G_starts.append(count)
        G_class.append(G_cur)
    #we are not at the beginning of the aln
    else:
        if count == py_aln_len:
            count = count
            G_ends.append(count)
            break
        #Gblocks found a bad section
        if G_cur != G_last:
            G_ends.append(count)
            G_starts.append(count)
            G_class.append(G_cur)
#            #checking if we are starting at the beginning of a codon
#            if (codon[count] == 3 or codon[count] == 1):
#                G_ends.append(count)
#                G_starts.append(count)
#                G_class.append(G_cur)
#            #we are not at beg of codon
#            else:
#                tmp_count = count + 2
#                #put the chunk that does not start at beg of codon in
#                #its own section and re label as bad (".")
#                G_ends.append(count)
#                G_starts.append(count)
#                G_class.append(".")
#                G_ends.append(tmp_count)
#                G_starts.append(tmp_count)
#                G_class.append(G_cur)


#editing starts and stops so that they all start and end with codons
#aln len of 3
updated_G_start = G_starts
updated_G_end = G_ends
updated_G_class = G_class
for i in range(0, len(G_starts)):
#    print("-------")
#    print("i", i)
#    print("G_class", G_class[i])
    codon = codon_class[ref_taxa]
    start = G_starts[i]
    end = G_ends[i]
    codon_s = codon[start]
    codon_e = codon[end]
    last_end = end -1
    next_start = start +1
#    print("start", start, "next_start", next_start, "end", end, "last_end", last_end)
#    print("codon_s", codon_s, "codon_e", codon_e, "codon[last_end]", codon[last_end], "codon[next_start]", codon[next_start])
    #only change things if we are in a good part of the aln
    if (G_class[i] == "#"):
        #if we start with codon pos 2 forward (so next should be 3)
        if (codon_s == 2 and codon[next_start] == 3):
            new_start = start +2
            updated_G_start[i] = new_start
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we start with codon pos 2 reverse (so next should be 1)
        if (codon_s == 2 and codon[next_start] == 1):
            new_start = start +2
            updated_G_start[i] = new_start
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we start with codon pos 3 forward (so next should be 1)
        if (codon_s == 3 and codon[next_start] == 1):
            new_start = start +1
            updated_G_start[i] = new_start
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we start with codon pos 1 reverse (so next should be 3)
        if (codon_s == 1 and codon[next_start] == 3):
            new_start = start +1
            updated_G_start[i] = new_start
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we end with codon pos 2 forward (so last should be 1)
        if (codon_e == 2 and codon[last_end] == 1):
            new_end = end -2
            updated_G_end[i] = new_end
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we end with codon pos 2 reverse (so last should be 3)
        if (codon_e == 2 and codon[last_end] == 3):
            new_end = end -2
            updated_G_end[i] = new_end
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we end with codon pos 1 forward (so last should be 3)
        if (codon_e == 1 and codon[last_end] == 3):
            new_end = end -1
            updated_G_end[i] = new_end
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
        #if we end with codon pos 3 reverse (so last should be 1)
        if (codon_e == 3 and codon[last_end] == 1):
            new_end = end -1
            updated_G_end[i] = new_end
            tmp_end = updated_G_end[i] + 1
            gene_len = tmp_end - updated_G_start[i]
#            print("gene_len", gene_len)
#        print("updated_G_start", updated_G_start[i], "updated_G_end", updated_G_end[i])


#making sure everything is minimum 100bp in length
for i in range(0,len(G_class)):
    if G_ends[i] - G_starts[i] <100:
        G_class[i] = "."

sys.stderr.write("got Gblocks section starts and ends!\n")


#make list of lists with each genomic pos listed at each site in the
#aln for each taxa
genome_pos = []
for t in range(0,num_taxa):
    tmp_pos_list = []
    # if seq is reversed
    if rev_comp[t] == 1:
        start = block_ends[t]
        end = block_starts[t]
        tmp_list = createList(start,end)
        genome_pos.append(tmp_list)
    #seq is not reversed
    else:
        start = block_starts[t]
        end = block_ends[t]
        tmp_list = createList(start,end)
        genome_pos.append(tmp_list)

#print out the good sections of the alignment!
#opening files so we can write to them separatly
good_sec_file_name = block_name + "_good.txt"
bad_sec_file_name = block_name + "_bad.txt"
good_file = open(good_sec_file_name, "w+")
bad_file = open(bad_sec_file_name, "w+")
#actually printing it out
for i in range(0,len(updated_G_start)):
    diff = updated_G_end[i] - updated_G_start[i] + 1
    if diff <= 102:
        #section is too short and therefore is bad
        for j in range(0,num_taxa):
            geom_pos = genome_pos[j]
            #rev comp
            if rev_comp[j] == 1:
                #write header to bad file
                bad_file.write("\n>")
                bad_file.write(taxa_names[j])
                bad_file.write(":")
                bad_file.write(str(geom_pos[int(updated_G_end[i])]))
                bad_file.write("-")
                bad_file.write(str(geom_pos[int(updated_G_start[i])]))
                bad_file.write(":Reversed\n")
                #write seq to bad file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    bad_file.write(align[j,k])
            #not rev comp
            else:
                #write header to bad file
                bad_file.write("\n>")
                bad_file.write(taxa_names[j])
                bad_file.write(":")
                bad_file.write(str(geom_pos[int(updated_G_start[i])]))
                bad_file.write("-")
                bad_file.write(str(geom_pos[int(updated_G_end[i])]))
                bad_file.write("\n")
                #write seq to bad file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    bad_file.write(align[j,k])
        #print separator for each section
        good_file.write("\n=\n")
        bad_file.write("\n=\n")
        continue
    #if we have a good section
    if G_class[i] == "#":
        for j in range(0,num_taxa):
            geom_pos = genome_pos[j]
            #rev comp
            if rev_comp[j] == 1:
                #write header to good file
                good_file.write("\n>")
                good_file.write(taxa_names[j])
                good_file.write(":")
                good_file.write(str(geom_pos[int(updated_G_end[i])]))
                good_file.write("-")
                good_file.write(str(geom_pos[int(updated_G_start[i])]))
                good_file.write(":Reversed\n")
                #write seq to good file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    good_file.write(align[j,k])
            #not rev comp
            else:
                #write header to good file
                good_file.write("\n>")
                good_file.write(taxa_names[j])
                good_file.write(":")
                good_file.write(str(geom_pos[int(updated_G_start[i])]))
                good_file.write("-")
                good_file.write(str(geom_pos[int(updated_G_end[i])]))
                good_file.write("\n")
                #write seq to good file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    good_file.write(align[j,k])
    # we have a bad section
    else:
        for j in range(0,num_taxa):
            geom_pos = genome_pos[j]
            #rev comp
            if rev_comp[j] == 1:
                #write header to bad file
                bad_file.write("\n>")
                bad_file.write(taxa_names[j])
                bad_file.write(":")
                bad_file.write(str(geom_pos[int(updated_G_end[i])]))
                bad_file.write("-")
                bad_file.write(str(geom_pos[int(updated_G_start[i])]))
                bad_file.write(":Reversed\n")
                #write seq to bad file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    bad_file.write(align[j,k])
            #not rev comp
            else:
                #write header to bad file
                bad_file.write("\n>")
                bad_file.write(taxa_names[j])
                bad_file.write(":")
                bad_file.write(str(geom_pos[int(G_start[i])]))
                bad_file.write("-")
                bad_file.write(str(geom_pos[int(updated_G_end[i])]))
                bad_file.write("\n")
                #write seq to bad file
                for k in range(updated_G_start[i],updated_G_end[i]+1):
                    bad_file.write(align[j,k])
    #print separator for each section
    good_file.write("\n=\n")
    bad_file.write("\n=\n")

sys.stderr.write("finito!\n")
