#!/usr/local/bin/python3
#
#program will go through *_coding_start_end.txt files and grab the gene
#start, end, and name info.
#will also go through the *.mafft alnignment file and grab the
#sequences and use the ref seq to split up the aln into multiple aln
#files, one for each gene
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
# FOR DEBUGGING ONLY pdb.set_trace()

# function to get unique values
def unique(list1):

    # intilize a null list
    unique_list = []
    unique_pos = []

    # traverse for all elements
    i = 0
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
        else:
            unique_pos.append(i)
        i += 1
    return unique_list

#get index of duplicates for each itm in a list
def duplicates(lst, item):
    dups_list = []
    return [i for i, x in enumerate(lst) if x == item]

def unique_pos(my_list):
    #initalize vec of zeros
    unique_pos = [0] * len(my_list)
    #loop through each itm in list
    for x in my_list:
        #use above function to get index of duplicates of each item
        dups = duplicates(my_list, x)
        unique_pos[dups[-1]] = 1
    
    return unique_pos


##below will open the coding sequence start and end stuff
ref_taxa = int(sys.argv[2])
print("ref_tax num", ref_taxa)
print("ref_tax numclass", type(ref_taxa))
max_genome_len = int(sys.argv[3])
print("max_genome_len num", max_genome_len)
print("max_genome_len num class", type(max_genome_len))



#below will open all the coding seq start and end files and keep track
#of which taxa they come from (hopefully)
start_array = []
end_array = []
name_array = []
bac_name_array = []
#keeping track of if the gene is a complement or not
complement_array = []
for filepath in glob.glob(os.path.join('', '*_coding_start_end.txt')):
    with open(filepath) as f:
        bac_name = os.path.basename(filepath)
        #store 'a' which is the '_' character so we can sub it out later
        a = re.compile("\_coding\_start\_end")
        bac_name = a.sub(" ",bac_name)
        #split the ID so we separate the actual ID from the rest of
        #the file name
        split_id = re.split(" ",bac_name)
        #append start and stop and taxa name to list so we can use it later
        bac_name_array.append(split_id[0])
        data_gbk = pd.read_csv(f, sep="\t", header = None)
        data_gbk.columns = ["Start", "End", "Midpoint","Name", "Complement"]
        start_array.append(data_gbk[['Start']])
        end_array.append(data_gbk[['End']])
        name_array.append(data_gbk[['Name']])
        complement_array.append(data_gbk[['Complement']])

print(bac_name_array)

#below will open the mafft file and get the aln info
align = AlignIO.read(sys.argv[1], "fasta")
aln_len = len(align[ref_taxa].seq)
ref_tax_test = align[ref_taxa].id
print("REF TAX ID",ref_tax_test)

#getting the block name so we can print to 2 separate cod and non-cod
#files later
tmp_var = sys.argv[1]
sub_var = re.compile("\.mafft")
tmp_var = sub_var.sub("",tmp_var)
block_name = tmp_var
sys.stderr.write(block_name)

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
print("taxa_names")
print(taxa_names)

sys.stderr.write("done reading in the files!\n")

print("aln_len",aln_len)
#making a list that will have the order of the taxa from the alignment
#file in the order they are in the bac_names_array (bac order in start and end
#files being read in) to be used later
indices = []
for value in taxa_names:
    indices.append(bac_name_array.index(value))

#go through each col in aln and keep track of if there are gaps or not
#IN ALL GENOMES by making an array that is 0 if its is a gap, and the
#genome pos if it is not a gap
gaps_list = [] #list of arrays with each taxa's gap info
genome_pos_list = [] #list of arrays with each taxa's genome pos info
for t in range(0,num_taxa):
    gaps_array = np.zeros((aln_len), dtype=int)
    genome_pos_array = np.zeros((aln_len), dtype=int)
    rev_c = rev_comp[t]
    #if the seq is a reverse complement
    if (rev_c == 1):
        #actual genome pos of the ref seq in this block
        #indexing at zero
        #using end pos bc the seq is a rev comp
        genome_pos = block_ends[t] + 1
        print("genome_pos",genome_pos)
        for i in range(0,aln_len):
            if (i == 0):#if we are at the beginning of the aln
                if align[t,i] == "-":
                    gaps_array[i] = 0
                    genome_pos_array[i] = genome_pos
                else:
                    #working backwards bc it is a rev comp seq
                    genome_pos -= 1
                    gaps_array[i] = genome_pos
                    genome_pos_array[i] = genome_pos
            else: #if we are not a the beginning of the aln   
                last_pos = genome_pos +1
                if align[t,i] == "-":
                    gaps_array[i] = 0
                    genome_pos_array[i] = genome_pos
                else:
                    #working backwards bc it is a rev comp seq
                    genome_pos -= 1
                    gaps_array[i] = genome_pos
                    genome_pos_array[i] = genome_pos
    #if the seq is not a reverse complement
    else:
        #actual genome pos of the ref seq in this block
        genome_pos = block_starts[t] -1
        for i in range(0,aln_len):
            if (i == 0):#if we are at the beginning of the aln
                if align[t,i] == "-":
                    gaps_array[i] = 0
                    genome_pos_array[i] = genome_pos
                else:
                    genome_pos += 1
                    gaps_array[i] = genome_pos
                    genome_pos_array[i] = genome_pos
            else: #if we are not a the beginning of the aln   
                last_pos = genome_pos -1
                if align[t,i] == "-":
                    gaps_array[i] = 0
                    genome_pos_array[i] = genome_pos
                else:
                    genome_pos += 1
                    gaps_array[i] = genome_pos
                    genome_pos_array[i] = genome_pos
    gaps_list.append(gaps_array)
    genome_pos_list.append(genome_pos_array)

print("genome_pos_list")
print(genome_pos_list)
print(len(gaps_list))
print(gaps_list)

#tot_gaps is an array that has a 1 if there is a gap and a 0 if there
#is not a gap ANYWHERE in the anl. so minimum one gap in any col
tot_gaps = np.zeros((aln_len), dtype=int)
for i in range(0,aln_len):
    for j in range(0,num_taxa):
        if align[j,i] == "-":
            tot_gaps[i] = 1

#up to this point takes maybe 30 seconds to run
#the part below takes about 12 min to run

#dealing with the fact that the complement array is out of order
ord_rev = ["0"] * num_taxa
i = 0
for t in indices:
    ord_rev[t] = rev_comp[i]
    i += 1

sys.stderr.write("done dealing with gaps! yaaas queen!\n")

#set up same array as before that says coding and non coding except it
#now puts in the gene names instead of a 1
#this will NOT put the cod_list and gng_list into the same order
#as the alignment by using the indices array...for some weird reason
#specify the absoloute max of the genome (aka biggest genome in data)
cod_list = []
gng_list = []
gene_name_list = []
gene_comp_list = []
nex_gene_codon_pos = 0
for t in indices:
    print("t", t, bac_name_array[t])
    j = 0
    last_gene = 0
    cod_array = ["0"] * max_genome_len
    #making an array to keep track of gene names to be used later
    gene_name_array = ["0"] * max_genome_len
    gene_comp_array = ["0"] * max_genome_len
    gene_cod_array = ["0"] * max_genome_len
    #make arrays to keep track of genes within genes that need to be
    #accounted for later. this array will tell us what pos in the actual
    #gene start and end arrays need to be accounted for later
    gng = []
    s_array = np.asarray(start_array[t])
    e_array = np.asarray(end_array[t])
    n_array = np.asarray(name_array[t])
    c_array = np.asarray(complement_array[t])
    length = len(e_array)
    end_len = length
    codon_pos = 0
    alt_codon_pos = 0
    overlap_end = 0
    rev_c = ord_rev[t]
    #for each site in the genome
    for i in range(0,max_genome_len):
        #below 4 lines take care of resetting the codon_pos after it
        #hits 3 within a gene. so after it is done with one codon
        if (codon_pos == 3):
            codon_pos = 0
        if (nex_gene_codon_pos == 3):
            nex_gene_codon_pos = 0
        if (j >= length):
            break
        #dealing with j being bigger than the # of genes/rows in start/end file
        #site within coding region!
        #making sure that the start and end pos are on a 0 array start
        #system
        start = s_array[j] -1
        end = e_array[j]-1
        name = n_array[j]
        comple = c_array[j]
        #if the next gene overlaps this gene
        if (s_array[j+1] >= s_array[j] and s_array[j+1] <= e_array[j]):
            overlap_diff = e_array[j+1] - e_array[j]
            #if difference in overlap between these 2 genes is negative, then we have only overlap (no gene in gene)
            if (overlap_diff < 0):
                overlap_end = e_array[j+1] - 1
            #if difference in overlap between these 2 genes is positive, then we have a gene in gene situation
            else:
                overlap_end = e_array[j] - 1
        if (j == 0):
            if (i >= start and i <= end):
                codon_pos += 1
                if (comple == 1):
                    #for complement genes switching 1->3 and 3->1
                    gene_name_array[i] = n_array[j] 
                    gene_comp_array[i] = c_array[j] 
                    if (codon_pos == 1):
                        cod_array[i] = 3
                    else:
                        if (codon_pos == 3):
                            cod_array[i] = 1
                        else:
                            cod_array[i] = codon_pos
                else:
                    #if the gene is not a complement then codon_pos is
                    #unchanged
                    cod_array[i] = codon_pos
                    gene_name_array[i] = n_array[j] 
                    gene_comp_array[i] = c_array[j] 
                if i == end:
                    j += 1
                    codon_pos = 0
                    if j == end_len:
                        break
        else:
            nex = j + 1
            if (nex == len(e_array)):
                if (i >= start and i <= end):
                    codon_pos += 1
#                    if (comple == 1 or rev_c == 1):
                    if (comple == 1):
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                        #for complement genes switching 1->3 and 3->1
                        if (codon_pos == 1):
                            cod_array[i] = 3
                        else:
                            if (codon_pos == 3):
                                cod_array[i] = 1
                            else:
                                cod_array[i] = codon_pos
                    else:
                        #if the gene is not a complement then codon_pos is
                        #unchanged
                        cod_array[i] = codon_pos
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                    if i == end:
                        j += 1
                        codon_pos = 0
                        if j == end_len:
                            break
            if (e_array[j+1] <= e_array[j]):#if the next gene is in this gene
                if (e_array[j+1] == 0):
                    last_gene += 1
                if (i >= start and i <= end):
                    codon_pos += 1
                    if (comple == 1):
#                    if (comple == 1 or rev_c == 1):
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                        #for complement genes switching 1->3 and 3->1
                        if (codon_pos == 1):
                            cod_array[i] = 3
                        else:
                            if (codon_pos == 3):
                                cod_array[i] = 1
                            else:
                                cod_array[i] = codon_pos
                    else:
                        #if the gene is not a complement then codon_pos is
                        #unchanged
                        cod_array[i] = codon_pos
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                    if i == end:
                        #to skip the next gene bc its in the range of the
                        #previous gene
                        j += 2
                        gng.append(j-1)
                        if j == end_len:
                            break
                        if i == overlap_end: #if we are at the end of the overlap to reset alt_codon_pos for next overlap
                            alt_codon_pos = 4
                            gene_name_array[i] = 4 
                            gene_comp_array[i] = 4 
                            cod_array[i] = alt_codon_pos
                            nex_gene_codon_pos += 1
                            #making sure the codon pos is now following the actual codon pos of the next gene
                            codon_pos = nex_gene_codon_pos
                            nex_gene_codon_pos = 0
                            alt_codon_pos = 0
            else:#else if the next gene is not in this gene
                if (i >= start and i <= end):
                    codon_pos += 1
                    if (comple == 1):
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                        #for complement genes switching 1->3 and 3->1
                        if (codon_pos == 1):
                            cod_array[i] = 3
                        else:
                            if (codon_pos == 3):
                                cod_array[i] = 1
                            else:
                                cod_array[i] = codon_pos
                    else:
                        #if the gene is not a complement then codon_pos is
                        #unchanged
                        gene_name_array[i] = n_array[j] 
                        gene_comp_array[i] = c_array[j] 
                        cod_array[i] = codon_pos
                    if i == end:
                        if (i >= s_array[j+1]-1 and i <= e_array[j+1]-1):
                            #keeping track of what the codon pos should be
                            nex_gene_codon_pos += 1
                            codon_pos = nex_gene_codon_pos
                            alt_codon_pos = 4
                            cod_array[i] = alt_codon_pos
                            gene_name_array[i] = 4 
                            gene_comp_array[i] = 4 
                        #making sure the codon pos is now following the actual codon pos of the next gene
                        codon_pos = nex_gene_codon_pos
                        alt_codon_pos = 0
                        nex_gene_codon_pos = 0
                        j += 1
                        if j >= end_len:
                            last_gene += 1
                        if last_gene == 2:
                            break
                #if there is overlap at this site with the next gene
                if (i >= s_array[j+1]-1 and i <= e_array[j+1]-1):
                    #keeping track of what the codon pos should be
                    nex_gene_codon_pos += 1
                    alt_codon_pos = 4
                    cod_array[i] = alt_codon_pos
                    gene_name_array[i] = 4 
                    gene_comp_array[i] = 4 
                    if i == overlap_end: #if we are at the end of the overlap to reset alt_codon_pos for next overlap
                        alt_codon_pos = 4
                        cod_array[i] = alt_codon_pos
                        gene_name_array[i] = 4 
                        gene_comp_array[i] = 4 
                        #making sure the codon pos is now following the actual codon pos of the next gene
                        codon_pos = nex_gene_codon_pos
                        alt_codon_pos = 0
                        nex_gene_codon_pos = 0
    cod_list.append(cod_array)
    gene_name_list.append(gene_name_array)
    gene_comp_list.append(gene_comp_array)
    gng_list.append(gng)


#now to go through the gng array (which has a list of all genes that
#completely fall within other genes) and change the value of those
#positions in the cod_list to 4 (denoting overlap)
gng_lens = []
for item in gng_list:
    length_gng = len(item)
    gng_lens.append(length_gng)
    print("gng_len",length_gng)
countt = 0
if gng_list:
    for last in gng_list:
        if (gng_lens[countt] != 0):
            del last[-1]
        countt += 1
#dealing with the fact that the start and end arrays are out of order
ord_s_array = []
ord_e_array = []
for t in indices:
    s_array = np.asarray(start_array[t])
    e_array = np.asarray(end_array[t])
    ord_s_array.append(s_array)
    ord_e_array.append(e_array) 

#now dealing with gene in gene lists and making sure that they are all
#classified as codon pos 4
updated_cod_list = []
updated_gene_name_list = []
updated_gene_comp_list = []
for i in range(0,num_taxa):
    gng_array = gng_list[i]
    cod_array = cod_list[i]
    gene_name_array = gene_name_list[i]
    gene_comp_array = gene_comp_list[i]
    s_array = np.asarray(ord_s_array[i])
    e_array = np.asarray(ord_e_array[i])
    if (gng_lens[i] != 0):
        for nums in gng_array:
            start = np.asscalar(s_array[nums]) - 1
            end = np.asscalar(e_array[nums]) - 1
            for j in range(start, end):
                cod_array[j] = 4
                gene_name_array[j] = 4
                gene_comp_array[j] = 4
    else:
        cod_array = cod_list[i]
        gene_name_array = gene_name_list[i]
        gene_comp_array = gene_comp_list[i]
    updated_cod_list.append(cod_array)
    updated_gene_name_list.append(gene_name_array)
    updated_gene_comp_list.append(gene_comp_array)

for item in updated_cod_list:
    print("updated cod len", len(item))
for item in updated_gene_name_list:
    print("updated gene_name len", len(item))


sys.stderr.write("get those genome classifications!\n")

##now go through the gaps array to see when it is within a gene or not
#if it is not in a gene then it is a 0, otherwise it has the gene name
#this is the part that is actually accounting for the gaps in the ref
#genome! so gene_name already has gaps accounted for!
#does the same thing for the gene name array that is used later to get
#the aln start and stop pos for each gene. If there is a gap, then it
#will have the previous site's classification/gene name (only for the
#gene_name list)
codon_pos_list = []
aln_gene_name_list = []
aln_gene_comp_list = []
for t in range(0,num_taxa):
#for t in range(5,num_taxa):
    gene_name = ["0"] * aln_len
    aln_gene_name = ["0"] * aln_len
    aln_gene_comp = ["0"] * aln_len
    k = 0
    gaps_array = gaps_list[t]
    cod_array = updated_cod_list[t]
    gn_array = updated_gene_name_list[t]
    gc_array = updated_gene_comp_list[t]
    pos_array = genome_pos_list[t]
    #getting the reversed complement info
    print("t", t, "taxa", taxa_names[t])
    print("gaps_array", len(gaps_array))
    print("cod_array", len(cod_array))
#    #below will deal with a weird error that I think happens because the array is too big?
    for x in np.nditer(gaps_array):
        print("-------------------")
        print("x",x)
        print("k",k)
        tmp_val = x -1
        print("tmp_val",tmp_val)
        print("aln",align[t,k])
        if x != 0:#there is no gap at this site
            #to make sure that the indexing is right we minus 1
            tmp_x = x -1
            print("tmp_x",tmp_x)
            #if the site is within a gene
            if (cod_array[tmp_x] != '0'): #we are in a coding section of that genome
                gene_name[k] = cod_array[tmp_x]
                aln_gene_name[k] = gn_array[tmp_x]
                aln_gene_comp[k] = gc_array[tmp_x]
            else: #we are in a non-coding section of that genome
                gene_name[k] = cod_array[tmp_x]
                aln_gene_name[k] = gn_array[tmp_x]
                aln_gene_comp[k] = gc_array[tmp_x]
        else: #there is a gap at that site
            print("gap", x)
            gene_name[k] = "-"
            print("craps out here 1")
            print("pos_array[k]",pos_array[k])
            print("len(gn_array", len(gn_array))
            print("gn_array[pos_array[k]]",gn_array[pos_array[k]])
            print("aln_gene_name[k]",aln_gene_name[k])
            aln_gene_name[k] = gn_array[pos_array[k]]
            print("craps out here 2")
            aln_gene_comp[k] = gc_array[tmp_val]
        print("RECORDED codon_pos", gene_name[k])
        print("RECORDED aln_gene_name", aln_gene_name[k])
        k += 1
    codon_pos_list.append(gene_name)
    aln_gene_name_list.append(aln_gene_name)
    aln_gene_comp_list.append(aln_gene_comp)

#checking if each col in the alignment actually has all the same codon
#pos
#and storing the info for the overall class of each site in the aln: 
#x = not the same class
#1,2,3 = all taxa coding class
#0 = all taxa non-coding class
#4 = all taxa misc
#- = at least one gap in column of alignment
overall_class = []
for i in range(0,aln_len):
    print("--------------------")
    print("aln site:", i)
    tot_gap = tot_gaps[i]
    column = []
    print("aln nucs", align[:,i])
    #creating function to check if all items are the same in a list
    def all_same(items):
        return all(x == items[0] for x in items)
    #storing info about each nuc in that col of the alignment
    for t in range(0,num_taxa):
        print(taxa_names[t])
        tax_aln = codon_pos_list[t]
        print("tax_aln[i]",tax_aln[i])
        column.append(tax_aln[i])
    print("col",column, "tot_gap", tot_gap)
    #if there is a gap
    if (tot_gap == 1):
        overall_class.append("-")
    else:
        #checking if all sites @ the col are the same
        if all_same(column):
            overall_class.append(column[0])
        else:
            overall_class.append("x")



#a list of arrays that hold section info for each taxa, so if a gene
#is separated by gaps, it keeps that info
#columns we are tossing (-, not the same class in that col, 4) = x
#coding = gene name (like b0010)
#non-coding = 0
section_name_list = []
section_gene_comp_list = []
#combining the info from the overall_class list and the gene name
# to be used to get the aln gene starts and stops
for t in range(0,num_taxa):
    section_name = ["0"] * aln_len
    section_gene_comp = ["0"] * aln_len
    gn_array = aln_gene_name_list[t]
    gc_array = aln_gene_comp_list[t]
    #getting the reversed complement info
#    #below will deal with a weird error that I think happens because the array is too big?
    for x in range(0, aln_len):
        print("-------------------")
        print("x",x)
        o_class = overall_class[x]
        print("o_class", o_class)
        #if we are in coding section
        if (o_class == 1):
            #then our section_name has the same name as our gene
            section_name[x] = gn_array[x]
            section_gene_comp[x] = gc_array[x]
        if (o_class == 2):
            #then our section_name has the same name as our gene
            section_name[x] = gn_array[x]
            section_gene_comp[x] = gc_array[x]
        if (o_class == 3):
            #then our section_name has the same name as our gene
            section_name[x] = gn_array[x]
            section_gene_comp[x] = gc_array[x]
        #if we are at a site we will toss
        if (o_class == "-"):
            section_name[x] = "x"
            section_gene_comp[x] = gc_array[x]
        #if we are at a site we will toss
        if (o_class == 4):
            section_name[x] = "x"
            section_gene_comp[x] = gc_array[x]
        #if we are at a site we will toss
        if (o_class == "x"):
            section_name[x] = "x"
            section_gene_comp[x] = gc_array[x]
        #if we are in non-coding
        if (o_class == 0):
            section_name[x] = o_class
            section_gene_comp[x] = gc_array[x]
    section_name_list.append(section_name)
    section_gene_comp_list.append(section_gene_comp)


sys.stderr.write("combining class and gaps is done\n")

#getting the aln start and end positions for each gene for each genome
aln_start_list = []
aln_end_list = []
gene_class_list = []
gene_len_list = []
gene_comp_list = []
for t in range(0, num_taxa):
    print("taxa", taxa_names[t])
    j = 0
    gene_name = section_name_list[t]
    gene_comp_arr = section_gene_comp_list[t]
    aln_start = []
    aln_end = []
    gene_class = []
    gene_comp = []
    gene_len_arr = []
    print("len of gene_name",len(gene_name))
    for site in gene_name:
        print("-----------")
        print("site gene name:",site,"j:",j)
        if j == 0:
            next_site = j + 1
            start_v = 0
            if (site != '0'):
                if (site != gene_name[next_site]):
                    print("aln start",j)
                    print("aln end",j)
                    if (len(aln_start) != 0):
                        aln_end.append(j)
                    if (start_v == 0): #first site of aln
                        gene_len = j - start_v
                        gene_len_arr.append(gene_len)
                        print("gene_len_arr val", gene_len)
                        if (len(aln_start) == 0):
                            aln_start.append(0)
                            aln_end.append(j)
                            gene_class.append(str(site))
                            gene_comp.append(int(gene_comp_arr[j]))
                    else:
                        tmp_j = j +1
                        gene_len = tmp_j - start_v
                        gene_len_arr.append(gene_len)
                        print("gene_len_arr val", gene_len)
                        if (len(aln_start) == 0):
                            aln_start.append(0)
                            aln_end.append(j)
                            gene_class.append(str(site))
                            gene_comp.append(int(gene_comp_arr[j]))
                else:
                    print("aln start a",j, "gene_class", site)
                    gene_class.append(str(site))
                    gene_comp.append(int(gene_comp_arr[j]))
                    start_v = j
                    aln_start.append(j)
        else:
            next_site = j + 1
            last_site = j - 1
            print("j", j,"next", next_site, "last", last_site)
            if (next_site == aln_len):
                if (len(aln_start) != 0):
                    aln_end.append(j)
                if (start_v == 0): #first site of aln
                    gene_len = j - start_v
                    gene_len_arr.append(gene_len)
                    print("gene_len_arr val", gene_len)
                    if (len(aln_start) == 0):
                        aln_start.append(0)
                        aln_end.append(j)
                        gene_class.append(str(site))
                        gene_comp.append(int(gene_comp_arr[j]))
                else:
                    tmp_j = j +1
                    gene_len = tmp_j - start_v
                    gene_len_arr.append(gene_len)
                    print("gene_len_arr val", gene_len)
                    if (len(aln_start) == 0):
                        aln_start.append(0)
                        aln_end.append(j)
                        gene_class.append(str(site))
                        gene_comp.append(int(gene_comp_arr[j]))
                print("aln end c", j)
                break
            if (site != gene_name[last_site]):
                print("aln start b",j, "gene_class", site)
                gene_class.append(str(site))
                gene_comp.append(int(gene_comp_arr[j]))
                start_v = j
                aln_start.append(j)
            if (site != gene_name[next_site]):
                print("aln end b",j)
                if (len(aln_start) != 0):
                    aln_end.append(j)
                if (start_v == 0): #first site of aln
                    gene_len = j - start_v
                    gene_len_arr.append(gene_len)
                    print("gene_len_arr val", gene_len)
                    if (len(aln_start) == 0):
                        aln_start.append(0)
                        aln_end.append(j)
                        gene_class.append(str(site))
                        gene_comp.append(int(gene_comp_arr[j]))
                else:
                    tmp_j = j +1
                    gene_len = tmp_j - start_v
                    gene_len_arr.append(gene_len)
                    print("gene_len_arr val", gene_len)
                    if (len(aln_start) == 0):
                        aln_start.append(0)
                        aln_end.append(j)
                        gene_class.append(str(site))
                        gene_comp.append(int(gene_comp_arr[j]))
        j += 1
    aln_start_list.append(aln_start)
    aln_end_list.append(aln_end)
    gene_class_list.append(gene_class)
    gene_comp_list.append(gene_comp)
    gene_len_list.append(gene_len_arr)


#editing the aln starts and stops to make sure the aln len is a
#multiple of 3
updated_aln_start_list = aln_start_list
updated_aln_end_list = aln_end_list
updated_gene_len_list = gene_len_list
#do this for each taxa
for t in range(0,num_taxa):
    print("taxa", taxa_names[t])
    aln_start_arr = aln_start_list[t]
    aln_end_arr = aln_end_list[t]
    updated_start = updated_aln_start_list[t]
    updated_end = updated_aln_end_list[t]
    gene_class = gene_class_list[t]
    updated_gene_len = updated_gene_len_list[t]
    #go through each start and stop
    for i in range(0,len(aln_start_arr)):
        print("-----")
        print("i", i)
        print("gene_class", gene_class[i])
        print("gene_len", updated_gene_len[i])
        start = aln_start_arr[i]
        end = aln_end_arr[i]
        o_class_start = overall_class[start]
        o_class_end = overall_class[end]
        last_end = end -1
        next_start = start +1
        print("start", start, "next_start", next_start,"end", end,"last_end", last_end)
        print("oclass[start]", o_class_start, "oclass[end]", o_class_end, "oclass[next_start", overall_class[next_start], "oclass[last_end]", overall_class[last_end])
        #if we start with codon pos 2 forward (so next should be 3)
        if (o_class_start == 2 and overall_class[next_start] == 3):
            new_start = start +2
            updated_start[i] = new_start
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we start with codon pos 2 reverse (so next should be 1)
        if (o_class_start == 2 and overall_class[next_start] == 1):
            new_start = start +2
            updated_start[i] = new_start
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we start with codon pos 3 forward (so next should be 1)
        if (o_class_start == 3 and overall_class[next_start] == 1):
            new_start = start +1
            updated_start[i] = new_start
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we start with codon pos 1 reverse (so next should be 3)
        if (o_class_start == 1 and overall_class[next_start] == 3):
            new_start = start +1
            updated_start[i] = new_start
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we end with codon pos 2 forward (so last should be 1)
        if (o_class_end == 2 and overall_class[last_end] == 1):
            new_end = end -2
            updated_end[i] = new_end
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we end with codon pos 2 reverse (so last should be 3)
        if (o_class_end == 2 and overall_class[last_end] == 3):
            new_end = end -2
            updated_end[i] = new_end
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we end with codon pos 1 forward (so last should be 3)
        if (o_class_end == 1 and overall_class[last_end] == 3):
            new_end = end -1
            updated_end[i] = new_end
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        #if we end with codon pos 3 reverse (so last should be 1)
        if (o_class_end == 3 and overall_class[last_end] == 1):
            new_end = end -1
            updated_end[i] = new_end
            tmp_end = updated_end[i] + 1
            gene_len = tmp_end - updated_start[i]
            updated_gene_len[i] = gene_len
        print("updated_start", updated_start[i], "updated_end", updated_end[i])
        print("updated_gene_len", updated_gene_len[i])
    updated_aln_start_list[t] = updated_start
    updated_aln_end_list[t] = updated_end
    updated_gene_len_list[t] = updated_gene_len


#now to print out each gene alignment separated by =. Only issue is
#that each time you print something it puts it as a new line or with a
#space. So I will have to use a pearl script to remove those. So I put
#a < at the end of each taxa title so it can be easily switched into a
#new line.
#Right now the start and end pos of each gene is the whole block
#start and end for each taxa.Not the actual location of that gene in
#the genome. May have to fix this eventually if I do
#any ancestral reconstruction
#for each aln start (should be the same as number of genes in this
#aln)
#creating unique gene_class_list
gene_class_array = gene_class_list[ref_taxa]
gene_class_unique = unique(gene_class_array)
print("unique gene_class len", len(gene_class_unique))
print(gene_class_unique)
if "0" in gene_class_unique:
    gene_class_unique.remove("0")
if "x" in gene_class_unique:
    gene_class_unique.remove("x")
gene_class_unique_pos = unique_pos(gene_class_array)
#opening to 2 files so we can print to them separatly
cod_file_name = block_name + "_coding.txt"
non_cod_file_name = block_name + "_noncoding.txt"
cod_file = open(cod_file_name,"w+")
non_cod_file = open(non_cod_file_name,"w+")
#using the ref taxa as the aln start and ends
ref_aln_start = updated_aln_start_list[ref_taxa]
ref_aln_end = updated_aln_end_list[ref_taxa]
current_gene = []
gene_count = -1
q = -1
for i in range(0,len(ref_aln_start)):
    print("-----------")
    non_div = 0
    gene_gapped = 0
    gene_len_arr = updated_gene_len_list[ref_taxa]
    gene_len = gene_len_arr[i]
    gene_class_arr = gene_class_list[ref_taxa]
    gene_class = gene_class_arr[i]
    print("gene_class", gene_class)
    print("gene_len", gene_len)
    #looping through each of the taxa
    for j in range(0,num_taxa):
        gene_aln_len = 0
        genome_pos_arr = genome_pos_list[j] 
        gene_comp_arr = gene_comp_list[j]
        gene_comp = gene_comp_arr[i]
        #if we are NOT a gap, 4, or col where not everyone is the same
        if (gene_class != "x"):
            #checking if the gene is a mult of 3
            #if we are in non-coding region
            if (gene_class == "0"):
                #looping through each aln column between the aln gene start
                #and aln gene end
                #seq is reverse complement
                if (rev_comp[j] == 1):
                    non_cod_file.write("\n>")
                    non_cod_file.write(taxa_names[j])
                    non_cod_file.write(" ")
                    star = str(genome_pos_arr[int(ref_aln_start[i])])
                    endd = str(genome_pos_arr[int(ref_aln_end[i])])
                    non_cod_file.write(endd)
                    non_cod_file.write("-")
                    non_cod_file.write(star)
                    non_cod_file.write(" Reversed:")
                    non_cod_file.write("\n")
                    print(">TEST",taxa_names[j],":",endd,"-",star," Reversed:")
                else:
                    if (gene_comp == 1):
                        non_cod_file.write("\n>")
                        non_cod_file.write(taxa_names[j])
                        non_cod_file.write(" ")
                        star = str(genome_pos_arr[int(ref_aln_start[i])])
                        endd = str(genome_pos_arr[int(ref_aln_end[i])])
                        non_cod_file.write(star)
                        non_cod_file.write("-")
                        non_cod_file.write(endd)
                        non_cod_file.write(" Reversed:")
                        non_cod_file.write("\n")
                        print(">TEST",taxa_names[j],":",star,"-",endd," Reversed:")
                    else:
                        non_cod_file.write("\n>")
                        non_cod_file.write(taxa_names[j])
                        non_cod_file.write(":")
                        star = str(genome_pos_arr[int(ref_aln_start[i])])
                        endd = str(genome_pos_arr[int(ref_aln_end[i])])
                        non_cod_file.write(star)
                        non_cod_file.write("-")
                        non_cod_file.write(endd)
                        non_cod_file.write("\n")
                        print(">TEST",taxa_names[j],":",star,"-",endd)
                for k in range(int(ref_aln_start[i]),int(ref_aln_end[i])+1):
                    non_cod_file.write(align[j,k])
            #or we are in a coding region
            else:
                if (gene_len <= 1):
                    print("gene len <= 1")
                    continue;
                #seq is reverse complement
                if (gene_len % 3 != 0):
                    print("a NOT a mult of 3")
                if (rev_comp[j] == 1):
                    cod_file.write("\n>")
                    cod_file.write(taxa_names[j])
                    cod_file.write(" ")
                    star = str(genome_pos_arr[int(ref_aln_start[i])])
                    endd = str(genome_pos_arr[int(ref_aln_end[i])])
                    cod_file.write(endd)
                    cod_file.write("-")
                    cod_file.write(star)
                    cod_file.write(" Reversed:")
                    cod_file.write("\n")
                    print(">TEST",taxa_names[j],":",endd,"-",star)
                else:
                    #gene itself is a complement!
                    if (gene_comp == 1):
                        cod_file.write("\n>")
                        cod_file.write(taxa_names[j])
                        cod_file.write(" ")
                        star = str(genome_pos_arr[int(ref_aln_start[i])])
                        endd = str(genome_pos_arr[int(ref_aln_end[i])])
                        cod_file.write(star)
                        cod_file.write("-")
                        cod_file.write(endd)
                        cod_file.write(" Reversed:")
                        cod_file.write("\n")
                        print(">TEST",taxa_names[j],":",star,"-",endd, " Reversed")
                    else:
                        cod_file.write("\n>")
                        cod_file.write(taxa_names[j])
                        cod_file.write(":")
                        star = str(genome_pos_arr[int(ref_aln_start[i])])
                        endd = str(genome_pos_arr[int(ref_aln_end[i])])
                        cod_file.write(star)
                        cod_file.write("-")
                        cod_file.write(endd)
                        cod_file.write("\n")
                        print(">TEST",taxa_names[j],":",star,"-",endd)
                #looping through each aln column between the aln gene start
                #and aln gene end
                current_gene.append(gene_class)
                gene_count += 1
                for k in range(int(ref_aln_start[i]),int(ref_aln_end[i])+1):
                    cod_file.write(align[j,k])
    #end of section
    cod_file.write("\n=\n")
    non_cod_file.write("\n=\n")
    #checking if the gene we are looking at is in our unique gene list
    if (gene_class in gene_class_unique):
        #places where we have the end of a gene
        if (gene_class_unique_pos[i] == 1):
            cod_file.write("\n&\n")
