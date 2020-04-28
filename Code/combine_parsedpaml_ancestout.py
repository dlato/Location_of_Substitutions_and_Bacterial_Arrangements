#!/usr/local/bin/python3
#
#program will combin the substitution information from a parsed paml
#output with the output from the ancestral reconstruction of positions
####################################

import sys, re
import Bio
from Bio import SeqIO
import csv
import pandas as pd
import numpy as np

branch = []
from1 = []
to1 = []
site = []
from2 = []
to2 = []
prob = []
#reading in parsedpaml file
f = open(sys.argv[1], "r")
for line in f:
    mo = re.search('^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([ACTG-])\s+([ATCG-])\s+(\d+\.\d+)', line)
    if mo:
        branch.append(mo.group(1))
        from1.append(mo.group(2))
        to1.append(mo.group(3))
        site.append(int(mo.group(4)))
        from2.append(mo.group(5))
        to2.append(mo.group(6))
        prob.append(mo.group(7))

#reading in ancestor_out file
ancest_nod = []
with open(sys.argv[2], 'rU') as infile:
    ancest = pd.read_csv(infile, sep="\t", header = None)
    for i in range(0,len(list(ancest.columns))):
        ancest_nod.append(ancest[i].tolist())
del ancest_nod[-1]
aln_len = len(ancest_nod[0])

#read in the order that paml branches are associated with ancest_out
#nodes
#zero means that it is the root of the ancest file for which paml has
#no branch for, therefore not used
p_nodes = []
pnf = open(sys.argv[3], "r")
for line in pnf:
    tmp_nod = re.split("\t",line)
    p_nodes = tmp_nod
del p_nodes[-1]
num_nodes = len(p_nodes)

sys.stderr.write("finished reading files!\n")

#NOTE: THE FROM1 AND TO1 ARE NOT CORRECT IF THERE IS NO SUB AT THAT
#SITE! *SHOULDNT* matter for R analysis
##go through each site in the alignment 
sys.stdout.write("#branch\tfrom1\tto1\tsite\tchange\tfrom2\tto2\tprob\n")
for s in range(0,aln_len):
    #make s the same as site num in paml file
    s2 = s + 1
    if s2 in site: #the site has at least one sub
        loc_of_subs = [i for i, x in enumerate(site) if x == s2]
        node_done = [0] * num_nodes
        for sub in loc_of_subs:
            #location of branch in ancestor out file
            a_nod_loc = p_nodes.index(branch[sub])
            a_nod1 = ancest_nod[a_nod_loc]
            a_pos1 = a_nod1[s]
            sys.stdout.write(branch[sub])
            sys.stdout.write("\t")
            sys.stdout.write(from1[sub])
            sys.stdout.write("\t")
            sys.stdout.write(to1[sub])
            sys.stdout.write("\t")
            sys.stdout.write(str(a_pos1))
            sys.stdout.write("\t1\t")
            sys.stdout.write(from2[sub])
            sys.stdout.write("\t")
            sys.stdout.write(to2[sub])
            sys.stdout.write("\t")
            sys.stdout.write(prob[sub])
            sys.stdout.write("\n")
            node_done[a_nod_loc] = 1
        #go throught the remainder of the nodes for this site
        #and print out null info to show no changes at any of these sites
        for i in range(0,num_nodes):
            paml_pos = int(p_nodes[i])
            if (node_done[i] == 0):
                if (paml_pos != 0) :
                    a_nod = ancest_nod[i]
                    a_pos = a_nod[s]
                    sys.stdout.write(str(paml_pos))
                    sys.stdout.write("\t1\t2\t")
                    sys.stdout.write(str(a_pos))
                    sys.stdout.write("\t0\tNA\tNA\tNA\n") 
    else: #the site does not have at least one sub
        #so print out null info to show no changes at any site
        for n in range(0,num_nodes):
            paml_pos = int(p_nodes[n])
            if (paml_pos != 0) :
                a_nod = ancest_nod[n]
                a_pos = a_nod[s]
                sys.stdout.write(str(paml_pos))
                sys.stdout.write("\t1\t2\t")
                sys.stdout.write(str(a_pos))
                sys.stdout.write("\t0\tNA\tNA\tNA\n") 



