#!/usr/local/bin/python3
#
# Program takes an rst output file from a PAML run and
# parses the output to make it easier for spreadsheets/R
#
# run as "parse_rst_file.py rst_file" 
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
# FOR DEBUGGING ONLY pdb.set_trace()
#
# file = open(sys.argv[1], 'r')
# using 'file' masks the builtin object 'file'
# using 'with' as below creates a context manager and closes automatically 
with open(sys.argv[1],mode="r",encoding='utf-8') as infile:
    for line in infile:
        ma=re.search("^Ancestral reconstruction by BASEML",line,re.IGNORECASE)
        if ma:
            line = next(infile) # blank line before tree
            print("# Original treefile used by PAML")
#            pdb.set_trace()
            while(1):  # infinite loop
                line = next(infile)
                print("# " + line.rstrip())
                if(re.search(";$",line)): break
        ma=re.search("^tree with node labels for Rod Page's TreeView",line,re.IGNORECASE)
        if ma:
            print("#\n# Newick tree with internal nodes labelled by PAML")
            while(1):  # infinite loop
                line = next(infile)
                print("# " + line.rstrip());
                if(re.search(";$",line)): break
            print("#")
            print("#Branch From  To       Site From   To  Prob")
        ma=re.search("^Branch (\d+):\s+(\d+)..(\d+)",line,re.IGNORECASE)
        if ma:
            edgeNumber = int(ma.group(1))
            branchFrom = int(ma.group(2))
            branchTo   = int(ma.group(3))
            line = next(infile)  # should be a blank line before changes
            while(1):  # infinite loop
                line = next(infile)
                if(re.search("^\s*$",line)): break # an empty line
                mb=re.search("^\s+(\d+) (\S+) (\S+) -> (\S+)",line)
                if mb:
                    site=int(mb.group(1))
                    char1=str(mb.group(2))
                    char2=str(mb.group(4))
                    prob =str(mb.group(3))
#                   print("{0:4d} {1:4d} {2:4d} {3:10d} {4:4s} {5:4s} {6:6s}".format(edgeNumber,branchFrom,branchTo,site,char1,char2,prob))  # This works ... new syntax
                    print("%4d %4d %4d %10d %4s %4s %6s" % (edgeNumber,branchFrom,branchTo,site,char1, char2,prob)) # This works ... old syntax
