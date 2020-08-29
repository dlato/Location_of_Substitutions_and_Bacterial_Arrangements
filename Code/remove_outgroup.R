#remove outliers and branchlengts from tree files
# order of args:
# 1) tree file (newick)
# 2) outgroup taxa ID

##################
# libraries
library("ape")


#load in file:
#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
outgroup <- as.character(args[2])
#mytree <- read.tree('strep_chrom_superseq_bootstrapped_tree_outtree')
mytree <- read.tree(file_name)
#mytree <- read.tree('just_topology')
#plot(mytree)
no_outgroup <- drop.tip(mytree, outgroup)
#plot(no_outgroup)
write.tree(no_outgroup, "outgroup_tree_topology")
