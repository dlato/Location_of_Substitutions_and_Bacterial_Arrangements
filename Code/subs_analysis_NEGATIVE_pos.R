#will do all the R analysis for subs data
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME GENOME_POSITION_DATAFRAME_FILE LENGTH_OF_GENOME ORIC_POS TERMINUS_POS
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME LENGTH_OF_GENOME ORIC_POS TERMINUS_POS
# order of args:
# 1) expression file (csv)
# 2) bacteria name
# 3) replicon name 
# 4) max genome position
# 5) origin of replication location
# 6) terminus of replication location
# 7) length of chunks
# 8) out file bac name
# 8) cutoff for bp near and far from ori
library("RColorBrewer")
loadNamespace("ggplot2")
library("gplots")
library("ggplot2")
library(ggplot2)
library(plyr)
library(dplyr)
library(proto)
library("proto")
is.proto.proto <- is.proto
library(reshape)
library(testthat)
library(gridExtra)
library(grid)
library(lattice)

#set theme
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            #change size of facet header text
            theme(strip.text = element_text(size =10.49)) +
            theme(
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_line(colour = "grey90", size = 0.2),
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5),
                  panel.spacing = unit(0.25, "lines"),
                  # plot title
                  plot.title = element_text(hjust = 0.5, size = 18),
                  # axis info
                  axis.text.y = element_text(),
                  axis.title.y = element_text(),
                  axis.ticks.y = element_line(),
                  axis.text=element_text(size=18),
                  axis.title = element_text(size = 18),
                  axis.text.x= element_text(),
                  axis.title.x = element_text(),
                  axis.ticks.x = element_line(),
                  #plot margins
                  #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(size=18),
                  #legend
                  legend.title = element_blank(),
                  legend.position = c(0.9,0.9),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend.position="top")
)
)




#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
bac_name <- as.character(args[2])
replicon <- as.character(args[3])
out_file_name <- as.character(args[8])

#read in the data, which includes NA's and spots where there are no
#subs
options("scipen"=100, "digits"=5)
datafile <- file_name
merged_data <- read.table(datafile, sep = "\t")
colnames(merged_data) <- c("branch", "from1", "to1", "position", "change", "from2", "to2", "prob")
head(merged_data)
which(merged_data$from2 == "G")
which(merged_data$from2 == "N")
which(merged_data$to2 == "N")

#removing any substitutions that involve ambiguous nucleotides
ambig_nucs_dat <- merged_data[grepl("[A,C,T,G]", merged_data[["from2"]]) | grepl("[A,C,T,G]", merged_data[["from2"]]),]
ambig_nucs_dat
which(ambig_nucs_dat$from2 == "N")

mytitle <- substitute(italic(bac_name)~replicon, list(bac_name=bac_name, replicon=replicon))
options("scipen"=100, "digits"=10)
################################################################################
#ORIGIN SCALING AND BIDIRECTIONALITY
################################################################################
#first scaling things to the origin (if necessary)
max_pos <- as.numeric(args[4])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[5])
print("oriC")
oriC_pos
terminus <- as.numeric(args[6])
print("ter")
terminus
new_pos <- merged_data$position
tmp_pos <- merged_data$position
print("MIN POS")
min(merged_data$position)
head(merged_data)
 
if (bac_name == "E.coli" | replicon == "pSymA") {
  to_shift_ter <- max_pos - oriC_pos
  shifted_ter <-terminus + to_shift_ter
  terminus <- shifted_ter
}

if (replicon == "pSymB") {
  shifted_ter <- terminus - oriC_pos
  terminus <- shifted_ter
}

if (bac_name == "E.coli" | replicon == "pSymA" | replicon == "pSymB") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    } else {
      tmp_end <- max_pos - oriC_pos
      new_pos[i] <- tmp_pos[i] + tmp_end
    }
  }
  tmp_pos <- new_pos
}
 
 
#now accounting for the bidirectionality. if things are between the start pos and
#the terminus then they will stay as the same position. If not, then they will be
#changed to a new position starting at 1 and going to the terminus
new_pos2 <- tmp_pos
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] > terminus) {
      new_pos2[i] <- max_pos - tmp_pos[i]
    } else {
    }
  }
  tmp_pos <- new_pos2
   
  print("max tmp_pos")
  max(tmp_pos)
}
 

if (bac_name == "Streptomyces") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    }#if btwn origin and end of genome 
    if (tmp_pos[i] <= oriC_pos) {
      new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
    }#if btwn origin and beginning of genome
    if (tmp_pos[i] == oriC_pos) {
      new_pos[i] <- 0
    }#if equal to origin
  }
  tmp_pos <- new_pos
}


merged_data$bidir_pos <- tmp_pos
colnames(merged_data) <- c("branch", "from1", "to1", "position", "change", "from2", "to2", "prob", "bidir_pos")
head(merged_data)
print("bidir_pos max")
max(merged_data$bidir_pos)
print("bidir_pos min")
min(merged_data$bidir_pos)
 
#################################################################################
#################################################################################
## SUBS GRAPH
#################################################################################
#################################################################################

#new min and max for the scaled positions (rounded to chunk len
#specified in args)
chunklen <- as.numeric(args[7])
print("chunklen")
chunklen
nmin_pos <- round_any(min(merged_data$bidir_pos), chunklen, f=floor)
print("min (rounded)")
nmin_pos
nmax_pos <- max(merged_data$bidir_pos)
print("max (rounded)")
nmax_pos

#expty vec to hold the rows that split the dat into X kb chunks
subs_rows_to_split_dat <- vector()
#create vector of chunks
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling) 
chunk_len_of_genome
chunks <- seq(nmin_pos, chunk_len_of_genome, chunklen)
chunks
##checking if the lengths are equal for the later bind
#if (length(unlist(tot_gene_subs_10kb)) == length(chunks)) {
#    chunks_pos_NOzero <- chunks   
#} else { 
    chunks_pos_NOzero <- chunks[which(chunks != 0)]
#}
if (bac_name == "Streptomyces") {
    chunks_pos_NOzero <- chunks
}
chunks_pos_NOzero

#add in fake data for chunks so it will split properly
len_of_chunks <- length(chunks_pos_NOzero)
fbranch <- rep(as.numeric(8), len_of_chunks)
ffrom1 <- rep(as.numeric(9), len_of_chunks)
fto1 <- rep(as.numeric(3), len_of_chunks)
fchange <- rep(as.numeric(0), len_of_chunks)
ffrom2 <- rep("T", len_of_chunks)
fto2 <- rep("C", len_of_chunks)
fprob <- rep(as.numeric(1.00), len_of_chunks)
fake_dat <- cbind(fbranch, ffrom1, fto1, chunks_pos_NOzero, fchange, ffrom2, fto2, fprob, chunks_pos_NOzero)
colnames(fake_dat) <- c("branch", "from1", "to1", "position", "change", "from2", "to2", "prob", "bidir_pos")
head(fake_dat)
head(merged_data)
#add fake data to orig data
merged_data <- rbind(merged_data,fake_dat)
head(merged_data)
tail(merged_data)
merged_data$branch <- as.numeric(merged_data$branch)
merged_data$from1 <- as.numeric(merged_data$from1)
merged_data$to1 <- as.numeric(merged_data$to1)
merged_data$position <- as.numeric(merged_data$position)
merged_data$change <- as.numeric(merged_data$change)
merged_data$prob <- as.numeric(merged_data$prob)
merged_data$bidir_pos <- as.numeric(merged_data$bidir_pos)


#order the data by positions, only if its NOT strep
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti"| bac_name == "Streptomyces") {
	ord_pos <- order(merged_data$bidir_pos)
	merged_data$bidir_pos <- merged_data$bidir_pos[ord_pos]
	merged_data$branch <- merged_data$branch[ord_pos]
	merged_data$from1 <- merged_data$from1[ord_pos]
	merged_data$to1 <- merged_data$to1[ord_pos]
	merged_data$position <- merged_data$position[ord_pos]
	merged_data$change <- merged_data$change[ord_pos]
	merged_data$from2 <- merged_data$from2[ord_pos]
	merged_data$to2 <- merged_data$to2[ord_pos]
	merged_data$prob <- merged_data$prob[ord_pos]
}
print("SAVED BIDIRECTIONAL DATA TO FILE")
write.table(merged_data, 'bidirectional_data.csv', sep = "\t")

print("checking order")
head(merged_data)
tail(merged_data)

for (i in chunks) {
  subs_rows <-
which(abs(merged_data$bidir_pos-i)==min(abs(merged_data$bidir_pos-i)))
# finding the closest number to each 10kb without going over it
  subs_max_row <- max(subs_rows)
  subs_rows_to_split_dat <- c(subs_rows_to_split_dat, subs_max_row)
}#for

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_subs <- merged_data$change
list_dat_sets_just_subs <- split(data_just_subs, findInterval(1:nrow(merged_data), subs_rows_to_split_dat))
print("split subs into 10kb chunks")


##########
#add up all subs in each 10kb section
##########
list_tot_add_subs_10kb <- lapply(list_dat_sets_just_subs, sum)
tot_gene_subs_10kb <- as.data.frame(matrix(unlist(list_tot_add_subs_10kb), byrow = F))
print("head/tail of total_subs")
head(tot_gene_subs_10kb)
tail(tot_gene_subs_10kb)
nrow(tot_gene_subs_10kb)
length(chunks_pos_NOzero)
class(tot_gene_subs_10kb)
#remove last element if lengths do not match (for cbind later)
if (nrow(tot_gene_subs_10kb) == length(chunks_pos_NOzero)) {
    tot_gene_subs_10kb <- tot_gene_subs_10kb 
} else { 
    tot_gene_subs_10kb <- tot_gene_subs_10kb[-nrow(tot_gene_subs_10kb),] 
    tot_gene_subs_10kb <- tot_gene_subs_10kb[-1] 
}

tot_gene_subs_10kb <- as.data.frame(cbind(tot_gene_subs_10kb, chunks_pos_NOzero))
head(tot_gene_subs_10kb)
tail(tot_gene_subs_10kb)
print("tot_gene_subs_10kb")
write.csv(tot_gene_subs_10kb)
#################
##########
#get total number of sites in each 10kb section
##########
list_freq <- lapply(list_dat_sets_just_subs, table)
list_freq
list_freq_sum <- lapply(list_freq, sum)
list_freq_sum
tot_num_sites <- as.data.frame(matrix(unlist(list_freq_sum), byrow=F))
tot_num_sites
length(tot_num_sites)
#remove last element if lengths do not match (for cbind later)
if (nrow(tot_num_sites) == length(chunks_pos_NOzero)) {
    tot_num_sites <- tot_num_sites 
} else { 
    tot_num_sites <- tot_num_sites[-nrow(tot_num_sites),] 
    tot_num_sites <- tot_num_sites[-1] 
}
freq_dat <- as.data.frame(cbind(tot_num_sites,chunks_pos_NOzero))
print("freq_dat")
head(freq_dat)
tail(freq_dat)
write.csv(freq_dat,"freq_dat.csv")
##########
#weighted dataframe
##########
weighted_data <- cbind(tot_gene_subs_10kb,freq_dat)
colnames(weighted_data) <- c("tot_gene_subs_10kb", "chunks_pos_NOzero", "tot_num_sites", "chunks_pos_NOzero")
head(weighted_data)
tail(weighted_data)
length(weighted_data$tot_gene_subs_10kb)
length(weighted_data$tot_num_sites)
weighted_data$weighted_subs <- weighted_data$tot_gene_subs_10kb / weighted_data$tot_num_sites
head(weighted_data)
tail(weighted_data)

weighted_data <- weighted_data[-2]
head(weighted_data)
tail(weighted_data)
write.csv(weighted_data,"weighted_dat.csv")

######################################################
# OUTLIERS function
######################################################
#calculating outliers (function)
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  #response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
 # if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
#  } else{
 #   cat("Nothing changed", "n")
  #  return(invisible(var_name))
  #}
}

print("##################################")
print("## CALCULATING OUTLIERS FOR WEIGHTED  ##")
print("##################################")
#setting new df to remove outliers
weighted_outliers <- weighted_data
#to run function:
#outlierKD(data_frame, column_with_data)
#finding outlier bars in the weighted subs data
outlierKD(weighted_outliers, weighted_subs)
print("checking removed outliers df")
head(weighted_outliers)
tail(weighted_outliers)

weighted_outlier_rows <- which(is.na(weighted_outliers$weighted_subs))
print("rows where the bars are outliers")
weighted_outlier_rows


weighted_out_pos <- vector()
for (i in weighted_outlier_rows) {
  weighted_out_pos <- c(weighted_out_pos, weighted_outliers$chunks_pos_NOzero[i])
}
print("outlier genomic positions (bar position)")
weighted_out_pos

print("##########################################")
print("##### 10kb lin reg WEIGHTED")
print("##########################################")
### COD ###
lm_exp_10kb_weighted_cod <- lm(weighted_subs ~abs(chunks_pos_NOzero), data=weighted_outliers)
summary(lm_exp_10kb_weighted_cod)
#########################################################################
print("##################################")
print("## CALCULATING OUTLIERS FOR NON-WEIGHTED  ##")
print("##################################")
#setting new df to remove outliers
nonweighted_outliers <- weighted_data
outlierKD(nonweighted_outliers, tot_gene_subs_10kb)
print("checking removed outliers df")
head(nonweighted_outliers)
tail(nonweighted_outliers)

nonweighted_outlier_rows <- which(is.na(nonweighted_outliers$tot_gene_subs_10kb))
print("rows where the bars are outliers")
nonweighted_outlier_rows


nonweighted_out_pos <- vector()
for (i in nonweighted_outlier_rows) {
  nonweighted_out_pos <- c(nonweighted_out_pos, nonweighted_outliers$chunks_pos_NOzero[i])
}
print("outlier genomic positions (bar position)")
nonweighted_out_pos

print("##########################################")
print("##### 10kb lin reg NON WEIGHTED")
print("##########################################")
### COD ###
lm_exp_10kb_nonweighted_cod <- lm(tot_gene_subs_10kb ~abs(chunks_pos_NOzero), data=nonweighted_outliers)
summary(lm_exp_10kb_nonweighted_cod)
print("##########################################")
print("##### 10kb lin reg on number of sites")
print("##########################################")
lm_num_sites <- lm(tot_num_sites ~ abs(chunks_pos_NOzero), data=weighted_data)
summary(lm_num_sites)
print("##################################")
print("## REMOVING OUTLIERS FOR log reg on subs per site  ##")
print("##################################")
print("weighted_out_pos")
weighted_out_pos
#setting new df to remove outliers
merged_outliers <- merged_data
merged_outliers$Outlier <- rep(0,length(merged_outliers$bidir_pos))
head(merged_outliers)
tail(merged_outliers)
print("outliers = 1")
merged_outliers[which(merged_outliers$Outliers == 1),]
for (i in weighted_out_pos) {
low_val <- i - 10000
high_val <- i
merged_outliers$Outlier[which(merged_outliers$bidir_pos > low_val & merged_outliers$bidir_pos <= high_val)] <- 1
}
print("test outliers removed")
merged_no_outliers <- merged_outliers[which(merged_outliers$Outlier==0),]
merged_no_outliers[which(merged_no_outliers$Outlier ==1),]
print("##########################################")
print("##### LOGISTIC REGRESSION on subs per site")
print("##########################################")
logistic_reg_model=glm(change~abs(bidir_pos),family=binomial,merged_no_outliers, control = list(maxit = 1000)) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
print(summary(logistic_reg_model))
print("##################################")
print("## CALCULATING AVERAGE NUMBER OF SUBS PER SITE  ##")
print("##################################")
#first making other values NA if they are in an outlier bar
weighted_outliers$tot_gene_subs_10kb[weighted_outlier_rows] <- NA
weighted_outliers$tot_num_sites[weighted_outlier_rows] <- NA
print("test weighted_outliers")
weighted_outliers

avg_num_of_sites_per_10kb <- sum(weighted_outliers$tot_gene_subs_10kb, na.rm=TRUE)
print("total subs")
avg_num_of_sites_per_10kb
print("total sites")
tot_sites <- sum(weighted_outliers$tot_num_sites, na.rm=TRUE)
print("average subs per site")
avg_num_of_sites_per_10kb / tot_sites

print("#####################################################")
print("# log reg on 20kb near origin (supplementary test)")
print("#####################################################")
#kb away from beginning or end
cutoff <- as.numeric(args[9])
print("cutoff value")
cutoff
#calculate average number of data pts per 10Kbp
avg_pts <- mean(weighted_data$tot_num_sites)
print("avg_pts")
avg_pts
avg_pts_cutoff <- avg_pts * (cutoff/10000)
print("avg_pts_cutoff")
avg_pts_cutoff
#make sure we are dealing with the absolute pos so that it works for
#strep
merged_no_outliers$abs_pos <- abs(merged_no_outliers$bidir_pos)
print("abs value column")
head(merged_no_outliers)
#getting row num of min value
min_row_num <- min(which(merged_no_outliers$abs_pos == min(merged_no_outliers$abs_pos, na.rm=TRUE)))
print("min_row_num")
min_row_num
#getting seq of row nums for points
cutoff_row_nums <- seq(min_row_num,(min_row_num + avg_pts_cutoff),1)
head(cutoff_row_nums)
tail(cutoff_row_nums)
#getting actual data pts that are in this range
kb_near_ori <- merged_no_outliers[cutoff_row_nums,]

colnames(kb_near_ori) <- c("branch", "from1", "to1", "position", "change", "from2", "to2", "prob", "bidir_pos","Outlier", "abs_pos")
print("kb_near_orii head,tail")
head(kb_near_ori)
tail(kb_near_ori)
print("kb_near_orii max,min")
max(kb_near_ori$abs_pos)
min(kb_near_ori$abs_pos)
near_ori_log_reg=glm(change~abs_pos,family=binomial(link="logit"),kb_near_ori, control = list(maxit = 1000)) # run a logistic regression model (in #this case, generalized linear model with logit link). see ?glm
print(summary(near_ori_log_reg))
######### total number of subs in 10kb region
head(kb_near_ori)
tail(kb_near_ori)
length(which(kb_near_ori$change == 1))
length(kb_near_ori$change)
tot_subs <- sum(kb_near_ori$change)
tot_subs
######### total number of data pts in 20kb region
tot_dat_pts <- length(kb_near_ori$change)
print("######## total number of subs / number of sites (20kb)")
tot_subs / tot_dat_pts

print("#####################################################")
print("# log reg on 20kb FAR origin (supplementary test)")
print("#####################################################")
head(merged_no_outliers)
end_of_genome <- max(max(merged_no_outliers$abs_pos,na.rm = TRUE))
end_of_genome
#getting row num of min value
max_row_num <- which(merged_no_outliers$abs_pos == end_of_genome)
#getting seq of row nums for points
cutoff_row_nums <- seq((max_row_num - avg_pts_cutoff),max_row_num,1)
head(cutoff_row_nums)
tail(cutoff_row_nums)
#getting actual data pts that are in this range
kb_far_ori <- merged_no_outliers[cutoff_row_nums,]

colnames(kb_near_ori) <- c("branch", "from1", "to1", "position", "change", "from2", "to2", "prob", "bidir_pos","Outlier", "abs_pos")
print("kb_far_ori")
head(kb_far_ori)
tail(kb_far_ori)
max(kb_far_ori$bidir_pos)
min(kb_far_ori$bidir_pos)
far_ori_log_reg=glm(change~abs_pos,family=binomial(link="logit"),kb_far_ori, control = list(maxit = 1000)) # run a logistic regression model (in #this case, generalized linear model with logit link). see ?glm
print(summary(far_ori_log_reg))
######### total number of subs in 10kb region
head(kb_far_ori)
tail(kb_far_ori)
length(which(kb_far_ori$change == 1))
length(kb_far_ori$change)
tot_subs <- sum(kb_far_ori$change)
tot_subs
######### total number of data pts in 20kb region
tot_dat_pts <- length(kb_far_ori$change)
print("######## total number of subs / number of sites (20kb)")
tot_subs / tot_dat_pts

print("####################")
print("# DOING OTHER END OF CHROM ARM FOR STREP (FAR FROM ORI) #")
print("####################")
if ( any(merged_no_outliers$bidir_pos < 0)) {
    print("short arm terminus pos")
    end_of_short_arm <- min(merged_no_outliers$bidir_pos)
    print(end_of_short_arm)
#getting row num of min value
min_row_num <- min(which(merged_no_outliers$bidir_pos ==end_of_short_arm))
#getting seq of row nums for points
cutoff_row_nums <- seq(min_row_num,(min_row_num + avg_pts_cutoff),1)
head(cutoff_row_nums)
tail(cutoff_row_nums)
#getting actual data pts that are in this range
kb_far_short_arm <- merged_no_outliers[cutoff_row_nums,]
kb_far_short_arm
    print("data near short arm end (head, tail)")
    print(head(kb_far_short_arm))
    print(tail(kb_far_short_arm))
    far_short_arm_log_reg=glm(change~bidir_pos,family=binomial(link="logit"),kb_far_short_arm, control = list(maxit = 1000)) # run a logistic regression model (in #this case, generalized linear model with logit link). see ?glm 
    print(summary(far_short_arm_log_reg))
}


#colours for histogram (outliers and reg bars)
#out_pos_colours <- rep("#BD93BD", length(rm_outliers$weighted_subs))
#out_pos_colours <- rep("#CFAFCF", length(rm_outliers$weighted_subs))
weighted_out_pos_colours <- rep("black", length(weighted_outliers$weighted_subs))
for (i in weighted_outlier_rows) {
    #colour for an outlier
#  out_pos_colours[i] <- "#F2EDEB"
#  out_pos_colours[i] <- "#CAD2C5"
  weighted_out_pos_colours[i] <- "#D1D1D1"
#  out_pos_colours[i] <- "#E8D7F1"
#  out_pos_colours[i] <- "#C1CAD6"
}

weighted_out_pos_colours


######################################################
# PLOT
######################################################
# WEIGHTED SUBS PER SITE
###########################
options(scipen=3)
weighted_data$chunks_pos_NOzero <- weighted_data$chunks_pos_NOzero / 1000000
subs_hist <- (ggplot(weighted_data, aes(x = chunks_pos_NOzero, y = weighted_subs)) 
  + geom_histogram(stat = "identity", fill= weighted_out_pos_colours) 
#  geom_histogram(stat = "identity", fill= "#FE5F55") 
  + labs(x = "Distance from the Origin of Replication (Mbp)", y = "Substitutions per 10Kbp") 
  + ggtitle(mytitle) 
   + geom_vline(xintercept = 0, colour = "red")
 + scale_x_continuous(labels = function(x) ifelse(x ==0, "0", x))
 + scale_y_continuous(labels = function(x) ifelse(x ==0,"0", x))
)
###########################
pdf(paste(out_file_name,"_weighted_subs_bidirectionality_colour.pdf", sep=""))
subs_hist
#grid.newpage()
#grid.draw(rbind(ggplotGrob(gene_num_hist), ggplotGrob(exp_bar_top), size = "last"))
## next line adds border
#grid.rect(width = 0.99, height = 0.99, gp = gpar(lwd = 2, col = "black", fill = NA))
dev.off()
###########################

#####################
# distribution of total number of sites (to check for potential
# missing data)
############################
options(scipen=3)
weighted_data$chunks_pos_NOzero <- weighted_data$chunks_pos_NOzero / 1000000
subs_hist <- (ggplot(weighted_data, aes(x = chunks_pos_NOzero, y = tot_num_sites)) 
  + geom_histogram(stat = "identity", fill= weighted_out_pos_colours) 
  + labs(x = "Distance from the Origin of Replication (Mbp)", y = "Total Number of Data Points per 10Kbp") 
  + ggtitle(mytitle) 
   + geom_vline(xintercept = 0, colour = "red")
 + scale_x_continuous(labels = function(x) ifelse(x ==0, "0", x))
 + scale_y_continuous(labels = function(x) ifelse(x ==0,"0", x))
)
###########################
pdf(paste(out_file_name,"_total_num_sites_graph.pdf", sep=""))
subs_hist
#grid.newpage()
#grid.draw(rbind(ggplotGrob(gene_num_hist), ggplotGrob(exp_bar_top), size = "last"))
## next line adds border
#grid.rect(width = 0.99, height = 0.99, gp = gpar(lwd = 2, col = "black", fill = NA))
dev.off()
###########################
