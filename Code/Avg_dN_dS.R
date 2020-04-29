#will calculate the weighted mean for dN and dS
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME GENOME_POSITION_DATAFRAME_FILE LENGTH_OF_GENOME ORIC_POS TERMINUS_POS
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
library(tibble)
library(tidyr)

#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
bac_name <- as.character(args[2])
replicon <- as.character(args[3])
genome_pos_file <- as.character(args[4])
chunk <- as.numeric(args[9])

mytitle <- substitute(italic(bac_name)~replicon, list(bac_name=bac_name, replicon=replicon))
#set theme
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            #change size of facet header text
            #theme(strip.text = element_text(size =10.49)) +
            theme(strip.text = element_text(size =18)) +
            theme(plot.title = element_text(hjust = 0.5, size = 18),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_line(colour = "grey90", size = 0.2),
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5),
                  panel.spacing = unit(0.25, "lines"),
                  axis.text=element_text(size=18),
                  axis.title = element_text(size = 18),
                  #plot margins
                  #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(size=18),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  legend.position="none")
                  #legend.position="top")
)


options("scipen"=100, "digits"=5)
datafile <- file_name
data <- read.table(datafile, sep = "\t", header = TRUE)
#now incorpreting the genome position info
genome_pos_datafile <- genome_pos_file
genome_pos_data <- read.table(genome_pos_datafile, sep = "\t", header = TRUE)
#merge the dataframes
merged_data <- merge(genome_pos_data, data, by=c("block","gene","sec","sec_two"))
print("merged_data")
head(merged_data)
#dealing with dS = 0.0001 values that should actually be zero (PAML
#does weird things when there are no subs)
merged_data$dS[which(merged_data$dS == 0.0001)] <- 0
#dealing with weird omega values
merged_data$omega[which(merged_data$dS == 0)] <- NA
merged_data$omega[which(merged_data$dN== 0 & merged_data$dS == 0)] <- NA
merged_data$omega[which(merged_data$dN== 0)] <- 0
no_NA_dat <- na.omit(merged_data)
no_NA_dat[which(no_NA_dat$omega >= 900),]
merged_ordered <- merged_data[order(merged_data$midpoint),]
print("head merged_ordered")
head(merged_ordered)
print("head no_NA_dat")
head(no_NA_dat)
#save dataframe so that this can be put into a supplementary table for
#paper
table_dat <- select(merged_ordered, gene_name,dN,dS,omega)
print("head table_dat")
head(table_dat)
##########################
#supplementary file table
##########################
##calculate the average per gene so data frame has only one pt per gene
test_dat <- table_dat %>%
  group_by(gene_name) %>%
  summarise_at(.vars = c("dN", "dS", "omega"), .funs = mean)
file_name <- as.character(args[8])
class(test_dat)
write.csv(test_dat, paste(file_name,"_supp_table2.csv",sep=""),row.names= FALSE, quote = FALSE)
######################################################
#NOT FOR SUPPLEMENT: to associate pos and gene name
##calculate the average per gene so data frame has only one pt per gene
sel_dat <- select(merged_ordered, gene_name,dN,dS,omega, midpoint)
head(sel_dat)
name_pos_dat <- sel_dat %>%
  group_by(gene_name) %>%
  summarise_at(.vars = c("dN", "dS", "omega", "midpoint"), .funs = mean)
#################################################################################
#
#
options("scipen"=100, "digits"=10)
################################################################################
#ORIGIN SCALING AND BIDIRECTIONALITY
################################################################################
#first scaling things to the origin (if necessary)
max_pos <- as.numeric(args[5])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[6])
print("oriC")
oriC_pos
terminus <- as.numeric(args[7])
print("ter")
terminus
new_pos <- merged_data$midpoint
tmp_pos <- merged_data$midpoint
print("MIN POS")
min(merged_data$midpoint)
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


merged_data$midpoint <- tmp_pos
colnames(merged_data) <- c("block", "gene", "sec","sec_two", "start", "end", "tmp_pos","gene_name", "dS", "dN", "omega", "sec_len")
head(merged_data)
max(merged_data$tmp_pos)
min(merged_data$tmp_pos)
 
#now do the same for the proper omega values
max_pos <- as.numeric(args[5])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[6])
print("oriC")
oriC_pos
terminus <- as.numeric(args[7])
print("ter")
terminus
new_pos <- no_NA_dat$midpoint
tmp_pos <- no_NA_dat$midpoint
 
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


no_NA_dat$midpoint <- tmp_pos
colnames(no_NA_dat) <- c("block", "gene", "sec","sec_two", "start", "end", "tmp_pos","gene_name", "dS", "dN", "omega", "sec_len")
head(no_NA_dat)
max(no_NA_dat$tmp_pos)
min(no_NA_dat$tmp_pos)
################################################################################
# now doing this for selection data set that has the position and gene
# name in it:
new_pos <- name_pos_dat$midpoint
tmp_pos <- name_pos_dat$midpoint
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
   
}
 
print("max tmp_pos")
max(tmp_pos)

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


print("tibble?")
class(name_pos_dat)
as.data.frame(name_pos_dat)
class(name_pos_dat)
name_pos_dat <- add_column(name_pos_dat, tmp_pos)
print("max after bidir")
max(name_pos_dat$tmp_pos)
print("head name_pos_dat")
head(name_pos_dat)
 
class(name_pos_dat)
as.data.frame(name_pos_dat)
class(name_pos_dat)
write.csv(name_pos_dat, paste(file_name,"_selection_data.csv",sep=""),row.names= FALSE, quote = FALSE)


pdf("missing_data_check_dot_plot.pdf")  
#omega graph
(dot_plot <- ggplot(data = name_pos_dat, mapping = aes(x=midpoint,y=dN)) 
  +geom_point() 
  +xlab("Distance From the Origin of Replication (bp)")
  + ggtitle(mytitle)
)
dev.off()


print("test_data_no_zeros")
test_data_no_zeros <- name_pos_dat %>%
			filter(omega !=0)
test_data_no_zeros[(test_data_no_zeros$omega ==0),]
####################################################
print("number of data pts that have omega =0")
head(merged_data)
dat_zero <- merged_data %>%
		filter(omega ==0)
omeg_zero <- dat_zero$omega
head(omeg_zero)
length(omeg_zero)
print("total data pts")
length(merged_data$omega)
####################################################


####################################################
# OUTLIERS
####################################################
#calculate outliers function
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  #m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  #boxplot(var_name, main="With outliers")
  #hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  #mo <- mean(outlier,na.rm = T)
  #var_name <- ifelse(var_name %in% outlier, NA, var_name)
  var_name <- ifelse(var_name %in% outlier, "outlier", var_name)
  #boxplot(var_name, main="Without outliers")
  #hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  #title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  #cat("Mean of the outliers:", round(mo, 2), "n")
  #m2 <- mean(var_name, na.rm = T)
  #cat("Mean without removing outliers:", round(m1, 2), "n")
  #cat("Mean if we remove outliers:", round(m2, 2), "n")
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

print("boxplot stats")
boxplot.stats(merged_data$omega)

####################################################
print("calculating outliers for each selection category")
test_outliers_df_dN <- merged_data
head(test_outliers_df_dN)
print("outliers for dN")
#outlierKD(test_outliers_df_dN, dN)
print("outliers for dS")
#outlierKD(test_outliers_df_dN, dS)
#right now we are ONLY using omega as the classification for outliers
print("outliers for omega")
outlierKD(test_outliers_df_dN, omega)
print("TEST OUTLIERS")
head(test_outliers_df_dN)
write.csv(test_outliers_df_dN, paste(file_name,"test_data_zeros_outliers",sep=""))
#test_outliers_df_dN

#creating new column that treats an outlier in any row as an outlier overall
test_outliers_df_dN$Outlier <- ifelse(test_outliers_df_dN$dN == "outlier", 1, 0)
test_outliers_df_dN <- within(test_outliers_df_dN, {
				f <- dS == 'outlier'
				Outlier[f] <- 1
})
test_outliers_df_dN <- within(test_outliers_df_dN, {
				f <- omega == 'outlier'
				Outlier[f] <- 1
})

print("OUTLIERS df")
merged_data <- cbind(merged_data,test_outliers_df_dN$Outlier)
head(merged_data)
colnames(merged_data)[13] <- "Outlier"
colnames(merged_data)[7] <- "pos"

print("summary of outliers")
print("total data points")
length(merged_data$Outlier)
print("total outliers")
length(which(merged_data$Outlier == 1))
print("percent of total data that is outliers")
length(which(merged_data$Outlier == 1)) /length(merged_data$Outlier)

print("summary of dN outliers")
print("total dN = 0")
length(which(merged_data$dN == 0))
print("total dN < 0")
length(which(merged_data$dN < 0))
print("percent of total data that is dN = 0")
length(which(merged_data$dN == 0)) /length(merged_data$Outlier)
print("percent of total data that is dS = 0")
length(which(merged_data$dS == 0)) /length(merged_data$Outlier)
print("percent of total data that is omega = 0")
length(which(merged_data$omega == 0)) /length(merged_data$Outlier)
####################################################
#reshape data into a format that ggplot will recognize
merged_data <- merged_data %>%
		gather(class,value, dS:omega)
print("gathered")
head(merged_data)
merged_data$bac <- rep(file_name,length(merged_data$value))
head(merged_data)

#save dataframe so we can use it to put the box plots all together!
print("MIN POS")
min(merged_data$pos)
write.csv(merged_data, paste(file_name,"_box_plot_dat2.csv",sep=""))
print("just printed csv")


################################################################################
#DN, DS, OMEGA ANALYSIS
################################################################################
print("MERGED DAT")
merged_data
no_NA_dat <- na.omit(merged_data)
print("MERGED DAT NO NA")
no_NA_dat
#below will make a separate column for dN, dS, omega by separating the
#class and values columns
no_NA_dat <- no_NA_dat %>%
		spread(class,value)
no_NA_dat <- na.omit(no_NA_dat)
hilight_outliers <- no_NA_dat[no_NA_dat$Outlier == 1,]
print("hilight_outliers")
head(hilight_outliers)
no_NA_dat <- no_NA_dat[no_NA_dat$Outlier == 0,]
print("total dN")
length(no_NA_dat$dN)
print("total dN = 0")
length(which(no_NA_dat$dN == 0))
print("TEST SPREAD")
no_NA_dat
print("TEST NA")
no_NA_dat[which(is.na(no_NA_dat)),]

print("WEIGHTED average")
#dS whole genome weighted average
print("gene/whole genome dS weighted average")
weighted.mean(no_NA_dat$dS, no_NA_dat$sec_len)

#dN whole genome weighted average
print("gene/whole genome dN weighted average")
weighted.mean(no_NA_dat$dN, no_NA_dat$sec_len)

#omega whole genome weighted average
print("gene/whole genome omega weighted average")
weighted.mean(no_NA_dat$omega, no_NA_dat$sec_len)
print("########################################################")
print("NON-WEIGHTED average")
#dS whole genome non-weighted average
print("gene/whole genome dS average")
mean(no_NA_dat$dS)

#dN whole genome non-weighted average
print("gene/whole genome dN average")
mean(no_NA_dat$dN)

#omega whole genome non-weighted average
print("gene/whole genome omega average")
mean(no_NA_dat$omega)

################################################################################
print("GENOME POSITION AND DN, DS, OMEGA DISTRIBUTION")
################################################################################
head(no_NA_dat)
tail(no_NA_dat)
min_pos <- min(no_NA_dat$pos)
print("min_pos")
min_pos
max_pos <- max(no_NA_dat$pos)
print("max_pos")
max_pos
round_max_pos <- round_any(max_pos, chunk, f=ceiling)
print("round_max_pos")
round_max_pos
class(no_NA_dat)
class(min_pos)
class(round_max_pos)
chunks_len <- length(seq(min_pos,round_max_pos,chunk))
#order the data by position
print("order the data by position")
no_NA_dat2 <- no_NA_dat[order(no_NA_dat$pos),]
head(no_NA_dat2)
tail(no_NA_dat2)
all_dat_no_bins <- no_NA_dat2
# creating the bins
print("breaks")
break_seq <- c(seq(min_pos -1,round_max_pos+1,chunk),round_max_pos)
break_seq
no_NA_dat2$group<-cut(no_NA_dat2$pos, breaks = break_seq,labels = seq(1:chunks_len))
no_NA_dat2 <- na.omit(no_NA_dat2)
# Finding the mean value for all bins
print("finding mean of each bin")
print("max_pos bin data")
max(no_NA_dat2$pos)
# dN
no_NA_dat2$dNavg = rep(0,length(no_NA_dat2$pos))
for (i in 1:chunks_len) {
y = no_NA_dat2$dN[(no_NA_dat2$group == i)]
no_NA_dat2$dNavg[(no_NA_dat2$group == i)] = mean(no_NA_dat2$dN[(no_NA_dat2$group == i)])
}
# dS
no_NA_dat2$dSavg = rep(0,length(no_NA_dat2$pos))
for (i in 1:chunks_len) {
y = no_NA_dat2$dS[(no_NA_dat2$group == i)]
no_NA_dat2$dSavg[(no_NA_dat2$group == i)] = mean(no_NA_dat2$dS[(no_NA_dat2$group == i)])
}
# omega
no_NA_dat2$omegaavg = rep(0,length(no_NA_dat2$pos))
for (i in 1:chunks_len) {
y = no_NA_dat2$omega[(no_NA_dat2$group == i)]
no_NA_dat2$omegaavg[(no_NA_dat2$group == i)] = mean(no_NA_dat2$omega[(no_NA_dat2$group == i)])
}
# position (mdeian)
no_NA_dat2$pos_mid = rep(0,length(no_NA_dat2$pos))
for (i in 1:chunks_len) {
y = no_NA_dat2$pos[(no_NA_dat2$group == i)]
no_NA_dat2$pos_mid[(no_NA_dat2$group == i)] = mean(no_NA_dat2$pos[(no_NA_dat2$group == i)])
}
print("all no_NA_dat2")
write.csv(no_NA_dat2, "no_NA_dat2_TEST.csv")
no_NA_dat2
# The plot:
#dN graph
print("dN_dist.pdf")  
#pdf("dN_dist.pdf")  
#(dN_graph <- ggplot(data = no_NA_dat2, mapping = aes(x=group,y=dN)) 
#  +geom_jitter(aes(color="#F3D9B1"),alpha=0.7, width = 0.2) 
#  +geom_point(aes(y = dNavg, color = "#C29979")) 
#  +geom_smooth(aes(x = as.numeric(group), y = dNavg,color = "#6494AA"),span=0.5, method = "loess") 
#  +guides(color=FALSE)
#  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
#  + scale_x_discrete(breaks = c(1,10,20,30,40,50))
#  +scale_color_manual(values = c("#6494AA","#C29979","#F3D9B1"))
#  +ylab("dN") 
#  +xlab("Distance From the Origin of Replication (100Kbp)")
#  + ggtitle(mytitle)
#)
#dev.off()
#dS graph
print("dS_dist.pdf")  
pdf("dS_dist.pdf")  
(dS_graph <- ggplot(data = no_NA_dat2, mapping = aes(x=group,y=dS)) 
  +geom_jitter(aes(color="#F3D9B1"),alpha=0.7, width = 0.2) 
  +geom_point(aes(y = dSavg, color = "#C29979")) 
  +geom_smooth(aes(x = as.numeric(group), y = dSavg,color = "#6494AA"),span=0.5, method = "loess") 
  +guides(color=FALSE)
  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
  + scale_x_discrete(breaks = c(1,10,20,30,40,50))
  +scale_color_manual(values = c("#6494AA","#C29979","#F3D9B1"))
  +ylab("dS") 
  +xlab("Distance From the Origin of Replication (100Kbp)")
  + ggtitle(mytitle)
)
dev.off()
print("omega_dist.pdf")  
#pdf("omega_dist.pdf")  
##omega graph
#(omega_graph <- ggplot(data = no_NA_dat2, mapping = aes(x=group,y=omega)) 
#  +geom_jitter(aes(color="#F3D9B1"),alpha=0.7, width = 0.2) 
#  +geom_point(aes(y = omegaavg, color = "#C29979")) 
#  +geom_smooth(aes(x = as.numeric(group), y = omegaavg,color = "#6494AA"),span=0.5, method = "loess") 
#  +guides(color=FALSE)
#  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
#  + scale_x_discrete(breaks = c(1,10,20,30,40,50))
#  +scale_color_manual(values = c("#6494AA","#C29979","#F3D9B1"))
#  +ylab(expression(omega)) 
#  +xlab("Distance From the Origin of Replication (100Kbp)")
#  + ggtitle(mytitle)
#+ theme(
#      axis.title.y = element_text(size = 25),
#)
#)
#dev.off()
#####################
# facet graph of dN, dS, omega
#####################
print("facet")
print("get data in to long format")
#do this for outliers
hilight_outliers <- gather(hilight_outliers, class,value, dN:omega, factor_key=TRUE)
print("hilight_outliers")
head(hilight_outliers)
#do this for data without bins and averages
all_dat_no_bins <- gather(all_dat_no_bins, class,value, dN:omega, factor_key=TRUE)
print("all_dat_no_bins")
head(all_dat_no_bins)
#now doing it with data with averages
print("READ MEall_dat")
all_dat <- subset(no_NA_dat2,select=c(block,gene,sec,sec_two,start,end,pos,gene_name,sec_len,Outlier,bac,group,dNavg,dSavg,omegaavg,pos_mid))
head(all_dat)
names(all_dat)[names(all_dat) == 'dNavg'] <- 'dN'
names(all_dat)[names(all_dat) == 'dSavg'] <- 'dS'
names(all_dat)[names(all_dat) == 'omegaavg'] <- 'omega'
all_dat <- gather(all_dat, class, avg, dN:omega, factor_key=TRUE)
head(all_dat)
print("all_dat")
write.csv(all_dat, "all_dat_TEST.csv")
head(all_dat)
tail(all_dat)



print("scaling the position values so they look prettier")
all_dat$pos <- all_dat$pos / 1000000
all_dat$pos_mid <- all_dat$pos_mid / 1000000
hilight_outliers$pos <- hilight_outliers$pos / 1000000
all_dat_no_bins$pos <- all_dat_no_bins$pos / 1000000
print("max_pos graph dat")
max(all_dat_no_bins$pos)
# filter dataframe to get data to be highligheted
highlight_df <- subset(all_dat,select= c(group,class,pos_mid,avg))
head(highlight_df)
highlight_df <- subset(highlight_df,!duplicated(highlight_df))
write.csv(highlight_df, "highlight_df_TEST.csv")
highlight_df
print("changing positions so that they are in the middle of the bin")
names(highlight_df)[names(highlight_df) == 'pos_mid'] <- 'pos'
print("HILIGHTED DATA")
head(highlight_df)
tail(highlight_df)
print("hilight dat = 0")
highlight_df[which(highlight_df$avg ==0),]
pdf("dN_dS_omega_dist.pdf")  
#omega graph
(omega_graph <- ggplot(data = all_dat_no_bins, mapping = aes(x=pos,y=value)) 
#light beige
  +geom_point(aes(color="#F3D9B1"),alpha=0.7) 
#more of a mustard
  +geom_point(data = hilight_outliers,aes(y = value, color = "#BEBEBE"), alpha=0.7, shape=1) 
  +geom_point(data = highlight_df,aes(y = avg, color = "#C29979")) 
  +geom_point(data = highlight_df,aes(y = avg, color = "#C29979")) 
  +geom_smooth(data=highlight_df,aes(y = avg,color = "#6494AA"),span=0.5, method = "loess") 
  #adding vertical line for the origin
  + geom_vline(xintercept = 0, colour = "#939290")
+ facet_wrap(ncol = 1,~class, labeller=label_parsed, strip.position = "left")
#omega = 1 referene line only on omega graph
  + geom_hline(data = data.frame(yint=1,class="omega"), aes(yintercept =yint), linetype = "dashed", colour = "#939290")
  +guides(color=FALSE)
  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x ==0, "0", x),breaks=c(0.001,0.01,0.1, 1, 10))
  # make it so the x-axis ticks start right at zero
  #and remove trailing zeros
  + scale_x_continuous(expand = c(0.009, 0),labels = function(x) ifelse(x ==0, "0",x))
  +scale_color_manual(values = c("#6494AA","#BEBEBE","#C29979","#F3D9B1"))
  +xlab("Distance From the Origin of Replication (Mbp)")
  + ggtitle(mytitle)
  + theme(
      axis.title.y = element_blank(),
)
)
dev.off()
