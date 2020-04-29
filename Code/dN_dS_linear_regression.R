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

#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])


options("scipen"=100, "digits"=5)
datafile <- file_name
data <- read.table(datafile, sep = ",")
colnames(data) <- c("row","block","gene","sec","sec_two","start","end","pos","gene_name","sec_len","Outlier","class","value","bac")
head(data)
#data
split_bac <- split(data, data$bac)
split_bac_class <- lapply(split_bac, function(x) split(x, x$class))
length(split_bac_class)

print("################################")
print("omega > 1")
print("################################")
omeg <- lapply(split_bac_class, function(x) lapply(x, function(df){
							      df[which(df$class == "omega"),]
		}))
omeg_big_1 <- lapply(omeg, function(x) lapply(x, function(df){
							      df[which(df$value >= 1),]
		}))
omeg_big_1

####################################################################
print("creating df with no outliers")
data <- data[data$Outlier == 0,]
head(data)
print("check no outliers")
data[which(data$Outlier == 1),]

print("test sinoC df")
data[which(data$bac == "sinoC"),]

split_bac <- split(data, data$bac)
split_bac_class <- lapply(split_bac, function(x) split(x, x$class))



print("##########################")
print(file_name)
print("##########################")
##############LINEAR REGRESSION FOR EXP DAT###########
# run a reg linear regression on the exp "data
lm_results <- lapply(split_bac_class, function(x) lapply(x, function(df){
			 lm1 <- lm(value~pos, data=df)
			 summary(lm1)
		}))
lm_results

##################################
# ORDER BY POSITION SO THE NEXT PART ACTUALLY LOOKS AT THE 20 GENES
# NEAR THE ORI AND TER
##################################
ordered_dat <- lapply(split_bac_class, function(x) lapply(x, function(df){
							      df[order(df$pos),]
		}))

##################################
# LIN REG FOR 20kb near ori (supplementary test)
##################################
#kb away from beginning or end of genome
cutoff <- 40000

print("##########################")
print("NEAR ORI")
print("##########################")
# NEAR ORI #

kb_near_ori <- lapply(ordered_dat, function(x) lapply(x, function(df){
							     head(df, n=20) 
		}))
kb_near_ori[1]
kb_near_ori[6]
kb_near_ori[5]

# run a reg linear regression on the 20kb near ori data
kb_near_lm_results <- lapply(kb_near_ori, function(x) lapply(x, function(df){
			 lm1 <- lm(value~pos, data=df)
			 summary(lm1)
		}))
kb_near_lm_results

# FAR ORI #
print("##########################")
print("FAR ORI")
print("##########################")

kb_far_ori <- lapply(ordered_dat, function(x) lapply(x, function(df){
							     tail(df, n=20) 
		}))
kb_far_ori[1]
kb_far_ori[6]
kb_far_ori[5]
# run a reg linear regression on the 20kb far ori data
kb_far_lm_results <- lapply(kb_far_ori, function(x) lapply(x, function(df){
			 lm1 <- lm(value~pos, data=df)
			 summary(lm1)
		}))
kb_far_lm_results



print("##########################")
print("STREP NEAR ORI")
print("##########################")
closest<-function(xv,sv){
xv[which(abs(xv-sv)==min(abs(xv-sv)))] }
strep_kb_near_ori <- lapply(ordered_dat[6], function(x) lapply(x, function(df){
                                                             df[seq(max(which(df$pos== closest(df$pos,0))) -9,max(which(df$pos ==closest(df$pos,0)))+10,by=1),]
		}))
strep_kb_near_ori
# run a reg linear regression on the 20kb far ori data
strep_kb_near_lm_results <- lapply(strep_kb_near_ori, function(x) lapply(x, function(df){
			 lm1 <- lm(value~pos, data=df)
			 summary(lm1)
		}))
strep_kb_near_lm_results
# FAR ORI #
print("##########################")
print("STREP FAR ORI")
print("##########################")
strep_kb_far_ori <- lapply(ordered_dat, function(x) lapply(x, function(df){
							     rbind(head(df,n=10),tail(df,n=10)) 
		}))
strep_kb_far_ori[6]
# run a reg linear regression on the 20kb far ori data
strep_kb_far_lm_results <- lapply(strep_kb_far_ori, function(x) lapply(x, function(df){
			 lm1 <- lm(value~pos, data=df)
			 summary(lm1)
		}))
strep_kb_far_lm_results[6]




