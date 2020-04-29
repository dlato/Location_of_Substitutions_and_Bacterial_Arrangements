#will calculate the weighted mean for dN and dS
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME GENOME_POSITION_DATAFRAME_FILE LENGTH_OF_GENOME ORIC_POS TERMINUS_POS

library("RColorBrewer")
library("gplots")
library("ggplot2")
library(dplyr)
library("methods")
library(methods)
options("scipen"=100, "digits"=10)
library(lattice)

#set theme
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
	    #change size of facet header text
            theme(strip.text = element_text(size =10.49)) +
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

file_list <- list.files(pattern ="*_plot_dat.csv")
file_list
length(file_list)
class(file_list)
print("names")
names(file_list)
length(file_list)
class(file_list)
print("test")
all_data <- do.call(rbind,lapply(file_list,read.csv))
print("test")
head(all_data)
levels(all_data$bac)

num_of_plots <- length(levels(all_data$bac)) 
print("number of bac")
num_of_plots
colours_arr <- rep(c("#CBE896","#8EBEC8","#C6878F"),num_of_plots)
cols_l <- c("#CBE896","#8EBEC8","#C6878F")
colours_arr <- rep(c("#C6878F","#6494AA","#C29979"),num_of_plots)
cols_l <- c("#C6878F","#6494AA","#C29979")
# BLACK AND WHITE COLOURS
#colours_arr <- rep(c("#969696","#D9D9D9","#525252"),num_of_plots)
print("length of colours array")
length(colours_arr)
sinoC_dat <- all_data %>%
		filter(bac == "sinoC")
levels(all_data$bac) <- c("ecoli" = expression(paste(italic("E.coli"), " Chromosome")),
			  "bass" = expression(paste(italic("B.Subtilis"), " Chromosome")),
			  "strep" = expression(paste(italic("Streptomyces"), " Chromosome")),
			  "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
			  "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
			  "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))


h <- 1

#remove outliers from all_data dataframe
all_data2 <- all_data
# filter dataframe to get data to be highligheted
highlight_df <- all_data %>% 
             filter(Outlier==1)
print("HILIGHTED DATA")
head(highlight_df)

set.seed(1738);
p<-(ggplot(all_data2, aes(x=class, y=value, fill=class, color = class)) 
  + geom_jitter(position=position_jitter(0.2), size = 1,alpha = 0.4)
  + geom_violin(colour = "black")
  + facet_wrap(~bac, labeller=label_parsed)
  #omega = 1 referene line
  + geom_hline(yintercept= 1, linetype = "dashed",color = "#939290") 
  + xlab("") 
  + ylab("Value")
  #proper colours for strip part of plot
  + scale_color_manual(values=cols_l,labels = c(" dN", " dS", expression(omega)))
  + scale_fill_manual(values=cols_l,labels = c(" dN", " dS", expression(omega)))
  #log transform, remove trailing zeros, custom breaks
  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x ==0, "0", x),breaks=c(0.0001,0.001,0.01,0.1, 1, 10,100)) 
)

pdf("ALL_BAC_dN_dS_omega_violinplots.pdf")
p +
     scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
dev.off()







print("testing just sinoC")
head(sinoC_dat)

mybxp = function(x) {
  bxp = boxplot.stats(x)[["stats"]]
  names(bxp) = c("ymin","lower", "middle","upper","ymax")
  return(bxp)
}  

# Function to use boxplot.stats for the outliers
myout = function(x) {
  data.frame(y=boxplot.stats(x)[["out"]])
}

# filter dataframe to get data to be highligheted
highlight_df <- sinoC_dat %>% 
             filter(Outlier==1)
print("HILIGHTED DATA")
head(highlight_df)

set.seed(1738);
p<-(ggplot(sinoC_dat, aes(x=class, y=value, fill=class, color = class)) 
  + geom_violin()
  + facet_wrap(~bac, labeller=label_parsed)
  #omega = 1 referene line
  + geom_hline(yintercept= 1, linetype = "dashed",color = "#939290") 
  + xlab("") 
  + ylab("Value")
  #proper colours for strip part of plot
  + scale_color_manual(values=cols_l,labels = c(" dN", " dS", expression(omega)))
)

pdf("sinoC_violinplots.pdf")
p +
     scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
dev.off()


