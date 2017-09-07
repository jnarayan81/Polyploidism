#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# 20-06-2016
# works with R 2.15.2 and ggplot 0.9.3.1
# Check ggplot2 help forums or contact Jitendra Narayan jnarayan81@gmail.com if something doesn't run
# because of updated programs/packages
#Store the data to Alien variable
#Alien <- read.csv("TESTOUT2.csv")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt";
}

Alien <- read.table( args[1], sep="\t", header=TRUE)
names(Alien);
library("ggplot2")
# Basic box plot
pdf(args[2], useDingbats=FALSE)

p10 <- ggplot(Alien, aes(x = Ploidy, y = GC_per)) + 
        geom_boxplot(aes(group = cut_width(Ploidy, 0.25), fill = factor (Ploidy), outlier.colour = "red", outlier.shape = 1,notch = TRUE)) + geom_jitter(width = 0.2)
p10

dev.off()
