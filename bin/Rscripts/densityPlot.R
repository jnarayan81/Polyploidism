#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

A1 <- read.table( args[1], sep="\t", header=TRUE)
pdf(args[2], useDingbats=FALSE)

df <- data.frame(A1)
df
library(ggplot2)
# Basic density

hist(df$Ploidy, 
     main="Histogram for Ploidy", 
     xlab="Ploidy", 
     border="blue", 
     col="green",
     las=1, 
     breaks=5)

# Add mean line
p<-ggplot(df, aes(x=df$Ploidy)) + 
  geom_histogram(color="black", fill="white")

p+ geom_vline(aes(xintercept=mean(df$Ploidy)),
            color="blue", linetype="dashed", size=1)
# Histogram with density plot
ggplot(df, aes(x=df$Ploidy)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = NULL)+
 geom_density(alpha=.2, fill="#FF6666") 
