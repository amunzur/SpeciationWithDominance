# load packages ####
install.packages("styler")
library(styler)

install.packages("tidyverse")
library(tidyverse)

install.packages('gplots')
library('gplots')

install.packages('RColorBrewer')
library('RColorBrewer')

install.packages('circlize')
library('circlize')

library('ggplot2')

getwd()
summary <- read.csv("data/summary_averaged_h_fixed.csv")
summary

#store the row names, will add later to the matrix 
rownames <- summary[,1]
rownames

# convert the data table into a matrix
matsummary <- as.matrix(summary[,2:ncol(summary)])
matsummary

# add the row names to the matrix
row.names(matsummary) <- rownames

#add the colors
my_palette <- colorRampPalette(c("green", "white", "red"))(n = 299)

heatmap.2(matsummary, main = "u and sigma", notecol = "black", density.info = "none", trace = "none", margins = c(12, 9),  col = my_palette,  dendrogram = "none",  Rowv = FALSE, Colv = FALSE, xlab = "u", ylab = "sigma")

title(ylab="sigma", mgp=c(1,1,0), family="Calibri Light",cex.lab=1.2)
