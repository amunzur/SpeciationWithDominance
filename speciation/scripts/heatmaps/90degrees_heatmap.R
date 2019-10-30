# practice plotting

# load packages ####
install.packages("styler")
library(styler)

install.packages("tidyverse")
library(tidyverse)

install.packages('gplots')
library('gplots')

install.packages('RColorBrewer')
library('RColorBrewer')

degrees90 <- read.csv(file = 'data/degrees90.csv')
degrees90

#store the row names, will add later to the matrix 
rownames <- degrees90[,1]
rownames

# matrix_degrees90 <- degrees90.matrix(degrees90[,2:ncol(degrees90)])
# convert the data table into a matrix
matrix_degrees90 <- as.matrix(degrees90[,2:ncol(degrees90)])
matrix_degrees90

# add the row names to the matrix
row.names(matrix_degrees90) <- rownames

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

heatmap.2(matrix_degrees90, main = "90 degrees, u and sigma", notecol = "black", density.info = "none", trace = "none",  margins = c(12, 9),  col = my_palette,  dendrogram = "row",  Colv = "NA")

