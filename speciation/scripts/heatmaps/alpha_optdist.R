library(tidyverse)

library('gplots')

#add the colors
my_palette <- colorRampPalette(c("green", "white", "red"))(n = 299)

aug_30 = read.csv("data/summary_n[2]_N[700]_alpha[0.01, 0.025, 0.05, 0.075, 0.1, 0.2]_u[0.001]_sigma[1]_opt_dist[1, 2, 3, 4, 5]_dom['variable'].csv")

practice

practice <- aug_30 %>% 
  group_by(alpha_adapt, opt_dist.1) %>%   
  select(h_fixed_mean) %>% 
  summarise(h_fixed_mean = mean(h_fixed_mean))

#convert the csv into a matrix 
mat_prac <- as.matrix(practice[,1:ncol(practice)])
mat_prac

#this didnt work, i will manually add the data in 

#load the data in the correct format for the heatmap
map_data <- read.csv("opt_dist_and_alpha.csv")

#save the row names
rownames <- map_data[,1]
rownames

#convert the csv into a matrix 
mat_data <- as.matrix(map_data[,2:ncol(map_data)])
mat_data

# add the row names to the matrix
row.names(mat_data) <- rownames

heatmap.2(mat_data, main = 'alpha_optdist', notecol = "black", density.info = "none", trace = "none", margins = c(12, 9),  col = my_palette,  dendrogram = "none",  Rowv = FALSE, Colv = FALSE,)
