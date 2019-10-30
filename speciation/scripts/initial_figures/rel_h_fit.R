# this is a script to plot the hybrid fitness from simulations with various sigma values - sept 8th 

# load packages ####
library(styler)

library(tidyverse)

library('gplots')

library('RColorBrewer')

library('circlize')

library('ggplot2')

# load the data 
fit <- read.csv("data/hybrid_fit/HF in various sigma /fitness_n[2]_N[1000]_alpha[0.1]_u[0.001]_sigma[1, 3, 5]_opt_dist[1]_dom[9, 0.5].csv")

# find the relative hybrid fitness compared to both parents 
fit <- transform(fit, h_parent1 = hybrids / parent1)
fit <- transform(fit, h_parent2 = hybrids / parent2)

# find which relative fitness is greater and make it a new column 
fit <- transform(fit, rel_fit = pmax(fit$h_parent1, fit$h_parent2))

# lets work on the case where sigma = 1 
# select the columns we need 
practice <- fit %>% 
  select(angle, sigma, rel_fit) %>% 
  filter (sigma == 1) %>% 
  arrange(angle)

# plot the data 
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point()

lw1 <- loess(rel_fit ~ angle,data=practice)            
plot(rel_fit ~ angle, data=practice,pch=19,cex=0.1)               
j <- order(practice$angle)              
lines(practice$angle[j],lw1$fitted[j])  


# sigma = 3 
practice <- fit %>% 
  select(angle, sigma, rel_fit) %>% 
  filter (sigma == 3) %>% 
  arrange(angle)

# plot the data 
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point()

# sigma = 5 
practice <- fit %>% 
  select(angle, sigma, rel_fit) %>% 
  filter (sigma == 3) %>% 
  arrange(angle)

# plot the data 
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point()

