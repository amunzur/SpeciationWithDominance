# this script makes figures that compare hybrid fitness between parents, F1 and F2. 

library(tidyverse)
library(dplyr)
library(ggplot2)
fit_data = read.csv("data/fitness_n[2]_N[1000]_alpha[0.025, 0.045, 0.065, 0.085, 0.105]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

# select related columns
trial <- fit_data %>% 
  select(rep, angle, h, alpha, parent1, parent2, F1, F2)

# replace all the 9s in the h column with variable
trial$h[trial$h == "9"] <- "var"

# ANGLE = 0 ####
angle_0 <- trial %>% 
  filter(angle == 0) %>% 
  group_by(rep, angle, h, alpha)

#calculate rel fit in both parents' environment
F1_parent1_relfit = angle_0$F1 / angle_0$parent1
F1_parent2_relfit = angle_0$F1 / angle_0$parent2

F2_parent1_relfit = angle_0$F2 / angle_0$parent1
F2_parent2_relfit = angle_0$F2 / angle_0$parent2

#make a new column with the greater of the two relfit values 
angle_0["F1_relfit"] = pmax(F1_parent1_relfit, F1_parent2_relfit)
angle_0["F2_relfit"] = pmax(F2_parent1_relfit, F2_parent2_relfit)

angle_0 <- summarise(angle_0, parent1 = mean(parent1), parent2 = mean(parent2) , 
                     F1 = mean(F1), F2 = mean(F2),
                     F1_relfit = mean(F1_relfit),
                     F2_relfit = mean(F2_relfit))

F1_relfit0 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F1_relfit, color = h)) + geom_point()
F1_relfit0

F2_relfit0 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F2_relfit, color = h)) + geom_point()
F2_relfit0

# ANGLE = 90 ####
angle_90 <- trial %>% 
  filter(angle == 90) %>% 
  group_by(rep, angle, h, alpha)

#calculate rel fit in both parents' environment
F1_parent1_relfit = angle_90$F1 / angle_90$parent1
F1_parent2_relfit = angle_90$F1 / angle_90$parent2

F2_parent1_relfit = angle_90$F2 / angle_90$parent1
F2_parent2_relfit = angle_90$F2 / angle_90$parent2

#make a new column with the greater of the two relfit values 
angle_90["F1_relfit"] = pmax(F1_parent1_relfit, F1_parent2_relfit)
angle_90["F2_relfit"] = pmax(F2_parent1_relfit, F2_parent2_relfit)

angle_90 <- summarise(angle_90, parent1 = mean(parent1), parent2 = mean(parent2) , 
                      F1 = mean(F1), F2 = mean(F2),
                      F1_relfit = mean(F1_relfit),
                      F2_relfit = mean(F2_relfit))

F1_relfit90 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F1_relfit, color = h)) + geom_point()
F1_relfit90

F2_relfit90 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F2_relfit, color = h)) + geom_point()
F2_relfit90

# ANGLE = 180 ####

angle_180 <- trial %>% 
  filter(angle == 180) %>% 
  group_by(rep, angle, h, alpha)

#calculate rel fit in both parents' environment
F1_parent1_relfit = angle_180$F1 / angle_180$parent1
F1_parent2_relfit = angle_180$F1 / angle_180$parent2

F2_parent1_relfit = angle_180$F2 / angle_180$parent1
F2_parent2_relfit = angle_180$F2 / angle_180$parent2

#make a new column with the greater of the two relfit values 
angle_180["F1_relfit"] = pmax(F1_parent1_relfit, F1_parent2_relfit)
angle_180["F2_relfit"] = pmax(F2_parent1_relfit, F2_parent2_relfit)

angle_180 <- summarise(angle_180, parent1 = mean(parent1), parent2 = mean(parent2) , 
                       F1 = mean(F1), F2 = mean(F2),
                       F1_relfit = mean(F1_relfit),
                       F2_relfit = mean(F2_relfit))

F1_relfit180 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F1_relfit, color = h)) + geom_point()
F1_relfit180

F2_relfit180 = ggplot(data = angle_0, mapping = aes(x = alpha, y = F2_relfit, color = h)) + geom_point()
F2_relfit180

# HYBRID PHENOTYPES ####
data_phen <- read_csv("data/phenotypes_n[2]_N[1000]_alpha[0.025, 0.045, 0.065, 0.085, 0.105]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

data_plot <- data_phen %>% 
  filter(angle == 180) %>% 
  select(X1, F1_hybrid_phenos1:F2_hybrid_phenos2) %>% 
  gather(F1_hybrid_phenos1:F2_hybrid_phenos2, key = "gen_trait", value = "phenotype") %>% 
  separate(gen_trait, into = c("gen", "hybrid", "trait"), sep = "_") %>% 
  spread(key = trait, value = phenotype)

ggplot(data_plot, aes(x = phenos1, y = phenos2, color = gen)) +
  geom_point(alpha = 0.02)

# looking at hybrid breakdown
fit_mean = read_csv("data/fitMean_n[2]_N[1000]_alpha[0.01, 0.025, 0.075, 0.1]_u[0.001]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

fit_trial <- fit_mean %>% 
  filter(alpha_adapt == 0.1) %>% 
  mutate(hybrid_breakdown = F2_hybrids/F1_hybrids)

hybrid_breakdown_boxplot = ggplot(fit_trial, aes(x = factor(angle), y = hybrid_breakdown, color = factor(h))) +
  geom_boxplot()
hybrid_breakdown_boxplot
