# this script analyzes how mean h fixed changes under various alpha values. 

# load packages
library(tidyverse)

data <- read_csv("data/mutsummary_n[2]_N[1000]_alpha[0.025, 0.045, 0.065, 0.085, 0.105]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

trial <- data %>% 
  select(angle, reps, alpha, h_fixed_mean) %>% 
  group_by(as.factor(alpha)) %>% 
  na.omit() %>% 
  summarise(h_fixed = mean(h_fixed_mean), 
            se_h_fixed = stats::sd(h_fixed_mean) / sqrt(n()))

# use these data to make a supplementary figure showing that mean h fixed does not depend on alpha
trial


trial0.045 <- trial %>% 
  filter(alpha == 0.045)

trial0.065 <- trial %>% 
  filter(alpha == 0.065)

trial0.085 <- trial %>% 
  filter(alpha == 0.085)

trial0.105 <- trial %>% 
  filter(alpha == 0.105)

# GRAPH THEM ON THE SAME AXIS ####
alpha_h <- ggplot(data = trial, aes(x = angle, y = h_fixed_mean)) +
  geom_point(data = trial0.045, color = "black") +
  geom_point(data = trial0.065, color = "dodgerblue3") +
  geom_point(data = trial0.085, color = "deepskyblue1") +
  geom_point(data = trial0.105, color = "cadetblue1") +
  geom_smooth(data = trial0.045, method = "lm", color = "black", se = FALSE, size = 3) +
  geom_smooth(data = trial0.065, method = "lm", color = "dodgerblue3", se = FALSE, size = 3) +
  geom_smooth(data = trial0.085, method = "lm", color = "deepskyblue1", se = FALSE, size = 3) +
  geom_smooth(data = trial0.105, method = "lm", color = "cadetblue1", se = FALSE, size = 3) +
  labs(x = "angle", y = "mean of fixed h") +
  theme_KT

alpha_h


