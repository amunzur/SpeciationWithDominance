#this script compares F2 and F1 fitness under dominance and no dominance conditions to see if dominance affects hybrid breakdown 

#load packages 
library(tidyverse)

#load the data 
fit <- read_csv("data/fitness_n[2]_N[1000]_alpha[0.025, 0.045, 0.065, 0.085, 0.105]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

fittrial <- fit %>% 
  select(rep, angle, alpha, h, F1, F2) %>% 
  mutate(hy_breakdown = F2 - F1) %>% 
  group_by(angle, alpha, h) %>% 
  summarise(hy_breakdown = mean(hy_breakdown)) %>% 
  spread(key = h, value = hy_breakdown)

fittrial2 <- fit %>% 
  select(rep, angle, alpha, h, F1, F2) %>% 
  mutate(hy_breakdown = F2 - F1) %>% 
  group_by(angle, alpha, h) %>% 
  summarise(hy_breakdown = mean(hy_breakdown))

names(fittrial)[names(fittrial) == "0.5"] <- "no_dom"
names(fittrial)[names(fittrial) == "9"] <- "dom"

fittrial$hy_breakdown_difference <- fittrial$dom - fittrial$no_dom


# make 2 different variables for different h values in alpha == 0.045
fittrial0.045 <- fittrial %>% 
  filter(alpha == 0.045)

fittrial0.045_0.5 <- fittrial %>% 
  filter(alpha == 0.045, h == 0.5)

fittrial0.045_9 <- fittrial %>% 
  filter(alpha == 0.045, h == 9)

# make 2 variables
fittrial0.065 <- fittrial %>% 
  filter(alpha == 0.065)

fittrial0.065_0.5 <- fittrial %>% 
  filter(alpha == 0.065, h == 0.5)

fittrial0.065_9 <- fittrial %>% 
  filter(alpha == 0.065, h == 9 )

# make 2 variables
fittrial0.085 <- fittrial %>% 
  filter(alpha == 0.085)

fittrial0.085_0.5 <- fittrial %>% 
  filter(alpha == 0.085, h == 0.5)

fittrial0.085_9 <- fittrial %>% 
  filter(alpha == 0.085, h == 9)

# make 2 variables
fittrial0.105 <- fittrial %>% 
  filter(alpha == 0.105)

fittrial0.105_0.5 <- fittrial %>% 
  filter(alpha == 0.105, h == 0.5)

fittrial0.105_9 <- fittrial %>% 
  filter(alpha == 0.105, h == 9)

# following plots show hybrid breakdown in separate alpha values 
# PLOT THE DATA for alpha == 0.045 ####

alpha_0.045_hy_breakdown <- ggplot(data = fittrial0.045, aes(x = angle, y = hy_breakdown)) + 
  geom_point(data = fittrial0.045_0.5, color = "blue") +
  geom_point(data = fittrial0.045_9, color = "orange") +
  geom_smooth(data = fittrial0.045_0.5, method = "lm", color = "blue", se = FALSE, size = 2) +
  geom_smooth(data = fittrial0.045_9, method = "lm", color = "orange", se = FALSE, size = 2) +
  labs(title = "Hybrid breakdown observed when alpha = 0.045", x = "angle", y = "hybrid breakdown") +
  theme_KT

alpha_0.045_hy_breakdown

# PLOT THE DATA for alpha == 0.065 ####

alpha_0.065_hy_breakdown <- ggplot(data = fittrial0.065, aes(x = angle, y = hy_breakdown)) + 
  geom_point(data = fittrial0.065_0.5, color = "blue") +
  geom_point(data = fittrial0.065_9, color = "orange") +
  geom_smooth(data = fittrial0.065_0.5, method = "lm", color = "blue", se = FALSE, size = 2) +
  geom_smooth(data = fittrial0.065_9, method = "lm", color = "orange", se = FALSE, size = 2) +
  labs(title = "Hybrid breakdown observed when alpha = 0.065", x = "angle", y = "hybrid breakdown") +
  theme_KT

alpha_0.065_hy_breakdown

# PLOT THE DATA for alpha == 0.085 ####

alpha_0.085_hy_breakdown <- ggplot(data = fittrial0.085, aes(x = angle, y = hy_breakdown)) + 
  geom_point(data = fittrial0.085_0.5, color = "blue") +
  geom_point(data = fittrial0.085_9, color = "orange") +
  geom_smooth(data = fittrial0.085_0.5, method = "lm", color = "blue", se = FALSE, size = 2) +
  geom_smooth(data = fittrial0.085_9, method = "lm", color = "orange", se = FALSE, size = 2) +
  labs(title = "Hybrid breakdown observed when alpha = 0.085", x = "angle", y = "hybrid breakdown") +
  theme_KT

alpha_0.085_hy_breakdown

# PLOT THE DATA for alpha == 0.105 ####

alpha_0.105_hy_breakdown <- ggplot(data = fittrial0.105, aes(x = angle, y = hy_breakdown)) + 
  geom_point(data = fittrial0.105_0.5, color = "blue") +
  geom_point(data = fittrial0.105_9, color = "orange") +
  geom_smooth(data = fittrial0.105_0.5, method = "lm", color = "blue", se = FALSE, size = 2) +
  geom_smooth(data = fittrial0.105_9, method = "lm", color = "orange", se = FALSE, size = 2) +
  labs(title = "Hybrid breakdown observed when alpha = 0.105", x = "angle", y = "hybrid breakdown") +
  theme_KT

alpha_0.105_hy_breakdown

fittrial0.105_new <- fittrial0.105 %>% 
  mutate(h_factor = as.factor(h))

alpha_0.105_hy_breakdown_new <- ggplot(data = fittrial0.105_new, aes(x = angle, y = -hy_breakdown, colour = h_factor)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  # scale_colour_manual()
  labs(title = "Hybrid breakdown observed when alpha = 0.105", x = "angle", y = "hybrid breakdown") +
  theme_KT

alpha_0.105_hy_breakdown_new

# following data is used to graph the hy_breakdown change for different alpha values ####

alpha0.045 <- fittrial %>% 
  filter(alpha == 0.045)

alpha0.065 <- fittrial %>% 
  filter(alpha == 0.065)

alpha0.085 <- fittrial %>% 
  filter(alpha == 0.085)

alpha0.105 <- fittrial %>% 
  filter(alpha == 0.105)

# PLOT ####
hy_breakdown_difference <- ggplot(data = fittrial, aes(x = angle, y = hy_breakdown_difference)) + 
  geom_point(data = alpha0.045, color = "black") +
  geom_point(data = alpha0.065, color = "dodgerblue3") +
  geom_point(data = alpha0.085, color = "deepskyblue1") +
  geom_point(data = alpha0.105, color = "cadetblue1") +
  geom_smooth(data = alpha0.045, method = "lm", color = "black", se = FALSE, size = 3) +
  geom_smooth(data = alpha0.065, method = "lm", color = "dodgerblue3", se = FALSE, size = 3) +
  geom_smooth(data = alpha0.085, method = "lm", color = "deepskyblue1", se = FALSE, size = 3) +
  geom_smooth(data = alpha0.105, method = "lm", color = "cadetblue1", se = FALSE, size = 3) +
  theme_KT

hy_breakdown_difference

# DELETE LATER: 
fittrial0.045 <- fittrial %>% 
  filter(alpha == 0.045)

fittrial0.045_0.5 <- fittrial %>% 
  filter(alpha == 0.045, h == 0.5)

fittrial0.045_9 <- fittrial %>% 
  filter(alpha == 0.045, h == 9)