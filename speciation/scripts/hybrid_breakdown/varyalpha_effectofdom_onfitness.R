#this script is used to make a figure to show how different alpha values under dom and no dom conditions affect the F1 fitness

# load packages
install.packages("viridis")
library(viridis)


HB <- read_csv("data/fitness_n[2]_N[1000]_alpha[0.025, 0.045, 0.065, 0.085, 0.105]_sigma[1]_opt_dist[1]_dom[9, 0.5].csv")

theme_KT <-
  
  theme(
    aspect.ratio = 1.0,
    
    panel.background = element_blank(),
    
    panel.grid.major = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_blank(),
    
    axis.line = element_line(size = 1),
    
    axis.line.x = element_line(color = "black", size = 1),
    
    axis.line.y = element_line(color = "black", size = 1),
    
    axis.ticks = element_line(color = "black"),
    
    axis.text = element_text(color = "black"),
    
    axis.title = element_text(color = "black"),
    
    axis.title.y = element_text(vjust = 0.2, size = 18),
    
    axis.title.x = element_text(vjust = 0.1, size = 18),
    
    axis.text.x = element_text(size = 15),
    
    axis.text.y = element_text(size = 15),
    
    legend.position = "left",
    
    legend.justification = c("right", "top")
  )

diff_dataset <- HB %>% 
  select(rep, alpha, angle, h, alpha, parent1, parent2, F1)

F1_parent1_relfit = diff_dataset$F1 / diff_dataset$parent1
F1_parent2_relfit = diff_dataset$F1 / diff_dataset$parent2

diff_dataset["F1_relfit"] = pmax(F1_parent1_relfit, F1_parent2_relfit)

diff_dataset$h[diff_dataset$h == "9"] <- "var"
diff_dataset$h[diff_dataset$h == "0.5"] <- "no_dom"

diff_general <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  group_by(angle, h, alpha) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

# DEFINE SUBSETS OF DATA FOR EACH ALPHA VALUE ####
diff_dataset_alpha0.025 <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  filter(alpha == 0.045) %>% 
  group_by(angle, h) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

diff_dataset_alpha0.045 <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  filter(alpha == 0.045) %>% 
  group_by(angle, h) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

diff_dataset_alpha0.065 <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  filter(alpha == 0.065) %>% 
  group_by(angle, h) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

diff_dataset_alpha0.085 <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  filter(alpha == 0.085) %>% 
  group_by(angle, h) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

diff_dataset_alpha0.105 <- diff_dataset %>% 
  select(rep, h, angle, alpha, F1_relfit) %>% 
  filter(alpha == 0.105) %>% 
  group_by(angle, h) %>% 
  summarise(F1_relfit = mean(F1_relfit)) %>% 
  spread(key = h, value = F1_relfit) %>% 
  mutate(dom_effect_on_fitness = var - no_dom)

# GRAPH THEM ON THE SAME AXIS ####
ggplot(data = diff_general, aes(x = angle, y = dom_effect_on_fitness)) +
  geom_point(data = diff_dataset_alpha0.045, color = "deeppink") +
  geom_point(data = diff_dataset_alpha0.065, color = "darkturquoise") +
  geom_point(data = diff_dataset_alpha0.085, color = "darkolivegreen3") +
  geom_point(data = diff_dataset_alpha0.105, color = "sienna1") +
  geom_smooth(data = diff_dataset_alpha0.045, method = "lm", color = "deeppink", se = FALSE, size = 3) +
  geom_smooth(data = diff_dataset_alpha0.065, method = "lm", color = "darkturquoise", se = FALSE, size = 3) +
  geom_smooth(data = diff_dataset_alpha0.085, method = "lm", color = "darkolivegreen3", se = FALSE, size = 3) +
  geom_smooth(data = diff_dataset_alpha0.105, method = "lm", color = "sienna1", se = FALSE, size = 3) +
  geom_hline(yintercept = 0) +
  labs(x = "angle", y = "Effect of dominance on fitness") +
  theme_KT
