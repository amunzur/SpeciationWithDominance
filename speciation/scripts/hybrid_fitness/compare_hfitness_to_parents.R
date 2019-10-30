# this script is for comparing hybrid fitness to parent fitness with and without dominance. 

dominance <- read.csv("data/fitness_n[2]_N[1000]_alpha[0.01]_u[0.001]_sigma[1]_opt_dist[1]_dom[0.5].csv")

# find the average fitnesses when angle is 0 
practice <- dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 0) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))

# find the average fitnesses when angle is 90
practice <- dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 90) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))

# find the average fitnesses when angle is 180 
practice <- dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 180) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))

#load the no dominance data 
no_dominance <- read.csv("data/fitness_n[2]_N[1000]_alpha[0.01]_u[0.001]_sigma[1]_opt_dist[1]_dom['variable'].csv")

# find the average fitnesses when angle is 0 
practice <- no_dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 0) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))

# find the average fitnesses when angle is 90
practice <- no_dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 90) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))

# find the average fitnesses when angle is 180 
practice <- no_dominance %>% 
  select(parent1, parent2, hybrids, angle) %>% 
  filter(angle == 180) %>% 
  summarise(parent1 = mean(parent1), parent2 = mean(parent2), hybrids = mean(hybrids))
