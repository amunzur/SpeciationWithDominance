# this is a script to plot the hybrid fitness from simulations with various sigma values - sept 8th

# load packages ####
library(styler)

library(tidyverse)

library("gplots")

library("RColorBrewer")

library("circlize")

library("ggplot2")

# load the data
fit <- read.csv("data/hybrid_fit/HF in various sigma /fitness_n[2]_N[1000]_alpha[0.1]_u[0.001]_sigma[1, 3, 5]_opt_dist[1]_dom[9, 0.5].csv")
# theme
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

  axis.title.y = element_text(vjust = 0.2, size = 12),

  axis.title.x = element_text(vjust = 0.1, size = 12),

  axis.text.x = element_text(size = 10),

  axis.text.y = element_text(size = 10),

  legend.position = "none"
)


# temp: add grouping var every 1k rows.
fit_groups <- fit %>%
  mutate(rep = 1:nrow(.) %/% 1e3) %>%
  group_by(rep) %>%
  filter(row_number() != 1 & row_number() != n()) %>%
  mutate(h = ifelse(rep < 541, "variable", "no_dom"))

# create summary data to plot
fit_plottable <- fit_groups %>%
  group_by(rep, h, sigma, angle) %>%
  summarise(
    mean_fit_p1 = mean(parent1),
    mean_fit_p2 = mean(parent2),
    mean_fit_hy = mean(hybrids)
  ) %>%
  mutate(rel_hy_fit = mean_fit_hy / mean(mean_fit_p1, mean_fit_p2)) %>%
  filter(sigma == 1)

fit_fig <- ggplot(data = fit_plottable, mapping = aes(x = angle, y = rel_hy_fit, colour = h)) +
  geom_point() +
  geom_smooth() +
  theme_KT
fit_fig

# next time: need rep, need h column.

# find the relative hybrid fitness compared to both parents
# fit <- transform(fit, h_parent1 = hybrids / parent1)
# fit <- transform(fit, h_parent2 = hybrids / parent2)

# find which relative fitness is greater and make it a new column
# fit <- transform(fit, rel_fit = pmax(fit$h_parent1, fit$h_parent2))

# h is variable
# fit <- fit[540001:1080000, ]

# lets work on the case where sigma = 1
# select the columns we need
names
practice <- fit %>%
  select(angle, sigma, rel_fit) %>%
  filter(sigma == 1)

# plot the data
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point() + geom_smooth()

# sigma = 3
practice <- fit %>%
  select(angle, sigma, rel_fit) %>%
  filter(sigma == 3)

# plot the data
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point() + geom_smooth()

# sigma = 5
practice <- fit %>%
  select(angle, sigma, rel_fit) %>%
  filter(sigma == 5)

# plot the data
ggplot(data = practice, mapping = aes(x = angle, y = rel_fit)) + geom_point() + geom_smooth()
