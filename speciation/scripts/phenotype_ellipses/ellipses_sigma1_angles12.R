# this script graphs the phenotypes to make sense of the fitness graphs 

library(dplyr)
library(ggplot2)

pheno <- read.csv("data/hybrid_fit/HF in various sigma /phenotypes_n[2]_N[1000]_alpha[0.1]_u[0.001]_sigma[1, 3, 5]_opt_dist[1]_dom[9, 0.5].csv")
# SIGMA IS 1
# h is variable
pheno_sigma_1_h_var <- pheno[1:180000, ]

pheno_sigma_1_h_var $ h <- "variable"

pheno_sigma_1_h_var <- pheno_sigma_1_h_var %>%  
  select(reps, angle, hybrid_phenos1, hybrid_phenos2, h)

# h is always 0.5 
pheno_sigma_1_h_0.5 <- pheno[540001:720000, ]

pheno_sigma_1_h_0.5 $ h <- "0.5"

pheno_sigma_1_h_0.5 <- pheno_sigma_1_h_0.5 %>% 
  select(reps, angle, hybrid_phenos1, hybrid_phenos2, h)

pheno_sigma1 <- rbind(pheno_sigma_1_h_var, pheno_sigma_1_h_0.5)

# work with angle = 0 
pheno_sigma1_angle0 <- pheno_sigma1 %>% 
  filter(angle == 0, 
         reps == 1)

# separate the variable so that h values are loaded in different variables
pheno_sigma1_angle0_hvar <- pheno_sigma1_angle0[1:15000, ]
pheno_sigma1_angle0_h05 <- pheno_sigma1_angle0[15001:30000, ]

# express hvar and h05 at angle0 in the same graph
fit_fig_angle0 <- ggplot(data = pheno_sigma1_angle0, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle0") +
  geom_point(x = 1, y = 0, shape = 16, colour = "black", size = 3)
  
fit_fig_angle0

# work with angle = 16.36
pheno_sigma1_angle16.36 <- pheno_sigma1 %>% 
  filter(angle == 16.36,
         reps == 1)

# plot
fit_fig_angle16.36 <- ggplot(data = pheno_sigma1_angle16.36, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle16.36")

fit_fig_angle16.36

# work with angle = 32.73
pheno_sigma1_angle32.73 <- pheno_sigma1 %>% 
  filter(angle == 32.73)

# plot
fit_fig_angle32.73 <- ggplot(data = pheno_sigma1_angle32.73, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle32.73")

fit_fig_angle32.73

# work with angle = 49.09	
pheno_sigma1_angle49.09	 <- pheno_sigma1 %>% 
  filter(angle == 49.09	)

# plot
fit_fig_angle49.09	 <- ggplot(data = pheno_sigma1_angle49.09, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle49.09")

fit_fig_angle49.09	

# work with angle = 65.45	
pheno_sigma1_angle65.45	 <- pheno_sigma1 %>% 
  filter(angle == 65.45	)

# plot
fit_fig_angle65.45	 <- ggplot(data = pheno_sigma1_angle65.45, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle65.45")

fit_fig_angle65.45	

# work with angle = 81.82		
pheno_sigma1_angle81.82		 <- pheno_sigma1 %>% 
  filter(angle == 81.82, 
         reps == 1)

# plot
fit_fig_angle81.82 <- ggplot(data = pheno_sigma1_angle81.82, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  # coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle81.82	") +
  geom_point(x = 0, y = 0, shape = 16, colour = "black", size = 3)

fit_fig_angle81.82		

# work with angle = 98.18			
pheno_sigma1_angle98.18 <- pheno_sigma1 %>% 
  filter(angle == 98.18)

# plot
fit_fig_angle98.18 <- ggplot(data = pheno_sigma1_angle98.18, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle98.18")

fit_fig_angle98.1			

# work with angle = 114.55			
pheno_sigma1_angle114.55 <- pheno_sigma1 %>% 
  filter(angle == 114.55)

# plot
fit_fig_angle114.55 <- ggplot(data = pheno_sigma1_angle114.55, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle114.55")

fit_fig_angle114.55

# work with angle = 130.91			
pheno_sigma1_angle130.91 <- pheno_sigma1 %>% 
  filter(angle == 130.91)

# plot
fit_fig_angle130.91 <- ggplot(data = pheno_sigma1_angle130.91, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle130.91")

fit_fig_angle130.91

# work with angle = 147.27				
pheno_sigma1_angle147.27 <- pheno_sigma1 %>% 
  filter(angle == 147.27)

# plot
fit_fig_angle147.27 <- ggplot(data = pheno_sigma1_angle147.27	, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle147.27")

fit_fig_angle147.27	

# work with angle = 163.64					
pheno_sigma1_angle163.64 <- pheno_sigma1 %>% 
  filter(angle == 163.64, 
         reps == 1)

# plot
fit_fig_angle163.64	<- ggplot(data = pheno_sigma1_angle163.64	, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle163.64")

fit_fig_angle163.64	

# work with angle = 180						
pheno_sigma1_angle180 <- pheno_sigma1 %>% 
  filter(angle == 180.00, 
         reps == 1)

# plot
fit_fig_angle180	<- ggplot(data = pheno_sigma1_angle180	, mapping = aes(x = hybrid_phenos1, y = hybrid_phenos2, colour = h)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 1/2) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggtitle("sigma1_angle180") +
  geom_point(x = 0, y = 0, shape = 16, colour = "black", size = 3) +
  geom_point(x = -1, y = 0, shape = 16, colour = "black", size = 3) +
  geom_point(x = 1, y = 0, shape = 16, colour = "black", size = 3)

fit_fig_angle180	

#these scripts deal with the mean values only. 
pheno_sigma1_1 <- pheno_sigma1 %>% 
  select(reps, angle, hybrid_phenos1, hybrid_phenos2, h) %>% 
  group_by(angle, h) %>%
  summarise(mean_hy_phenos1 = mean(hybrid_phenos1), 
            mean_hy_phenos2 = mean(hybrid_phenos2))

fit_fig <- ggplot(data = pheno_sigma1_1, mapping = aes(x = mean_hy_phenos1, y = mean_hy_phenos2, colour = h)) +
  geom_point() +
  geom_smooth() +
  
fit_fig

pheno_sigma1_2 <- pheno_sigma1 %>% 
  select(reps, angle, hybrid_phenos1, hybrid_phenos2, h) %>% 
  group_by(angle, reps, h) %>%
  summarise(mean_hy_phenos1 = mean(hybrid_phenos1), 
            mean_hy_phenos2 = mean(hybrid_phenos2)) %>% 
  filter(angle == 98.18)

fit_fig_2 <- ggplot(data = pheno_sigma1_2, mapping = aes(x = mean_hy_phenos1, y = mean_hy_phenos2, colour = h)) +
  geom_point() +
  geom_smooth()
fit_fig_2