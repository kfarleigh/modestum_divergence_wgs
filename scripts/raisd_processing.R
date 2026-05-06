

library(tidyverse)
library(plotgardener)


# Read in statistics and name columns
north_raisd <- read.delim("north.10kb.raisd.txt", header = FALSE)
south_raisd <- read.delim("south.10kb.raisd.txt", header = FALSE)


colnames(north_raisd) <- colnames(south_raisd) <- c("chromosome", "start", "end", "mu")


# Make columns numeric
north_raisd$start <- as.numeric(north_raisd$start)
north_raisd$end <- as.numeric(north_raisd$end)
north_raisd$mu <- as.numeric(north_raisd$mu)

south_raisd$start <- as.numeric(south_raisd$start)
south_raisd$end <- as.numeric(south_raisd$end)
south_raisd$mu <- as.numeric(south_raisd$mu)

# How many windows do we have on each chromosome
table(north_raisd$chromosome)
table(south_raisd$chromosome)

# How many NAs
summary(north_raisd$mu)
summary(south_raisd$mu)


# Let's make a quick histogram of the data for each populatoin
hist(north_raisd$mu)
hist(south_raisd$mu)

### Calculate z-scores to identify outlier windows

north_raisd$z <- (north_raisd$mu - mean(north_raisd$mu, na.rm = TRUE)) / sd(north_raisd$mu, na.rm = TRUE)
south_raisd$z <- (south_raisd$mu - mean(south_raisd$mu, na.rm = TRUE)) / sd(south_raisd$mu, na.rm = TRUE)

# Create a column to identify candidates, 0 is no, 1 is yes
north_raisd$cand <- 0
north_raisd$cand[which(north_raisd$z >= 3)] <- 1 

# Which windows were outliers (z-score >=3)
north_raisd_candidates <- north_raisd[which(north_raisd$z >= 3),]
south_raisd_candidates <- south_raisd[which(south_raisd$z >= 3),]

# What proportion of windows are candidates?
nrow(north_raisd_candidates)/nrow(north_raisd)
nrow(south_raisd_candidates)/nrow(south_raisd)

### Do raisd windows overlap with irwin model windows?
# Read in candidate windows

north_dwgf <- read.delim("../pixy/Modestum_divergencewgeneflow_pop1windows.bed", header = FALSE)
north_sia <- read.delim("../pixy/Modestum_selectioninallopatry_pop1windows.bed", header = FALSE)

colnames(north_dwgf) <- colnames(north_sia) <- c("chromosome", "start", "end")

south_dwgf <- read.delim("../pixy/Modestum_divergencewgeneflow_pop2windows.bed", header = FALSE)
south_sia <- read.delim("../pixy/Modestum_selectioninallopatry_pop2windows.bed", header = FALSE)

colnames(south_dwgf) <- colnames(south_sia) <- c("chromosome", "start", "end")


# Add 1 to the raisd results to match the bed windows from pixy
north_raisd_candidates$start <- north_raisd_candidates$start + 1

south_raisd_candidates$start <- south_raisd_candidates$start + 1

# Create metadata columns to compare the two methods
north_dwgf$meta <- paste(north_dwgf$chromosome, north_dwgf$start, north_dwgf$end, sep = "_")
north_sia$meta <- paste(north_sia$chromosome, north_sia$start, north_sia$end, sep = "_")
north_raisd_candidates$meta <- paste(north_raisd_candidates$chromosome, north_raisd_candidates$start, north_raisd_candidates$end, sep = "_")

north_dwgf$type <- "dwgf"
north_sia$type <- "sia"

table(north_dwgf$meta %in% north_raisd_candidates$meta)
table(north_sia$meta %in% north_raisd_candidates$meta)

south_dwgf$meta <- paste(south_dwgf$chromosome, south_dwgf$start, south_dwgf$end, sep = "_")
south_sia$meta <- paste(south_sia$chromosome, south_sia$start, south_sia$end, sep = "_")
south_raisd_candidates$meta <- paste(south_raisd_candidates$chromosome, south_raisd_candidates$start, south_raisd_candidates$end, sep = "_")

south_dwgf$type <- "dwgf"
south_sia$type <- "sia"

table(south_dwgf$meta %in% south_raisd_candidates$meta)
table(south_sia$meta %in% south_raisd_candidates$meta)

# Or we can just bind them
north_raisd$start <- north_raisd$start + 1

north_raisd_dwgf_comb <- left_join(north_raisd, north_dwgf)
north_raisd_sia_comb <- left_join(north_raisd, north_sia)

north_raisd_dwgf_comb$type[is.na(north_raisd_dwgf_comb$type)] <- "neutral"

south_raisd$start <- south_raisd$start + 1

south_raisd_dwgf_comb <- left_join(south_raisd, south_dwgf)
south_raisd_sia_comb <- left_join(south_raisd, south_sia)

south_raisd_dwgf_comb$type[is.na(south_raisd_dwgf_comb$type)] <- "neutral"

### Is the distribution of mu higher in dwgf windows than a genomic background?
# We will do this the same way that we did the Tajima's D simulations 

## Read in pi data 

Pi_dat <- read.csv("/Users/keaka/Downloads/Modestum_pixy_pi_sign.csv")

colnames(Pi_dat)[4:5] <- c("start", "end")

# Split by pop
Pi_pop1 <- Pi_dat[which(Pi_dat$pop == "V1"),]
Pi_pop2 <- Pi_dat[which(Pi_dat$pop == "V2"),] 

# Calculate z-score 
Pi_pop1$z_score_bypop <- (Pi_pop1$avg_pi-mean(Pi_pop1$avg_pi))/sd(Pi_pop1$avg_pi)
Pi_pop2$z_score_bypop <- (Pi_pop2$avg_pi-mean(Pi_pop2$avg_pi))/sd(Pi_pop2$avg_pi)

## Start with the north population
# Get pi windows with z-scores in the divergence with gene flow range
Pi_pop1_sims <- Pi_pop1 %>% filter(z_score_bypop < -0.67 & z_score_bypop > -3)
Pi_pop2_sims <- Pi_pop2 %>% filter(z_score_bypop < -0.67 & z_score_bypop > -3)

North_Pi_pop1_sims <- left_join(Pi_pop1_sims, north_raisd)


# Set our seed
set.seed(041426)

north_raisd_candidates$comp <- "north_raisd_emp"
# Randomly select the same number of the windows as we saw empirically (North = 34, South = 27)

sims <- 10000

mu_north_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 34)
  
  mu_north_dwgf_sims <- c(mu_north_dwgf_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for dwgf windows
north_dwgf_mu <- north_raisd_dwgf_comb %>% filter(type == "dwgf")

north_dwgf_mu$comp <- "north_dwgf_emp"

mu_north_dwgf_sims_df <- data.frame(mu = mu_north_dwgf_sims, comp = "North_dwgf_sim")

mu_north_dwgfvssims <- rbind(north_dwgf_mu[,c(4,9)], mu_north_dwgf_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 2.8344, df = 1, p-value = 0.092
kruskal.test(mu_north_dwgfvssims$mu, mu_north_dwgfvssims$comp)

# 4.804842
mean(north_dwgf_mu$mu, na.rm = TRUE)

# 4.012542
mean(mu_north_dwgf_sims_df$mu, na.rm = TRUE)

# Make a boxplot
N_dwgf_raisdboxplot <- ggplot(data =mu_north_dwgfvssims) + geom_boxplot(aes(x = comp, y = mu, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c', '#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c', '#ababab')) +
  theme_classic() + ylab("mu") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,8))

## North, selection in allopatry 

sims <- 10000

mu_north_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 85)
  
  mu_north_sia_sims <- c(mu_north_sia_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for sia windows
north_sia_mu <- north_raisd_sia_comb %>% filter(type == "sia")

north_sia_mu$comp <- "north_sia_emp"

mu_north_sia_sims_df <- data.frame(mu = mu_north_sia_sims, comp = "North_sia_sim")

mu_north_siavssims <- rbind(north_sia_mu[,c(4,9)], mu_north_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 2.8344, df = 1, p-value < 2.2 e-16
kruskal.test(mu_north_siavssims$mu, mu_north_siavssims$comp)

# 5.998013
mean(north_sia_mu$mu, na.rm = TRUE)

# 4.015004
mean(mu_north_sia_sims_df$mu, na.rm = TRUE)

# Make a boxplot
N_sia_raisdboxplot <- ggplot(data =mu_north_siavssims) + geom_boxplot(aes(x = comp, y = mu, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#6bacd6', '#ababab'), 0.5)) + scale_color_manual(values = c('#6bacd6', '#ababab')) +
  theme_classic() + ylab("mu") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,10))

## South, divergence with gene flow

south_Pi_pop2_sims <- left_join(Pi_pop2_sims, south_raisd)


# Set our seed
set.seed(041426)

south_raisd_candidates$comp <- "south_raisd_emp"
# Randomly select the same number of the windows as we saw empirically (south = 34, South = 27)

sims <- 10000

mu_south_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- south_Pi_pop2_sims %>% slice_sample(n = 27)
  
  mu_south_dwgf_sims <- c(mu_south_dwgf_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for dwgf windows
south_dwgf_mu <- south_raisd_dwgf_comb %>% filter(type == "dwgf")

south_dwgf_mu$comp <- "south_dwgf_emp"

mu_south_dwgf_sims_df <- data.frame(mu = mu_south_dwgf_sims, comp = "south_dwgf_sim")

mu_south_dwgfvssims <- rbind(south_dwgf_mu[,c(4,8)], mu_south_dwgf_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 2.8344, df = 1, p-value = 0.001
kruskal.test(mu_south_dwgfvssims$mu, mu_south_dwgfvssims$comp)

# 2.444266
mean(south_dwgf_mu$mu, na.rm = TRUE)

# 1.701605
mean(mu_south_dwgf_sims_df$mu, na.rm = TRUE)

# Make a boxplot
S_dwgf_raisdboxplot <- ggplot(data =mu_south_dwgfvssims) + geom_boxplot(aes(x = comp, y = mu, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c', '#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c', '#ababab')) +
  theme_classic() + ylab("mu") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,4))


## South, selection in allopatry 

sims <- 10000

mu_south_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- south_Pi_pop2_sims %>% slice_sample(n = 90)
  
  mu_south_sia_sims <- c(mu_south_sia_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for sia windows
south_sia_mu <- south_raisd_sia_comb %>% filter(type == "sia")

south_sia_mu$comp <- "south_sia_emp"

mu_south_sia_sims_df <- data.frame(mu = mu_south_sia_sims, comp = "south_sia_sim")

mu_south_siavssims <- rbind(south_sia_mu[,c(4,8)], mu_south_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 38.147, df = 1, p-value = 6.562e-10
kruskal.test(mu_south_siavssims$mu, mu_south_siavssims$comp)

# 2.37906
mean(south_sia_mu$mu, na.rm = TRUE)

# 1.702572
mean(mu_south_sia_sims_df$mu, na.rm = TRUE)

# Make a boxplot
S_sia_raisdboxplot <- ggplot(data =mu_south_siavssims) + geom_boxplot(aes(x = comp, y = mu, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#6bacd6', '#ababab'), 0.5)) + scale_color_manual(values = c('#6bacd6', '#ababab')) +
  theme_classic() + ylab("mu") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,5))


### Plotting 

pageCreate(width = 11, height = 8.5)

plotGG(N_dwgf_raisdboxplot, width = 2, height = 2, x = 0.1, y = 0.1)
plotGG(S_dwgf_raisdboxplot, width = 2, height = 2, x = 2.2, y = 0.1)
plotGG(N_sia_raisdboxplot, width = 2, height = 2, x = 4.3, y = 0.1)
plotGG(S_sia_raisdboxplot, width = 2, height = 2, x = 6.1, y = 0.1)

pageGuideHide()


# Plot genome scans

Chrom_names <- unique(north_raisd$chromosome)

Chrom_names <- Chrom_names[c(1,3,17,18,6,8,10,2,4,5,7,9,11:16)]

north_raisd$chromosome <- factor(north_raisd$chromosome, levels = Chrom_names)
south_raisd$chromosome <- factor(south_raisd$chromosome, levels = Chrom_names)

north_raisd$Plottin.name <- north_raisd$chromosome

north_raisd %>%
  #mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
  #TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:18))) %>%
  ggplot(aes(x = (start + end)/2, y = mu, color = "#ebebeb"))+
  geom_line(linewidth = 1.25, alpha = 0.25)+
  facet_grid(mu  ~ Plottin.name,
             scales = "free", switch = "x", space = "free_x") +
  xlab("Chromsome")+
  ylab("Statistic Value")+
  #scale_color_manual(values = c("#9e0142", "#5e4fa2"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none", 
        strip.text.x = element_text(angle = 90))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))



