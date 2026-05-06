
library(tidyverse)
library(plotgardener)


# Read in statistics and name columns; this is the outlier masked file that has been sorted using sort -k1,1 -k2,2n
recomb <- read.delim("modestum.10kb.sorted.recomb.txt", header = FALSE)

colnames(recomb) <- c("chromosome", "start", "end", "rate")

# Make columns numeric
recomb$rate <- as.numeric(recomb$rate)
recomb$start <- as.numeric(recomb$start)
recomb$end <- as.numeric(recomb$end)


### Comparison with irwin windows
# Read in candidate windows

north_dwgf <- read.delim("../pixy/Modestum_divergencewgeneflow_pop1windows.bed", header = FALSE)
north_sia <- read.delim("../pixy/Modestum_selectioninallopatry_pop1windows.bed", header = FALSE)

colnames(north_dwgf) <- colnames(north_sia) <- c("chromosome", "start", "end")

south_dwgf <- read.delim("../pixy/Modestum_divergencewgeneflow_pop2windows.bed", header = FALSE)
south_sia <- read.delim("../pixy/Modestum_selectioninallopatry_pop2windows.bed", header = FALSE)

colnames(south_dwgf) <- colnames(south_sia) <- c("chromosome", "start", "end")


# Add 1 to the recombination to match the bed windows from pixy
recomb$start <- recomb$start + 1

north_dwgf$type <- "dwgf"
north_sia$type <- "sia"

south_dwgf$type <- "dwgf"
south_sia$type <- "sia"


north_recomb_dwgf_comb <- left_join(recomb, north_dwgf)
north_recomb_sia_comb <- left_join(recomb, north_sia)

north_recomb_dwgf_comb$type[is.na(north_recomb_dwgf_comb$type)] <- "neutral"

south_recomb_dwgf_comb <- left_join(recomb, south_dwgf)
south_recomb_sia_comb <- left_join(recomb, south_sia)

south_recomb_dwgf_comb$type[is.na(south_recomb_dwgf_comb$type)] <- "neutral"

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

North_Pi_pop1_sims <- left_join(Pi_pop1_sims, recomb)


# Set our seed
set.seed(041426)

# Randomly select the same number of the windows as we saw empirically (North = 34, South = 27)

sims <- 10000

recomb_north_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 34)
  
  recomb_north_dwgf_sims <- c(recomb_north_dwgf_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the recombination rate for dwgf windows
north_dwgf_recomb <- north_recomb_dwgf_comb %>% filter(type == "dwgf")

north_dwgf_recomb$comp <- "north_dwgf_emp"

recomb_north_dwgf_sims_df <- data.frame(rate = recomb_north_dwgf_sims, comp = "North_dwgf_sim")

recomb_north_dwgfvssims <- rbind(north_dwgf_recomb[,c(4,6)], recomb_north_dwgf_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 13.5, df = 1, p-value = 0.0002386
kruskal.test(recomb_north_dwgfvssims$rate, recomb_north_dwgfvssims$comp)

# 1.09e-6
mean(north_dwgf_recomb$rate, na.rm = TRUE)

# 1.21e-6
mean(recomb_north_dwgf_sims_df$rate, na.rm = TRUE)

# Make a boxplot
N_dwgf_recombboxplot <- ggplot(data =recomb_north_dwgfvssims) + geom_boxplot(aes(x = comp, y = rate, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c', '#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c', '#ababab')) +
  theme_classic() + ylab("rate") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,3e-6))

## North, selection in allopatry 

sims <- 10000

recomb_north_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 85)
  
  recomb_north_sia_sims <- c(recomb_north_sia_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for sia windows
north_sia_recomb <- north_recomb_sia_comb %>% filter(type == "sia")

north_sia_recomb$comp <- "north_sia_emp"

recomb_north_sia_sims_df <- data.frame(rate = recomb_north_sia_sims, comp = "North_sia_sim")

recomb_north_siavssims <- rbind(north_sia_recomb[,c(4,6)], recomb_north_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 2.8344, df = 1, p-value = 0.2668
kruskal.test(recomb_north_siavssims$rate, recomb_north_siavssims$comp)

# 1.30e-6
mean(north_sia_recomb$rate, na.rm = TRUE)

# 1.21e-6
mean(recomb_north_sia_sims_df$rate, na.rm = TRUE)

# Make a boxplot
N_sia_recombboxplot <- ggplot(data =recomb_north_siavssims) + geom_boxplot(aes(x = comp, y = rate, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#6bacd6', '#ababab'), 0.5)) + scale_color_manual(values = c('#6bacd6', '#ababab')) +
  theme_classic() + ylab("rate") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,3e-6))

## South, divergence with gene flow

south_Pi_pop2_sims <- left_join(Pi_pop2_sims, recomb)


# Set our seed
set.seed(041426)

# Randomly select the same number of the windows as we saw empirically (south = 34, South = 27)

sims <- 10000

recomb_south_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- south_Pi_pop2_sims %>% slice_sample(n = 34)
  
  recomb_south_dwgf_sims <- c(recomb_south_dwgf_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the recombination rate for dwgf windows
south_dwgf_recomb <- south_recomb_dwgf_comb %>% filter(type == "dwgf")

south_dwgf_recomb$comp <- "south_dwgf_emp"

recomb_south_dwgf_sims_df <- data.frame(rate = recomb_south_dwgf_sims, comp = "south_dwgf_sim")

recomb_south_dwgfvssims <- rbind(south_dwgf_recomb[,c(4,6)], recomb_south_dwgf_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 22.022, df = 1, p-value = 2.695e-6
kruskal.test(recomb_south_dwgfvssims$rate, recomb_south_dwgfvssims$comp)

# 9.02e-7
mean(south_dwgf_recomb$rate, na.rm = TRUE)

# 1.19e-6
mean(recomb_south_dwgf_sims_df$rate, na.rm = TRUE)

# Make a boxplot
S_dwgf_recombboxplot <- ggplot(data =recomb_south_dwgfvssims) + geom_boxplot(aes(x = comp, y = rate, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c', '#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c', '#ababab')) +
  theme_classic() + ylab("rate") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,3e-6))

## south, selection in allopatry 

sims <- 10000

recomb_south_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- south_Pi_pop2_sims %>% slice_sample(n = 85)
  
  recomb_south_sia_sims <- c(recomb_south_sia_sims, mean(sim_dat[,15]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

# Extract the mu for sia windows
south_sia_recomb <- south_recomb_sia_comb %>% filter(type == "sia")

south_sia_recomb$comp <- "south_sia_emp"

recomb_south_sia_sims_df <- data.frame(rate = recomb_south_sia_sims, comp = "south_sia_sim")

recomb_south_siavssims <- rbind(south_sia_recomb[,c(4,6)], recomb_south_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 0.11782, df = 1, p-value = 0.7314
kruskal.test(recomb_south_siavssims$rate, recomb_south_siavssims$comp)

# 1.30e-6
mean(south_sia_recomb$rate, na.rm = TRUE)

# 1.19e-6
mean(recomb_south_sia_sims_df$rate, na.rm = TRUE)

# Make a boxplot
S_sia_recombboxplot <- ggplot(data =recomb_south_siavssims) + geom_boxplot(aes(x = comp, y = rate, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#6bacd6', '#ababab'), 0.5)) + scale_color_manual(values = c('#6bacd6', '#ababab')) +
  theme_classic() + ylab("rate") + xlab(NULL) + theme(legend.position = "none") + coord_cartesian(ylim = c(0,3e-6))


### Plotting 

pageCreate(width = 11, height = 8.5)

plotGG(N_dwgf_recombboxplot, width = 2, height = 2, x = 0.1, y = 0.1)
plotGG(S_dwgf_recombboxplot, width = 2, height = 2, x = 2.2, y = 0.1)
plotGG(N_sia_recombboxplot, width = 2, height = 2, x = 4.3, y = 0.1)
plotGG(S_sia_recombboxplot, width = 2, height = 2, x = 6.1, y = 0.1)

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
