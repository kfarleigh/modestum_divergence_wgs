### Load your packages
library(tidyverse)
library(plotgardener)
library(scales)

setwd("/Users/kfarleigh/Downloads")

### Read in files; these were created by concatinating pixy output (cat *fst* > Modestum_pixy_fst_sign.csv)
Fst_dat <- read.csv("Modestum_pixy_fst_sign.csv")
Dxy_dat <- read.csv("Modestum_pixy_dxy_sign.csv")
Pi_dat <- read.csv("Modestum_pixy_pi_sign.csv")

# Split by pop
Pi_pop1 <- Pi_dat[which(Pi_dat$pop == "V1"),]
Pi_pop2 <- Pi_dat[which(Pi_dat$pop == "V2"),] 

# Calculate z-score 
Pi_pop1$z_score_bypop <- (Pi_pop1$avg_pi-mean(Pi_pop1$avg_pi))/sd(Pi_pop1$avg_pi)
Pi_pop2$z_score_bypop <- (Pi_pop2$avg_pi-mean(Pi_pop2$avg_pi))/sd(Pi_pop2$avg_pi)

Pi_comb <- rbind(Pi_pop1, Pi_pop2)

# Average pi between populations 
Pi_pop1$avg_pi_pop1_pop2 <- Pi_pop1$avg_pi+Pi_pop2$avg_pi/2

# Isolate columns we want in each file 
Fst_filt <- Fst_dat[,c(1,2,5:7,10)]
colnames(Fst_filt)[6] <- "z_score_fst"

Dxy_filt <- Dxy_dat[,c(1,2,5:7,13)]
colnames(Dxy_filt)[6] <- "z_score_dxy"

Fst_dxy_comb <- left_join(Fst_filt, Dxy_filt)

Fst_dxy_pi_comb <- left_join(Fst_dxy_comb, Pi_pop1)

Fst_dxy_pi_comb$z_score_pi <- (Fst_dxy_pi_comb$avg_pi_pop1_pop2-mean(Fst_dxy_pi_comb$avg_pi_pop1_pop2))/sd(Fst_dxy_pi_comb$avg_pi_pop1_pop2)

# Figure out which ones satisfy the criteria
Div_gen <- Fst_dxy_pi_comb[which(Fst_dxy_pi_comb$z_score_fst > 3 & Fst_dxy_pi_comb$z_score_dxy > 3 & Fst_dxy_pi_comb$z_score_pi < -0.67 & Fst_dxy_pi_comb$z_score_pi > -3),]

SA_gen <- Fst_dxy_pi_comb[which(Fst_dxy_pi_comb$z_score_fst > 3 & Fst_dxy_pi_comb$z_score_dxy > -0.67 & Fst_dxy_pi_comb$z_score_dxy < 0.67 &  Fst_dxy_pi_comb$z_score_pi < -0.67 & Fst_dxy_pi_comb$z_score_pi > -3),]

# Make sure they match

Div_gen$meta <- paste(Div_gen$chromosome, Div_gen$window_pos_1, Div_gen$window_pos_2, sep = "_")

Div_gen_pop1_winds$meta <- paste(Div_gen_pop1_winds$chromosome, Div_gen_pop1_winds$window_pos_1, Div_gen_pop1$window_pos_2, sep = "_")
Div_gen_pop2_winds$meta <- paste(Div_gen_pop2_winds$chromosome, Div_gen_pop2_winds$window_pos_1, Div_gen_pop2$window_pos_2, sep = "_")

DG_winds <- c(Div_gen_pop1_winds$meta, Div_gen_pop2_winds$meta)
DG_winds <- unique(DG_winds)

table(Div_gen$meta %in% Div_gen_pop1_winds$meta)
table(Div_gen$meta %in% Div_gen_pop2_winds$meta)

table(Div_gen$meta %in% DG_winds)

# Find windows in each of the comparisons to make sure they have already been identified
Div_gen_avg_p1_shared <- Div_gen[which(Div_gen$meta %in% Div_gen_pop1_winds$meta),]
Div_gen_avg_p2_shared <- Div_gen[which(Div_gen$meta %in% Div_gen_pop2_winds$meta),]

DG_shared_uniq <- unique(c(Div_gen_avg_p1_shared$meta, Div_gen_avg_p2_shared$meta))

# Do it for SA
SA_gen$meta <- paste(SA_gen$chromosome, SA_gen$window_pos_1, SA_gen$window_pos_2, sep = "_")

SA_pop1_winds$meta <- paste(SA_pop1_winds$chromosome, SA_pop1_winds$window_pos_1, SA_pop1$window_pos_2, sep = "_")
SA_pop2_winds$meta <- paste(SA_pop2_winds$chromosome, SA_pop2_winds$window_pos_1, SA_pop2$window_pos_2, sep = "_")

SA_winds <- c(SA_pop1_winds$meta, SA_pop2_winds$meta)
SA_winds <- unique(SA_winds)


table(SA_gen$meta %in% SA_pop1_winds$meta)
table(SA_gen$meta %in% SA_pop2_winds$meta)

table(SA_gen$meta %in% SA_winds)

Div_gen_avg_p1_shared <- Div_gen[which(Div_gen$meta %in% Div_gen_pop1_winds$meta),]
Div_gen_avg_p2_shared <- Div_gen[which(Div_gen$meta %in% Div_gen_pop2_winds$meta),]

DG_shared_uniq <- unique(c(Div_gen_avg_p1_shared$meta, Div_gen_avg_p2_shared$meta))

######### Find the peaks and troughs in our Fst and Dxy data 

# Find the Fst peaks (z_score > 3)
Fst_dxy_comb_filt <- Fst_dxy_comb[which(Fst_dxy_comb$z_score_fst > 3),]

# Find windows with Fst and Dxy peaks (z-score >  3)
Fst_peak_dxy_peak <- Fst_dxy_comb_filt[which(Fst_dxy_comb_filt$z_score_dxy >= 3),]

# Find windows where there is an Fst peak and Dxy is about average (z-score between -0.67 and 0.67)
Fst_peak_dxy_avg <- Fst_dxy_comb_filt[which(Fst_dxy_comb_filt$z_score_dxy >= -0.67 & Fst_dxy_comb_filt$z_score_dxy <= 0.67),]

# Find windows where there is an Fst peak and Dxy is less than average (z-score between -0.67 and -3)
Fst_peak_dxy_minitrough <- Fst_dxy_comb_filt[which(Fst_dxy_comb_filt$z_score_dxy > -3 & Fst_dxy_comb_filt$z_score_dxy < -0.67),]

# Find windows where there is an Fst peak and Dxy is significantly less than average (z-score < -3)
Fst_peak_dxy_trough <- Fst_dxy_comb_filt[which(Fst_dxy_comb_filt$z_score_dxy < -3),]

#### Pair with nucleotide diversity estimates 
# Divergence with gene flow, fst and dxy peaks, avg pi

Fst_dxy_peak_pi <- left_join(Fst_peak_dxy_peak, Pi_comb)

Fst_dxy_peak_pi_avg <- Fst_dxy_peak_pi[which(Fst_dxy_peak_pi$z_score_bypop > -3 & Fst_dxy_peak_pi$z_score_bypop < -0.67),]

# Split by population
Div_gen_pop1 <- Fst_dxy_peak_pi_avg[which(Fst_dxy_peak_pi_avg$pop == "V1"),]
Div_gen_pop2 <- Fst_dxy_peak_pi_avg[which(Fst_dxy_peak_pi_avg$pop == "V2"),]

# Remove to clear space in global environment
remove(Fst_dxy_peak_pi, Fst_dxy_peak_pi_avg)

# Selection in allopatry, high fst, average dxy, low pi
Fst_peak_dxy_avg_pi <- left_join(Fst_peak_dxy_avg, Pi_comb)

Fst_peak_dxy_avg_pi_minitrough <- Fst_peak_dxy_avg_pi[which(Fst_peak_dxy_avg_pi$z_score_bypop > -3 & Fst_peak_dxy_avg_pi$z_score_bypop < -0.67),]

# Split by population 
SA_pop1 <- Fst_peak_dxy_avg_pi_minitrough[which(Fst_peak_dxy_avg_pi_minitrough$pop == "V1"),]
SA_pop2 <- Fst_peak_dxy_avg_pi_minitrough[which(Fst_peak_dxy_avg_pi_minitrough$pop == "V2"),]

# Remove to clear space in the global environment
remove(Fst_peak_dxy_avg_pi, Fst_peak_dxy_avg_pi_minitrough)


# No situations follow the recurrent seleciton or geographic sweep before differentiation model (minimum z-score of pi was -2.12)

# Check to make sure that there are no overlapping windows between the Div_gen and SA objects
Div_gen_pop1_winds <- Div_gen_pop1[,1:4]
Div_gen_pop1_winds$comp <- "DG1"

Div_gen_pop2_winds <- Div_gen_pop2[,1:4]
Div_gen_pop2_winds$comp <- "DG2"

SA_pop1_winds <- SA_pop1[,1:4]
SA_pop1_winds$comp <- "SA1"

SA_pop2_winds <- SA_pop2[,1:4]
SA_pop2_winds$comp <- "SA2"

DG1_SA1 <- left_join(Div_gen_pop1_winds, SA_pop1_winds, by = join_by(chromosome, Plotting.name, window_pos_1, window_pos_2)) #0, good
DG2_SA1 <- left_join(Div_gen_pop2_winds, SA_pop1_winds, by = join_by(chromosome, Plotting.name, window_pos_1, window_pos_2)) #0, good

DG1_SA2 <- left_join(Div_gen_pop1_winds, SA_pop2_winds, by = join_by(chromosome, Plotting.name, window_pos_1, window_pos_2)) #0, good
DG2_SA2 <- left_join(Div_gen_pop2_winds, SA_pop2_winds, by = join_by(chromosome, Plotting.name, window_pos_1, window_pos_2)) #0, good

# Figure out which divergence with gene flow windows are shared
DG1_DG2 <- left_join(Div_gen_pop1_winds, Div_gen_pop2_winds, by = join_by(chromosome, Plotting.name, window_pos_1, window_pos_2))

Div_gen_shared <- DG1_DG2[which(!is.na(DG1_DG2$comp.y)),]

# Bind everything together so that we can isolate the neutral windows
Cands <- rbind(Div_gen_pop1_winds, Div_gen_pop2_winds, SA_pop1_winds, SA_pop2_winds)

Neutral <- anti_join(Fst_dxy_comb, Cands)

# write out results
write.csv(Div_gen_pop1_winds, "Modestum_divergencewgeneflow_pop1windows.csv", row.names = F)
write.csv(Div_gen_pop2_winds, "Modestum_divergencewgeneflow_pop2windows.csv", row.names = F)
write.csv(Div_gen_shared, "Modestum_divergencewgeneflow_sharedwindows.csv", row.names = F)

write.csv(SA_pop1_winds, "Modestum_selectioninallopatry_pop1windows.csv")
write.csv(SA_pop2_winds, "Modestum_selectioninallopatry_pop2windows.csv")

# write out results as bed files
write.table(Div_gen_pop1_winds[,c(1,3,4)], "Modestum_divergencewgeneflow_pop1windows.bed", row.names = F, col.names = F, sep = '\t', quote = F)
write.table(Div_gen_pop2_winds[,c(1,3,4)], "Modestum_divergencewgeneflow_pop2windows.bed", row.names = F, col.names = F, sep = '\t', quote = F)
write.table(Div_gen_shared[,c(1,3,4)], "Modestum_divergencewgeneflow_sharedwindows.bed", row.names = F, col.names = F, sep = '\t', quote = F)

write.table(SA_pop1_winds[,c(1,3,4)], "Modestum_selectioninallopatry_pop1windows.bed", row.names = F, col.names = F, sep = '\t', quote = F)
write.table(SA_pop2_winds[,c(1,3,4)], "Modestum_selectioninallopatry_pop2windows.bed", row.names = F, col.names = F, sep = '\t', quote = F)


# Make plots 

Fst_dat_toplot <- Fst_dat[,c(3,4,1,5,6,7)]
Dxy_dat_toplot <- Dxy_dat[,c(3,4,1,5,6,7)]

Pi_pop1$pop2 <- NA
colnames(Pi_pop1)[3] <- "pop1"
Pi_pop1_toplot <- Pi_pop1[,c(3,15,1,4,5,6)]
colnames(Pi_pop1_toplot)[6] <- "avg_pi_pop1"

Pi_pop2$pop2 <- NA
colnames(Pi_pop2)[3] <- "pop1"
Pi_pop2_toplot <- Pi_pop2[,c(3,15,1,4,5,6)]
colnames(Pi_pop2_toplot)[6] <- "avg_pi_pop2"

# Create a list to store everything 

Fst_dat_toplot <- Fst_dat_toplot %>% 
  gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
         key = "statistic", value = "value")

Dxy_dat_toplot <- Dxy_dat_toplot %>% 
  gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
         key = "statistic", value = "value")

Pi_pop1_toplot <- Pi_pop1_toplot %>% 
  gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
         key = "statistic", value = "value")

Pi_pop2_toplot <- Pi_pop2_toplot %>% 
  gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
         key = "statistic", value = "value")

Pixy_plot_df <- bind_rows(Fst_dat_toplot, Dxy_dat_toplot, Pi_pop1_toplot, Pi_pop2_toplot) %>% 
  arrange(pop1, pop2, chromosome, window_pos_1, statistic)

Pixy_plot_df_filt <- Pixy_plot_df[which(grepl("CM", Pixy_plot_df$chromosome) | Pixy_plot_df$chromosome == "JAIPUX010001232.1" | Pixy_plot_df$chromosome == "JAIPUX010001880.1"),]

# Reformat chromosomes to numeric 
Chrom_names <- unique(Pixy_plot_df$chromosome)

Chrom_names <- Chrom_names[c(1,3,17,18,6,8,10,2,4,5,7,9,11:16)]

Pixy_plot_df_filt$chromosome <- factor(Pixy_plot_df_filt$chromosome, levels = Chrom_names)

Pixy_plot_df_filt$Plottin.name <- Pixy_plot_df_filt$chromosome

Pixy_plot_df_filt$Plottin.color <- "no"

# Identify candidates 





# Only select statistics we want
Pixy_plot_df_filt <- Pixy_plot_df_filt[grepl("avg", Pixy_plot_df_filt$statistic),]

pixy_labeller <- as_labeller(c(avg_pi_pop1 = "pi[A]",
                               avg_pi_pop2 = "pi[B]",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)

Pixy_plot_df_filt %>%
  #mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
  #TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:18))) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = "#ebebeb"))+
  geom_point(size = 1.25, alpha = 0.25, stroke = 0)+
  facet_grid(statistic  ~ Plottin.name,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
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

# Get the candidate windows
setwd("/Users/kfarleigh/Desktop/Github/Modestum_WGS/Analyses/pixy")

pixy_cands <- list.files(pattern="*windows.csv")

cand_windows <- list()

for(i in 1:length(pixy_cands)){
  
  if(grepl("shared", pixy_cands[i])){
    tmp_dat <- read.csv(pixy_cands[i])
    tmp_dat$comp <- paste(tmp_dat$comp.x, tmp_dat$comp.y, sep = "_")
    tmp_dat <- tmp_dat[,c(1:4,7)]
    
    cand_windows[[i]] <- tmp_dat
    
    remove(tmp_dat)
  } else{
    
    cand_windows[[i]] <- read.csv(pixy_cands[i])
  }
}

names(cand_windows) <- gsub(".csv", "", pixy_cands)

# Bind into one data frame 

All_cand_windows <- do.call("rbind", cand_windows)


## Determine which windows are candidates 
Pixy_plot_df_filt2 <- Pixy_plot_df_filt

Pixy_plot_df_filt2 $index <- 1:nrow(Pixy_plot_df_filt2)

for(i in 1:nrow(All_cand_windows)){
  
  tmp_cands_dat <- All_cand_windows[i,]
  
  tmp_dat <- Pixy_plot_df_filt2 [which(Pixy_plot_df_filt2 $chromosome == tmp_cands_dat$chromosome),]
  
  
  if(length(which(between(tmp_dat$window_pos_1, tmp_cands_dat$window_pos_1, tmp_cands_dat$window_pos_2))) != 0){
    
    tmp_dat_filt <- tmp_dat[which(between(tmp_dat$window_pos_1, tmp_cands_dat$window_pos_1, tmp_cands_dat$window_pos_2)),]
    
    tmp_idx <- tmp_dat_filt$index
    
    Pixy_plot_df_filt2$Plottin.color[tmp_idx] <- tmp_cands_dat$comp
    
    remove(tmp_dat_filt, tmp_idx)
    
  }
  
  remove(tmp_cands_dat, tmp_dat)
  
  gc()
  
  print(paste("Finished with ", i, " of 257.", sep = ""))
  
}

table(Pixy_plot_df_filt2$Plottin.color)

Pixy_plot_df_filt2 <- Pixy_plot_df_filt2[grepl("avg", Pixy_plot_df_filt2$statistic),]


Pixy_plot_df_filt2 %>%
  mutate(chromosome = factor(chromosome, levels = c(1:18))) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = Plottin.color))+
  geom_point(data = Pixy_plot_df_filt2[Pixy_plot_df_filt2$Plottin.color == "no",] , size = 1.25, alpha = 0.25, stroke = 0) +
  geom_point(data = Pixy_plot_df_filt2[Pixy_plot_df_filt2$Plottin.color != "no",] , size = 1.25, alpha = 0.75, stroke = 0) +
  facet_grid(statistic  ~ Plottin.name,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
  scale_color_manual(values = c("#ef3e2c", "#ef3e2c", "#ef3e2c", "#686868", "#6bacd6", "#6bacd6")) + 
  xlab("Chromsome")+
  ylab("Statistic Value")+
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

# Try a line to reduce file size 

Pixy_plot_df_filt2 %>%
  mutate(chromosome = factor(chromosome, levels = c(1:18))) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = Plottin.color))+
  geom_line(alpha = 0.5, color = "#888888") +
  geom_point(data = Pixy_plot_df_filt2[which(Pixy_plot_df_filt2$Plottin.color == "DG1" | Pixy_plot_df_filt2$Plottin.color == "DG1_DG2" | 
                                               Pixy_plot_df_filt2$Plottin.color == "DG2" | Pixy_plot_df_filt2$Plottin.color =="SA1" | 
                                               Pixy_plot_df_filt2$Plottin.color == "SA2"),] , 
             size = 2.5, alpha = 0.75, stroke = 0) +
  facet_grid(factor(statistic, levels = c("avg_wc_fst", "avg_dxy", "avg_pi_pop2", "avg_pi_pop1"))  ~ Plottin.name,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
  scale_color_manual(values = c("#ef3e2c", "#ef3e2c", "#ef3e2c", "#6bacd6", "#6bacd6")) + 
  xlab("Chromsome")+
  ylab("Statistic Value")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none", 
        strip.text.x = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))


### Isolate end of chromosome 2
Pix_plot_df_chr2 <- Pixy_plot_df_filt2 %>% filter(chromosome == "CM034703.1")

Pix_plot_df_chr2 %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = Plottin.color)) +
  geom_line(alpha = 0.5, color = "#888888") +
  geom_point(data = Pix_plot_df_chr2[which(Pix_plot_df_chr2$Plottin.color == "DG1" | Pix_plot_df_chr2$Plottin.color == "DG1_DG2" | 
                                             Pix_plot_df_chr2$Plottin.color == "DG2" | Pix_plot_df_chr2$Plottin.color =="SA1" | 
                                               Pix_plot_df_chr2$Plottin.color == "SA2"),] , 
             size = 2.5, alpha = 0.75, stroke = 0) +
  facet_grid(factor(statistic, levels = c("avg_wc_fst", "avg_dxy", "avg_pi_pop2", "avg_pi_pop1"))  ~ Plottin.name,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
  scale_color_manual(values = c("#ef3e2c", "#ef3e2c", "#ef3e2c", "#6bacd6", "#6bacd6")) + 
  xlab("Chromsome")+
  ylab("Statistic Value")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none", 
        strip.text.x = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))

Pix_plot_df_chr2_fst <- Pix_plot_df_chr2[which(Pix_plot_df_chr2$statistic == "avg_wc_fst"),]

# Beginning of chromosome
Pix_plot_df_chr2_fst %>%
  filter(index < 93974) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value)) +
  geom_line(alpha = 0.5, color = "#888888", linewidth = 0.75) + 
  stat_smooth("method" = "lm", se = FALSE) +
  xlab("Chromsome 2")+
  ylab("Fst")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none", 
        strip.text.x = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))

# Is it statistically significant?
Pix_plot_df_chr2_fst %>%
       filter(index < 93974) %>% lm(value~window_pos_1, data = .) %>% summary()

# End of chromosome
Pix_plot_df_chr2_fst %>%
  filter(index > 144459) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value)) +
  geom_line(alpha = 0.5, color = "#888888", linewidth = 0.75) + 
  stat_smooth("method" = "lm", se = FALSE) + 
  xlab("Chromsome 2")+
  ylab("Fst")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none", 
        strip.text.x = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))

Pix_plot_df_chr2_fst %>%
  filter(index > 144459) %>% lm(value~window_pos_1, data = .) %>% summary()




######### Read in Tajima's D
setwd("C:/Users/keaka/Desktop/Github/Modestum_WGS/Analyses/tajimas_d")

North_TajD <- read.delim('north.TajD.10kb.txt')
South_TajD <- read.delim('south.TajD.10kb.txt')

colnames(North_TajD) <- c("chromosome", "window_pos_1", "n_snps", "Tajimas.D")
colnames(South_TajD) <- c("chromosome", "window_pos_1", "n_snps", "Tajimas.D")



# Make numeric and add 1 so that we can join them 
North_TajD$window_pos_1 <- as.numeric(North_TajD$window_pos_1) + 1
North_TajD$Tajimas.D <- as.numeric(North_TajD$Tajimas.D)

South_TajD$window_pos_1 <- as.numeric(South_TajD$window_pos_1) + 1
South_TajD$Tajimas.D <- as.numeric(South_TajD$Tajimas.D)

# Get the neutral windows
Neut_winds_north <- left_join(Neutral, North_TajD,  by = join_by(chromosome, window_pos_1))
Neut_winds_north$comp <- "Neutral"

Neut_winds_south <- left_join(Neutral, South_TajD, by = join_by(chromosome, window_pos_1))
Neut_winds_south$comp <- "Neutral"
# Join with the divergence and selection in allopatry windows, isolate windows that aren't (neutral windows)

North_TajD_DG1 <- left_join(Div_gen_pop1_winds, North_TajD , by = join_by(chromosome, window_pos_1))
North_TajD_SA1 <- left_join(SA_pop1_winds, North_TajD, by = join_by(chromosome, window_pos_1))


South_TajD_DG2 <- left_join(Div_gen_pop2_winds, South_TajD , by = join_by(chromosome, window_pos_1))
South_TajD_SA2 <- left_join(SA_pop2_winds,South_TajD, by = join_by(chromosome, window_pos_1))


### Compare them (after revisions)

North_DG1vsims <- rbind(North_TajD_DG1[,c(7,5)], td_north_dwgf_sims_df[,c(2,1)])

######################
#### Simulations #####
######################

setwd("~/Desktop/Github/Modestum_WGS/Analyses/tajimas_d")

North <- read.delim("north.TajD.10kb.txt")

colnames(North)[1:2] <- c("chromosome", "window_pos_1")

North$window_pos_1 <- as.numeric(North$window_pos_1)

North$TajimaD <- as.numeric(North$TajimaD)

North$window_pos_1 <- North$window_pos_1 + 1

# Get pi windows with z-scores in the divergence with gene flow range
Pi_pop1_sims <- Pi_pop1 %>% filter(z_score_bypop < -0.67 & z_score_bypop > -3)

North_Pi_pop1_sims <- left_join(North, Pi_pop1_sims)

North_Pi_pop1_sims <- North_Pi_pop1_sims[complete.cases(North_Pi_pop1_sims),]

# Set our seed
set.seed(010526)

# Randomly select the same number of the windows as we saw empirically (North = 34, South = 27)

sims <- 10000

td_north_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 34)
  
  td_north_dwgf_sims <- c(td_north_dwgf_sims, mean(sim_dat[,4]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

td_north_dwgf_sims_df <- data.frame(Tajimas.D = td_north_dwgf_sims, comp = "North_dwgf_sim")

### SIA simulations
# The z-score threshold is the same so we can reuse the one we created before

td_north_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- North_Pi_pop1_sims %>% slice_sample(n = 85)
  
  td_north_sia_sims <- c(td_north_sia_sims, mean(sim_dat[,4]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

td_north_sia_sims_df <- data.frame(Tajimas.D = td_north_sia_sims, comp = "North_sia_sim")


### Do it for the South 

setwd("~/Desktop/Github/Modestum_WGS/Analyses/tajimas_d")

South <- read.delim("south.TajD.10kb.txt")

colnames(South)[1:2] <- c("chromosome", "window_pos_1")

South$window_pos_1 <- as.numeric(South$window_pos_1)

South$TajimaD <- as.numeric(South$TajimaD)

South$window_pos_1 <- South$window_pos_1 + 1

# Get pi windows with z-scores in the divergence with gene flow range
Pi_pop2_sims <- Pi_pop2 %>% filter(z_score_bypop < -0.67 & z_score_bypop > -3)

South_Pi_pop2_sims <- left_join(South, Pi_pop2_sims)

South_Pi_pop2_sims <- South_Pi_pop2_sims[complete.cases(South_Pi_pop2_sims),]

# Set our seed
set.seed(010526)

# Randomly select the same number of the windows as we saw empirically (North = 34, South = 27)

sims <- 10000

td_south_dwgf_sims <- c()

for(j in 1:sims){
  
  sim_dat <- South_Pi_pop2_sims %>% slice_sample(n = 27)
  
  td_south_dwgf_sims <- c(td_south_dwgf_sims, mean(sim_dat[,4]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

td_south_dwgf_sims_df <- data.frame(Tajimas.D = td_south_dwgf_sims, comp = "South_dwgf_sim")

#SIA simulations

### Do it for the South
sims <- 10000

td_south_sia_sims <- c()

for(j in 1:sims){
  
  sim_dat <- South_Pi_pop2_sims %>% slice_sample(n = 90)
  
  td_south_sia_sims <- c(td_south_sia_sims, mean(sim_dat[,4]))
  
  remove(sim_dat)
  
  print(paste("Finished simulation ", j, " of ", sims, ".", sep = ""))
  
}

td_south_sia_sims_df <- data.frame(Tajimas.D = td_south_sia_sims, comp = "South_sia_sim")

# Kruskal-Wallis chi-squared = 40.031, df = 1, p-value < 2.5e-10
kruskal.test(North_DG1vsims$Tajimas.D, North_DG1vsims$comp)

# -1.343486
mean(North_TajD_DG1$Tajimas.D)

# -1.154151
mean(td_north_dwgf_sims_df$Tajimas.D)

North_SIA1vsims <- rbind(North_TajD_SA1[,c(7,5)], td_north_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 87.407, df = 1, p-value < 2.2e-10
kruskal.test(North_SIA1vsims$Tajimas.D, North_SIA1vsims$comp)

# -1.335911
mean(North_TajD_SA1$Tajimas.D)

# -1.153551
mean(td_north_sia_sims_df$Tajimas.D)

### South

South_DG2vsims <- rbind(South_TajD_DG2[,c(7,5)], td_south_dwgf_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 0.9124, df = 1, p-value = 0.3395
kruskal.test(South_DG2vsims$Tajimas.D, South_DG2vsims$comp)

# -0.1163855
mean(South_TajD_DG2$Tajimas.D)

# -0.1643019
mean(td_south_dwgf_sims_df$Tajimas.D)

South_SA2vsims <- rbind(South_TajD_SA2[,c(7,5)], td_south_sia_sims_df[,c(2,1)])

# Kruskal-Wallis chi-squared = 20.624, df = 1, p-value = 5.589e-06
kruskal.test(South_SA2vsims$Tajimas.D, South_SA2vsims$comp)

# -0.2124127
mean(South_TajD_SA2$Tajimas.D)

# -0.1652915
mean(td_south_sia_sims_df$Tajimas.D)



### Make plots

N_dwgf_boxplot <- ggplot(data = North_DG1vsims) + geom_boxplot(aes(x = comp, y = Tajimas.D, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c','#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c','#ababab')) +
  theme_classic() + ylab("Tajima's D") + xlab(NULL) + theme(legend.position = "none") + ylim(-2,-0.5)

N_sia_boxplot <- ggplot(data = North_SIA1vsims) + geom_boxplot(aes(x = comp, y = Tajimas.D, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ababab', '#6bacd6'), 0.5)) + scale_color_manual(values = c('#ababab', '#6bacd6')) +
  theme_classic() + ylab("Tajima's D") + xlab(NULL) + theme(legend.position = "none") + ylim(-2,-0.5)

S_dwgf_boxplot <- ggplot(data = South_DG2vsims) + geom_boxplot(aes(x = comp, y = Tajimas.D, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#ef3e2c','#ababab'), 0.5)) + scale_color_manual(values = c('#ef3e2c','#ababab')) +
  theme_classic() + ylab("Tajima's D") + xlab(NULL) + theme(legend.position = "none") + ylim(-1.5,1.5)

S_sia_boxplot <- ggplot(data = South_SA2vsims) + geom_boxplot(aes(x = comp, y = Tajimas.D, fill = comp, color = comp), width = 0.2, outliers = FALSE) + 
  scale_fill_manual(values = alpha(c('#6bacd6', '#ababab'), 0.5)) + scale_color_manual(values = c('#6bacd6', '#ababab')) +
  theme_classic() + ylab("Tajima's D") + xlab(NULL) + theme(legend.position = "none") + ylim(-1.5,1.5)

pageCreate(width = 8.5, height = 11)

plotGG(N_dwgf_boxplot, x = 0.1, y = 0.1, width = 3, height = 3)

plotGG(N_sia_boxplot, x = 3.2, y = 0.1, width = 3, height = 3)

plotGG(S_dwgf_boxplot, x = 0.1, y = 3.2, width = 3, height = 3)

plotGG(S_sia_boxplot, x = 3.2, y = 3.2, width = 3, height = 3)

pageGuideHide()
