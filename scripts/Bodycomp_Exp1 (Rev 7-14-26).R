# GLP-1 experiment 1
#EchoMRI data after BW loss for NZO mice
#Exp 1: compare 2 GLP-1 agonists. Started with 12 wk old female NZO mice
#BW loss occurred over 35 days

#Lean mass data requested in 7-1-26 lab meeting by Laurie and Cathy

#Libraries ####
library(mmand)
library(pacman)
library(this.path)
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels
library(slider)
library(lubridate)
library(lme4)

#Format plot (LM version) ####
format.plot_LM <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"), 
  panel.border = element_blank(),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black"),
  plot.title = element_text(size=17, hjust = 0.5, face="bold", vjust=2),
  legend.title=element_text(size=15, face="bold"),
  legend.text=element_text(size=13),
  axis.title.x = element_blank(),
  axis.text.x = element_text(size= 13, angle=25, vjust=0.5, hjust=0.7),
  axis.title.y = element_text(face="bold", size= 15),
  axis.text.y = element_text(size = 13))
# Define custom colors
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb","Survodutide" = "#c93618", "Vehicle" = "#403d3c", "Naive" ="darkgray")
custom_colors_GLPB <- c("Tirzepatide" = "lightblue","Survodutide" = "lightgreen", "Vehicle" = "darkgray")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#read in metadata
metadata <- read_csv("../data/META.csv") %>% 
select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# echo MRI
echomri_csv_files <- tibble(
  filepath = list.files("../data/echoMRI", full.names = TRUE)) %>% 
  filter(grepl("*.xlsx", filepath)) 
echomri_csv_files

echomri_open_files <- echomri_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map(., function(X){
    readxl::read_xlsx(X$filepath) %>% 
      select(Label, Fat, Lean, Weight, TimeDateDura) %>% 
      rename(ID = Label) %>% 
      separate_wider_delim(TimeDateDura, delim = ";", names = c("Date", "A", "B")) %>% 
      select(-A, -B) %>% 
      separate_wider_delim(Date, delim = " ", names = c("hms", "month", "day", "year")) %>% 
      mutate(day = gsub(",", "", day),
             Date = paste(year, month, day, sep = "-"),
             Date = lubridate::ymd(Date),
             ID =  as.factor(ID)) %>% 
      select(-hms, -month, -day, -year)}) %>% 
  bind_rows() %>% 
  left_join(., metadata, by = "ID")
echomri_open_files

# compare adiposity index = fat / lean
echomri_data <- echomri_open_files %>% 
  mutate(adiposity_index = Fat / Lean) %>% 
  group_by(ID) %>% 
  mutate(n_measurement = as.numeric(as.factor(Date)))
echomri_data

write_csv(x = echomri_data, "../data/echomri.csv")

##The code below will pull data from echoMRI.csv
#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_data_19 <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) 

echoMRI_data_BWloss <- echoMRI_data_19 %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                          ID %in% c(3730, 3731, 3738, 3741) ~ "Tirzepatide",
                          ID %in% c(3732, 3733, 3737, 3739) ~ "Survodutide"),
         DRUG = as.factor(DRUG),
         n_measurement = case_when(n_measurement == "1" ~ "Pre BW loss", # Day 0
                                   n_measurement == "2" ~ "During BW loss", # Day 21
                                   n_measurement == "3" ~ "End BW loss", # Day 34
                                   n_measurement == "4" ~ "Post BW regain")) %>%
  mutate(
    n_measurement = factor(n_measurement, 
                           levels = c("Pre BW loss", #pre-treatment
                                      "During BW loss", #During BW loss
                                      "End BW loss", #End BW loss
                                      "Post BW regain"))) #BW regain

# Calculations  ####
#change in lean, fat, total body mass, adiposity index
echmoMRI_Exp1_delta <- echoMRI_data_BWloss %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(FatplusLean = Fat + Lean) %>%
  mutate(delta_lean = Lean - first (Lean), # negative number means lean mass was lost
         delta_fat = Fat - first(Fat),
         delta_FatplusLean = FatplusLean - first(FatplusLean),
         delta_AI= adiposity_index - first(adiposity_index),
         Lean_FatplusLean = delta_lean/delta_FatplusLean, #I think this isn't useful
         delta_lean_vs_delta_BW = 100*(delta_lean/delta_FatplusLean), # I think this is the percent of BW loss which is lean mass
         BW_pct_change = 100*((FatplusLean - first(FatplusLean)) / first(FatplusLean)),
         Lean_pct_change = 100*((Lean - first(Lean)) / first(Lean)),
         Fat_pct_change = 100*((Fat - first(Fat)) / first(Fat)),
         Lean_pct_BW = 100*(Lean/FatplusLean),
         Fat_pct_BW = 100*(Fat/FatplusLean))

plot_echo_delta <-echmoMRI_Exp1_delta %>%
  filter_out(n_measurement=="Day 0") %>%
  ungroup() %>%
  group_by(DRUG, n_measurement) %>%
  summarise(Avg_Fat_G =mean(Fat),
            SD_Fat_G =sd(Fat),
            Avg_Lean_G=mean(Lean),
            SD_Lean_G=sd(Lean),
            Avg_BW_G=mean(FatplusLean),
            SD_BW_G=sd(FatplusLean),
    Avg_delta_Lean_G = mean(delta_lean),
            SD_delta_Lean_G = sd(delta_lean),
            Avg_delta_Fat_G = mean(delta_fat),
            SD_delta_Fat_G = sd(delta_fat),
            Avg_delta_FatplusLean = mean(delta_FatplusLean),
            SD_delta_FatplusLean = sd(delta_FatplusLean),
            Avg_Lean_FatplusLean = mean(Lean_FatplusLean),
            SD_Lean_FatplusLean = sd(Lean_FatplusLean),
            Avg_Lean_pct_change = mean(Lean_pct_change),
            SD_Lean_pct_change = sd(Lean_pct_change),
            Avg_Lean_pct_BW = mean(Lean_pct_BW),
            SD_Lean_pct_BW = sd(Lean_pct_BW)) %>%
  select(n_measurement, Avg_delta_Lean_G, SD_delta_Lean_G, 
         Avg_delta_FatplusLean, SD_delta_FatplusLean, 
         Avg_Lean_FatplusLean, SD_Lean_FatplusLean, Avg_Lean_pct_change, SD_Lean_pct_change, 
         Avg_Lean_pct_BW, SD_Lean_pct_BW)

#Summary of mass, lean, and fat in grams
summary_grams <-echmoMRI_Exp1_delta %>%
  filter_out(n_measurement=="Day 0") %>%
  ungroup() %>%
  group_by(DRUG, n_measurement) %>%
  summarise(Avg_Fat_G =mean(Fat),
            SD_Fat_G =sd(Fat),
            Avg_Lean_G=mean(Lean),
            SD_Lean_G=sd(Lean),
            #SE_Lean_G=(SD_Lean_G)/n(ID),
            Avg_BW_G=mean(FatplusLean),
            SD_BW_G=sd(FatplusLean))

#---
#---
#Body weight (g) ####
##Graph BW (g) ####
#x axis as n_measurement and y axis as Fat+Lean
BW_plot_G <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=FatplusLean, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Body weight (g)",
       title= "Body weight (g)",
       color="Treatment", fill="Treatment")
BW_plot_G

#Export plot to a folder on laptop called "figures" 
ggsave(BW_plot_G,
       filename="Exp1_BW_plot_G.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: BW (g) ####
#Build multiple linear regression model for BW_pct_change #
model_BW_G <- lmer(FatplusLean ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_BW_G)

#Calculate estimated marginal means #
emm_BW_G <- emmeans(model_BW_G, ~ n_measurement*DRUG, cov.reduce = mean)
emm_BW_G_df <- as.data.frame(emm_BW_G)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_BW_G <- contrast(emm_BW_G, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_BW_G_df <- as.data.frame(contrasts_by_DRUG_BW_G)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_BW_G <- contrast(emm_BW_G, method = "pairwise", by = "n_measurement")
contrasts_by_n_measurement_BW_G_df <- as.data.frame(contrasts_by_n_measurement_BW_G)


#Lean mass (g) ####
##Graph lean (g) ####
#Graph with x axis as n_measurement and y axis as Lean
Lean_plot_G <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Lean, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Lean mass (g)",
       title= "Lean mass (g)",
       color="Treatment", fill="Treatment")
Lean_plot_G

#Export plot to a folder on laptop called "figures" 
ggsave(Lean_plot_G,
       filename="Exp1_Lean_plot_G.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: Lean mass (g) ####
#Build multiple linear regression model for Lean #
model_Lean_G <- lmer(Lean ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Lean_G)

#Calculate estimated marginal means #
emm_Lean_G <- emmeans(model_Lean_G, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Lean_G_df <- as.data.frame(emm_Lean_G)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Lean_G <- contrast(emm_Lean_G, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Lean_G_df <- as.data.frame(contrasts_by_DRUG_Lean_G)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Lean_G <- contrast(emm_Lean_G, method = "pairwise", by = "n_measurement")
contrasts_by_n_measurement_Lean_G_df <- as.data.frame(contrasts_by_n_measurement_Lean_G)

#Fat mass (g) ####
##Graph fat mass (g) ####
#Graph with x axis as n_measurement and y axis as Fat
Fat_plot_G <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Fat, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Fat mass (g)",
       title= "Fat mass (g)",
       color="Treatment", fill="Treatment")
Fat_plot_G

#Export plot to a folder on laptop called "figures" 
ggsave(Fat_plot_G,
       filename="Exp1_Fat_plot_G.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: Fat mass(g) ####
#Build multiple linear regression model for Fat_pct_change #
model_Fat_G <- lmer(Fat ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Fat_G)

#Calculate estimated marginal means #
emm_Fat_G <- emmeans(model_Fat_G, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Fat_G_df <- as.data.frame(emm_Fat_G)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Fat_G <- contrast(emm_Fat_G, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Fat_G_df <- as.data.frame(contrasts_by_DRUG_Fat_G)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Fat_G <- contrast(emm_Fat_G, method = "pairwise", by = "n_measurement")
contrasts_by_n_measurement_Fat_G_df <- as.data.frame(contrasts_by_n_measurement_Fat_G)

#---
#---
# Adiposity index ####
##Graph AI ####
AI_plot_Exp1 <-ggplot(echoMRI_data_BWloss, aes(x=n_measurement, y=adiposity_index, 
                                          group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Adiposity index (Fat/Lean)",
       title= "Adiposity index",
       color="Treatment", fill="Treatment")
AI_plot_Exp1

#Export plot to a folder on laptop called "figures" 
ggsave(AI_plot_Exp1,
       filename="Exp1_AI_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: Adiposity index (from initial fat mass)fat/lean) ####
#Build multiple linear regression model for Fat_pct_change #
model_AI <- lmer(adiposity_index ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_AI)

#Calculate estimated marginal means #
emm_AI <- emmeans(model_AI, ~ n_measurement*DRUG, cov.reduce = mean)
emm_AI_df <- as.data.frame(emm_AI)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_AI <- contrast(emm_AI, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_AI_df <- as.data.frame(contrasts_by_DRUG_AI)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_AI <- contrast(emm_AI, method = "pairwise", by = "n_measurement")
contrasts_by_n_measurement_AI_df <- as.data.frame(contrasts_by_n_measurement_AI)

#---
#---
# Δ BW (%) Percent change in BW from initial BW ####
## Graph ΔBW (%) from initial BW ####
BW_pct_change_plot <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=BW_pct_change, 
                                                     group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "Δ Body weight (%)",
       title= "Δ Body weight (%)",
       color="Treatment", fill="Treatment")
BW_pct_change_plot

#Export plot to a folder on laptop called "figures" 
ggsave(BW_pct_change_plot,
       filename="Exp1_BW_plot1.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔBW (%) from initial BW ####
#Build multiple linear regression model for BW_pct_change #
model_BW_pct_change <- lmer(BW_pct_change ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_BW_pct_change)

#Calculate estimated marginal means #
emm_BW_pct_change <- emmeans(model_BW_pct_change, ~ n_measurement*DRUG, cov.reduce = mean)
emm_BW_pct_change_df <- as.data.frame(emm_BW_pct_change)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_BW_pct_change <- contrast(emm_BW_pct_change, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_BW_pct_change_df <- as.data.frame(contrasts_by_DRUG_BW_pct_change)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_BW_pct_change <- contrast(emm_BW_pct_change, method = "pairwise", by = "n_measurement")
contrasts_n_measurement_BW_pct_change_df <- as.data.frame(contrasts_by_n_measurement_BW_pct_change)


#---
#---
# ΔLean mass (%) Percent change in lean mass from initial lean mass ####
## Graph ΔLean mass (%) from initial lean mass ####
Lean_pct_change_plot <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Lean_pct_change, 
                                                       group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "Δ Lean mass (%)",
       title= "Δ Lean mass (%)",
       color="Treatment", fill="Treatment")
Lean_pct_change_plot

#Export plot to a folder on laptop called "figures" 
ggsave(Lean_pct_change_plot,
       filename="Exp1_lean_plot1.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔLean mass (%) from initial lean mass ####
#Build multiple linear regression model for Lean_pct_change #
model_Lean_pct_change <- lmer(Lean_pct_change ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Lean_pct_change)

#Calculate estimated marginal means #
emm_Lean_pct_change <- emmeans(model_Lean_pct_change, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Lean_pct_change_df <- as.data.frame(emm_Lean_pct_change)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Lean_pct_change <- contrast(emm_Lean_pct_change, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Lean_pct_change_df <- as.data.frame(contrasts_by_DRUG_Lean_pct_change)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Lean_pct_change <- contrast(emm_Lean_pct_change, method = "pairwise", by = "n_measurement")
contrasts_n_measurement_Lean_pct_change_df <- as.data.frame(contrasts_by_n_measurement_Lean_pct_change)

#---
#---
# ΔFat mass (%) Percent change in fat mass from initial fat mass ####
## Graph ΔFat mass (%) from initial fat mass ####
Fat_pct_change_plot <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Fat_pct_change, 
                                                      group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "ΔFat (%)",
       title= "Percent change in fat mass",
       color="Treatment", fill="Treatment")
Fat_pct_change_plot

#Export plot to a folder on laptop called "figures" 
ggsave(Fat_pct_change_plot,
       filename="Exp1_fat_plot1.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔFat mass (%) from initial fat mass ####
#Build multiple linear regression model for Fat_pct_change #
model_Fat_pct_change <- lmer(Fat_pct_change ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Fat_pct_change)

#Calculate estimated marginal means #
emm_Fat_pct_change <- emmeans(model_Fat_pct_change, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Fat_pct_change_df <- as.data.frame(emm_Fat_pct_change)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Fat_pct_change <- contrast(emm_Fat_pct_change, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Fat_pct_change_df <- as.data.frame(contrasts_by_DRUG_Fat_pct_change)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Fat_pct_change <- contrast(emm_Fat_pct_change, method = "pairwise", by = "n_measurement")
contrasts_n_measurement_Fat_pct_change_df <- as.data.frame(contrasts_by_n_measurement_Fat_pct_change)


#---
#---
# % of BW comprised of lean mass ####
## Graph % of BW comprised of lean mass ####
Lean_pct_of_BW_plot <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Lean_pct_BW, 
                                                      group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "% Lean mass (Lean/BW)",
       title= "% Lean mass",
       color="Treatment", fill="Treatment")
Lean_pct_of_BW_plot

#Export plot to a folder on laptop called "figures" 
ggsave(Lean_pct_of_BW_plot,
       filename="Exp1_lean_plot2.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

##LMM: % of BW comprised of lean mass ####
#Build multiple linear regression model for Lean_pct_of_BW #
model_Lean_pct_of_BW <- lmer(Lean_pct_BW ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Lean_pct_of_BW)

#Calculate estimated marginal means #
emm_Lean_pct_of_BW <- emmeans(model_Lean_pct_of_BW, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Lean_pct_of_BW_df <- as.data.frame(emm_Lean_pct_of_BW)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Lean_pct_of_BW <- contrast(emm_Lean_pct_of_BW, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Lean_pct_of_BW_df <- as.data.frame(contrasts_by_DRUG_Lean_pct_of_BW)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Lean_pct_of_BW <- contrast(emm_Lean_pct_of_BW, method = "pairwise", by = "n_measurement")
contrasts_n_measurement_Lean_pct_of_BW_df <- as.data.frame(contrasts_by_n_measurement_Lean_pct_of_BW)

#---

#% of BW comprised of fat mass ####
##Graph % of BW comprised of fat mass ####
Fat_pct_of_BW_plot <-ggplot(echmoMRI_Exp1_delta, aes(x=n_measurement, y=Fat_pct_BW, 
                                group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "% Fat mass (Fat/BW)",
       title= "% Fat mass",
       color="Treatment", fill="Treatment")
Fat_pct_of_BW_plot

#Export plot to a folder on laptop called "figures" 
ggsave(Fat_pct_of_BW_plot,
       filename="Exp1_fat_plot2.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

##LMM: % of BW comprised of fat mass ####

#Build multiple linear regression model for Fat_pct_of_BW #
model_Fat_pct_of_BW <- lmer(Fat_pct_BW ~ n_measurement*DRUG + (1 | ID), data=echmoMRI_Exp1_delta)
summary(model_Fat_pct_of_BW)

#Calculate estimated marginal means #
emm_Fat_pct_of_BW <- emmeans(model_Fat_pct_of_BW, ~ n_measurement*DRUG, cov.reduce = mean)
emm_Fat_pct_of_BW_df <- as.data.frame(emm_Fat_pct_of_BW)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Fat_pct_of_BW <- contrast(emm_Fat_pct_of_BW, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Fat_pct_of_BW_df <- as.data.frame(contrasts_by_DRUG_Fat_pct_of_BW)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_n_measurement_Fat_pct_of_BW <- contrast(emm_Fat_pct_of_BW, method = "pairwise", by = "n_measurement")
contrasts_n_measurement_Fat_pct_of_BW_df <- as.data.frame(contrasts_by_n_measurement_Fat_pct_of_BW)

