# GLP-1 experiment 2A and 2B

#EchoMRI data after BW loss for NZO mice
#Exp 2A: dose escalation 10nmol TZP -> 20nmol TZP, each for one week in NZO mice
#Exp 2B: dose escalation 10nmol TZP -> 15nmol TZP, each for one week in NZO mice
  #Exp 2B was an exploration and did not have a control group

#Started with obese female NZO mice
#BW loss occurred over 14 days

#Lean mass data requested in 7-1-26 lab meeting by Laurie and Cathy

#Exp 2A: two groups -> TZP and vehicle
  #1st echo was when they reached adulthood
  #2nd echo was pre sable acclimation
  #3rd echo was immediately before start of BW loss (Peak obesity)
  #4th echo was day before BW loss sable recording started for TZP mice. Midway through 
      #saline injections for vehicle mice ocurred on the same date as 4th echo for TZP mice though.
  #5th echo was BW regain for TZP mice. Vehicle mice had been moved on to exp 2B
  #so I don't have a BW regain for vehicle mice in exp 2A

#Exp 2B: one group -> TZP mice. Use mice that had vehicle during Exp 2A
    # Baseline could be considered the 4th echo
    #End BW loss could be considered the 5th echo (July 6th or 8th)
    #End regain in last week of July --> This will complete data set for exp 2b

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

#Format plot (LM version)
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
custom_colors_GLP_Exp2 <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

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


# Exp. 2A ####
echoMRI_exp2a <- echoMRI_data_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3744,3745,3746,3748,3749,3752) ~ "Vehicle",
                          ID %in% c(3742,3743,3747,3750,3751,3753) ~ "Tirzepatide"),
         DRUG = as.factor(DRUG),
         Measurement = case_when( n_measurement == "1" ~ "Baseline I", 
                                  n_measurement == "2" ~ "Baseline II", 
                                   n_measurement == "3" ~ "Peak obesity", 
                                   n_measurement == "4" ~ "BW loss",
                                   n_measurement == "5" & 
                                   ID %in% c(3742,3743,3747,3750,3751,3753)~ "BW regain")) %>%
  drop_na(Measurement) %>%
  mutate(Measurement = factor(Measurement, 
                           levels = c("Baseline I", 
                                      "Baseline II", 
                                      "Peak obesity", 
                                      "BW loss",      
                                      "BW regain"))) 

# Calculate change in BW, lean, fat, & adiposity index ####
echmoMRI_exp2a_delta <- echoMRI_exp2a %>%
  filter(Measurement %in% c("Peak obesity", "BW loss", "BW regain"))%>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(BW = Fat + Lean) %>%
  mutate(delta_lean = Lean - first (Lean), # negative number means lean mass was lost
         delta_fat = Fat - first(Fat),
         delta_BW = BW - first(BW),
         
         Lean_BW = delta_lean/delta_BW, 
         delta_lean_vs_delta_BW = 100*(Lean_BW), # % of BW loss which is lean mass
         
         BW_pct_change = 100*((BW - first(BW)) / first(BW)),         # % change in BW
         Lean_pct_change = 100*((Lean - first(Lean)) / first(Lean)), # % change in lean mass
         Fat_pct_change = 100*((Fat - first(Fat)) / first(Fat)),     # % change in fat mass

         Lean_pct_BW = 100*(Lean/BW), # % of BW comprised of lean mass
         Fat_pct_BW = 100*(Fat/BW))   # % of BW comprised of fat mass

#### BW (g) ####
#Graph with x axis as Measurement and y axis as BW (g)...BW=Fat+Lean
BW_plot_exp2a_data <- echoMRI_exp2a %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(BW = Fat + Lean)
  
BW_plot_exp2a <-ggplot(BW_plot_exp2a_data, aes(x=Measurement, y=BW, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary",fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary",fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary",fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP_Exp2) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Body weight (g)", title= "Body Weight (g)", color="Treatment", fill="Treatment")
BW_plot_exp2a
  

####Lean mass (g) ####
#Graph with x axis as Measurement and y axis as Lean
Lean_plot_exp2a <-ggplot(echoMRI_exp2a, aes(x=Measurement, y=Lean, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP_Exp2) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Lean mass (g)", title= "Lean mass (g)", color="Treatment", fill="Treatment")
Lean_plot_exp2a


####Fat mass (g) ####
#Graph with x axis as n_measurement and y axis as Fat
Fat_plot_exp2a <-ggplot(echoMRI_exp2a, aes(x=Measurement, y=Fat, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary",fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP_Exp2) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "Fat mass (g)", title= "Fat mass (g)", color="Treatment", fill="Treatment")
Fat_plot_exp2a

#---
#---
# Adiposity index ####
##Graph AI ####
AI_plot_exp2a <-ggplot(echoMRI_exp2a, aes(x=Measurement, y=adiposity_index, 
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
  labs(y= "Adiposity index (fat/lean)",
       title= "Exp 2a: Adiposity index",
       color="Treatment", fill="Treatment")
AI_plot_exp2a

#---
#---
# ΔBW (%) Percent change in BW from initial BW ####
## Graph ΔBW (%) ####
BW_pct_change_exp2a <-ggplot(echmoMRI_exp2a_delta, aes(x=Measurement, y=BW_pct_change, 
                                                  group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "Δ Body weight (%)", title= "Change in body weight (%)", color="Treatment", fill="Treatment")
BW_pct_change_exp2a

#Export plot to a folder on laptop called "figures" 
ggsave(BW_pct_change_exp2a,
       filename="Exp2a_BW_plot2.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔBW (%) ####
#Build multiple linear regression model for BW_pct_change #
model_BW_pct_change_2a <- lmer(BW_pct_change ~ Measurement*DRUG + (1 | ID), data=echmoMRI_exp2a_delta)
summary(model_BW_pct_change_2a)

#Calculate estimated marginal means #
emm_BW_pct_change_2a <- emmeans(model_BW_pct_change_2a, ~ Measurement*DRUG, cov.reduce = mean)
emm_BW_pct_change_2a_df <- as.data.frame(emm_BW_pct_change_2a)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_BW_pct_change_2a <- contrast(emm_BW_pct_change_2a, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_BW_pct_change_2a_df <- as.data.frame(contrasts_by_DRUG_BW_pct_change_2a)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_Measurement_BW_pct_change_2a <- contrast(emm_BW_pct_change_2a, method = "pairwise", by = "Measurement")
contrasts_by_Measurement_BW_pct_change_2a_df <- as.data.frame(contrasts_by_Measurement_BW_pct_change_2a)

#---
#---
# ΔLean mass (%) Percent change in lean mass from initial lean mass ####
## Graph ΔLean mass (%)  ####
Lean_pct_change_exp2a <-ggplot(echmoMRI_exp2a_delta, aes(x=Measurement, y=Lean_pct_change, 
                                                  group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "ΔLean mass (%)", title= "Percent change in lean mass", color="Treatment", fill="Treatment")
Lean_pct_change_exp2a

#Export plot to a folder on laptop called "figures" 
ggsave(Lean_pct_change_exp2a,
       filename="Exp2a_lean_plot2.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔLean mass (%) ####
#Build multiple linear regression model for Lean_pct_change #
model_Lean_pct_change_2a <- lmer(Lean_pct_change ~ Measurement*DRUG + (1 | ID), data=echmoMRI_exp2a_delta)
summary(model_Lean_pct_change_2a)

#Calculate estimated marginal means #
emm_Lean_pct_change_2a <- emmeans(model_Lean_pct_change_2a, ~ Measurement*DRUG, cov.reduce = mean)
emm_Lean_pct_change_2a_df <- as.data.frame(emm_Lean_pct_change_2a)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Lean_pct_change_2a <- contrast(emm_Lean_pct_change_2a, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Lean_pct_change_2a_df <- as.data.frame(contrasts_by_DRUG_Lean_pct_change_2a)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_Measurement_Lean_pct_change_2a <- contrast(emm_Lean_pct_change_2a, method = "pairwise", by = "Measurement")
contrasts_by_Measurement_Lean_pct_change_2a_df <- as.data.frame(contrasts_by_Measurement_Lean_pct_change_2a)

#---
#---
# ΔFat mass (%) Percent change in fat mass from initial fat mass ####
## Graph ΔFat mass (%) ####
Fat_pct_change_exp2a <-ggplot(echmoMRI_exp2a_delta, aes(x=Measurement, y=Fat_pct_change, 
                                                  group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  labs(y= "Δ Fat mass (%)", title= "Percent change in fat mass", color="Treatment", fill="Treatment")
Fat_pct_change_exp2a

#Export plot to a folder on laptop called "figures" 
ggsave(Fat_pct_change_exp2a,
       filename="Exp2a_fat_plot2.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

## LMM: ΔFat mass (%)  ####
#Build multiple linear regression model for Fat_pct_change #
model_Fat_pct_change_2a <- lmer(Fat_pct_change ~ Measurement*DRUG + (1 | ID), data=echmoMRI_exp2a_delta)
summary(model_Fat_pct_change_2a)

#Calculate estimated marginal means #
emm_Fat_pct_change_2a <- emmeans(model_Fat_pct_change_2a, ~ Measurement*DRUG, cov.reduce = mean)
emm_Fat_pct_change_2a_df <- as.data.frame(emm_Fat_pct_change_2a)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Fat_pct_change_2a <- contrast(emm_Fat_pct_change_2a, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Fat_pct_change_2a_df <- as.data.frame(contrasts_by_DRUG_Fat_pct_change_2a)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_Measurement_Fat_pct_change_2a <- contrast(emm_Fat_pct_change_2a, method = "pairwise", by = "Measurement")
contrasts_by_Measurement_Fat_pct_change_2a_df <- as.data.frame(contrasts_by_Measurement_Fat_pct_change_2a)

#---
#---
#Percent of BW comprised of lean mass ####
## Graph: Lean/BW (%) ####
Lean_pct_of_BW_exp2a <-ggplot(echmoMRI_exp2a_delta, aes(x=Measurement, y=Lean_pct_BW, 
                              group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4, 
             position = position_dodge(width = 0.15)) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1, 
            position = position_dodge(width = 0.15)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1,
                position = position_dodge(width = 0.15)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "% of total BW", title= "Percent of BW comprised of lean mass",
       color="Treatment", fill="Treatment")
Lean_pct_of_BW_exp2a

#Export plot to a folder on laptop called "figures" 
ggsave(Lean_pct_of_BW_exp2a,
       filename="Exp2a_lean_plot3.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

##LMM: Lean/BW (%) ####
#Build multiple linear regression model for Lean_pct_of_BW #
model_Lean_pct_of_BW_2a <- lmer(Lean_pct_BW ~ Measurement*DRUG + (1 | ID), data=echmoMRI_exp2a_delta)
summary(model_Lean_pct_of_BW_2a)

#Calculate estimated marginal means #
emm_Lean_pct_of_BW_2a <- emmeans(model_Lean_pct_of_BW_2a, ~ Measurement*DRUG, cov.reduce = mean)
emm_Lean_pct_of_BW_2a_df <- as.data.frame(emm_Lean_pct_of_BW_2a)

# Pairwise contrasts within each DRUG (treatment group)
contrasts_by_DRUG_Lean_pct_of_BW_2a <- contrast(emm_Lean_pct_of_BW_2a, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Lean_pct_of_BW_2a_df <- as.data.frame(contrasts_by_DRUG_Lean_pct_of_BW_2a)

# Pairwise contrasts within each n_measurement (time point)
contrasts_by_Measurement_Lean_pct_of_BW_2a <- contrast(emm_Lean_pct_of_BW_2a, method = "pairwise", by = "Measurement")
contrasts_by_Measurement_Lean_pct_of_BW_2a_df <- as.data.frame(contrasts_by_Measurement_Lean_pct_of_BW_2a)

#---
#---
#Percent of BW comprised of fat mass ####
## Graph: Fat/BW (%) ####
Fat_pct_of_BW_exp2a <-ggplot(echmoMRI_exp2a_delta, aes(x=Measurement, y=Fat_pct_BW, 
                              group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4, 
             position = position_dodge(width = 0.15)) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1, 
            position = position_dodge(width = 0.15)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1, 
            position = position_dodge(width = 0.15)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  labs(y= "% of total BW", title= "Percent of BW comprised of fat mass",
       color="Treatment", fill="Treatment")
Fat_pct_of_BW_exp2a

#Export plot to a folder on laptop called "figures" 
ggsave(Fat_pct_of_BW_exp2a,
       filename="Exp2a_fat_plot3.png", 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

##LMM: Fat/BW (%) ####
#Build multiple linear regression model for Fat_pct_of_BW #
model_Fat_pct_of_BW_2a <- lmer(Fat_pct_BW ~ Measurement*DRUG + (1 | ID), data=echmoMRI_exp2a_delta)
summary(model_Fat_pct_of_BW_2a)

#Calculate estimated marginal means #
emm_Fat_pct_of_BW_2a <- emmeans(model_Fat_pct_of_BW_2a, ~ Measurement*DRUG, cov.reduce = mean)
emm_Fat_pct_of_BW_2a_df <- as.data.frame(emm_Fat_pct_of_BW_2a)

# Pairwise contrasts within each DRUG (treatment group) #
contrasts_by_DRUG_Fat_pct_of_BW_2a <- contrast(emm_Fat_pct_of_BW_2a, method = "pairwise", by = "DRUG")
contrasts_by_DRUG_Fat_pct_of_BW_2a_df <- as.data.frame(contrasts_by_DRUG_Fat_pct_of_BW_2a)

# Pairwise contrasts within each n_measurement (time point) #
contrasts_by_Measurement_Fat_pct_of_BW_2a <- contrast(emm_Fat_pct_of_BW_2a, method = "pairwise", by = "Measurement")
contrasts_by_Measurement_Fat_pct_of_BW_2a_df <- as.data.frame(contrasts_by_Measurement_Fat_pct_of_BW_2a)


