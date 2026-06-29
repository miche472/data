#For APS poster

#Created:3-27-26
#Revised:3-27-26

#Figure 1: body weight and composition during weight cycling
  #Graphs of body composition and BW at each Sable time point
  #Statistical analysis using unpaired t test


#Libraries
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

#---------------------------------#
#Graphs of Mass (Fat, lean, and total) - one graph per time point
#unpaired t test comparing the two groups at each time point

#Prepare data
#Import echo data
#Process EchoMRI data 
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #tech issue with sable cages
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(day_rel = Date - first(Date),
         STATUS = case_when(
           n_measurement == 1 ~ "Baseline",
           Date %in% as.Date(c("2025-01-27", "2025-01-29", "2025-02-07")) ~ "Peak obesity",
           ID %in% c(3711, 3727) & Date == as.Date("2025-02-20") ~ "Peak obesity",
           Date == as.Date("2025-03-28") ~ "BW loss",
           Date == as.Date("2025-05-27") ~ "BW maintenance",
           Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                               "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
           TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3711 & Date == as.Date("2025-01-27"))) %>% #repeated
  filter(!(ID == 3727 & Date == as.Date("2025-02-07"))) %>% #repeated
  mutate(STATUS = factor(STATUS,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  rename(SABLE = STATUS) 


# Convert echoMRI_data to long format
echoMRI_data_long <- echoMRI_data %>%
  group_by(ID, SABLE, Date) %>%
  rename(mass_Fat= Fat,
         mass_Lean=Lean,
         mass_Total = Weight,
         mass_adiposity_index = adiposity_index) %>%
  select(ID, 
         GROUP, 
         SABLE,
         mass_Fat,
         mass_Lean, 
         mass_Total, 
         mass_adiposity_index) %>%
  pivot_longer(
    cols = starts_with("mass_"),
    names_to = "Component",
    values_to = "Mass_G") %>%
  mutate(
    Component = factor(
      Component,
      levels = c("mass_Fat",
                 "mass_Lean",
                 "mass_Total",
                 "mass_adiposity_index"),
      labels = c("Fat",         
                 "Lean",                 
                 "Total",          
                 "Adiposity index")))             

# Graphs ####
#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
custom_colors2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")

##Graph baseline ####
#Make dataframe with echoMRI taken at baseline --
echoMRI_Baseline <- echoMRI_data_long %>%
  filter(SABLE == "Baseline") %>%
  filter(!Component == "Adiposity index")  %>%
  group_by(GROUP, Component)

ggplot(echoMRI_Baseline, aes(x = Component, y = Mass_G, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        axis.title.x = element_blank()) +
  format.plot +
  ylim(0,60) +
  labs(
    title="Figure 1A: Body composition at Baseline",
    y = "Mass (g)",
    fill = "Treatment group",
    color = "Treatment group")

#---
##Graph Peak obesity ####
#Make dataframe with echoMRI taken at Peak obesity
echoMRI_obesity <- echoMRI_data_long %>%
  filter(SABLE == "Peak obesity") %>%
  filter(!Component == "Adiposity index")  %>%
  group_by(GROUP, Component)

ggplot(echoMRI_obesity, aes(x = Component, y = Mass_G, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        axis.title.x = element_blank()) +
  format.plot +
  ylim(0,60) +
  labs(
    title="Figure 1B: Body composition at Peak obesity",
    y = "Mass (g)",
    fill = "Treatment group",
    color = "Treatment group")

#---
##Graph BW loss ####
#Make dataframe with echoMRI taken at BW loss
echoMRI_BWloss <- echoMRI_data_long %>%
  filter(SABLE == "BW loss") %>%
  filter(!Component == "Adiposity index")  %>%
  group_by(GROUP, Component)

ggplot(echoMRI_BWloss, aes(x = Component, y = Mass_G, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        axis.title.x = element_blank()) +
  format.plot +
  ylim(0,60) +
  labs(
    title="Figure 1C: Body composition at BW loss",
    y = "Mass (g)",
    fill = "Treatment group",
    color = "Treatment group")

#---
##Graph BW maintenance ####
#Make dataframe with echoMRI taken at BW maintenance
echoMRI_BWmaintenance <- echoMRI_data_long %>%
  filter(SABLE == "BW maintenance") %>%
  filter(!Component == "Adiposity index")  %>%
  group_by(GROUP, Component)

ggplot(echoMRI_BWmaintenance, aes(x = Component, y = Mass_G, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        axis.title.x = element_blank()) +
  format.plot +
  ylim(0,60) +
  labs(
    title="Figure 1D: Body composition at BW maintenance",
    y = "Mass (g)",
    fill = "Treatment group",
    color = "Treatment group")

#---
##Graph BW regain ####
#Make dataframe with echoMRI taken at BW regain
echoMRI_BWregain <- echoMRI_data_long %>%
  filter(SABLE == "BW regain") %>%
  filter(!Component == "Adiposity index")  %>%
  group_by(GROUP, Component)

ggplot(echoMRI_BWregain, aes(x = Component, y = Mass_G, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        axis.title.x = element_blank()) +
  format.plot +
  ylim(0,60) +
  labs(
    title="Figure 1E: Body composition at BW regain",
    y = "Mass (g)",
    fill = "Treatment group",
    color = "Treatment group")


#-------------------------------------------------------------------#
#Stats for Fat, Lean, & BW -> T tests ####
#Is fat, lean, or total mass significantly different between control & WC mice
  #at each of the sable time points?

#T-test for Fat mass (between groups at each sable time point)
ttest_results_Fat <- echoMRI_data  %>%
  group_by(SABLE) %>%
  t_test(Fat ~ GROUP, var.equal = TRUE) %>%   # or var.equal = FALSE if not assumed
  adjust_pvalue(method = "bonferroni") %>%   # optional multiple-comparison correction
  add_significance("p.adj")                  # adds stars based on adjusted p-values
ttest_results_Fat

#T-test for Lean mass (between groups at each sable time point)
ttest_results_Lean <- echoMRI_data  %>%
  group_by(SABLE) %>%
  t_test(Lean ~ GROUP, var.equal = TRUE) %>%   
  adjust_pvalue(method = "bonferroni") %>%   
  add_significance("p.adj")                 
ttest_results_Lean

#T-test for total body mass (between groups at each sable time point)
ttest_results_Weight <- echoMRI_data  %>%
  group_by(SABLE) %>%
  t_test(Weight ~ GROUP, var.equal = TRUE) %>%   
  adjust_pvalue(method = "bonferroni") %>%   
  add_significance("p.adj")                  
ttest_results_Weight

#-------------------------------------------------------------------#
#-------------------------------------------------------------------#
#Graph of BW over the course of the entire experiment











