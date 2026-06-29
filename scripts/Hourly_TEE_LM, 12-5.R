#Hourly EE and light/dark EE (unadjusted)
#LM, 12-2-25

# * indicates that I am still working on the section
#I need to figure out how to calculate adjusted values for hourly EE (adjusted for lean)

#Libraries####
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

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

# This is a modified version of NZO_FIgure7b -RMR_CorrectedbyLean 10-30 (accessed on 12-2-25)
  #In the source script the first step creates a data frame called sable_RMR_data
  # I named it sable_hrly_data to avoid confusion between the scripts
# I used the RMR script as a starting point to get minute by minute data for EE ("parameter"= kcal_hr)

sable_hourly_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  
  # remove dead mice, keep both complete days, remove mice with cage issues
 # filter(is_complete_day == 1, complete_days == 2) %>% #In NZO_Figure7b -RMR_CorrectedbyLean 10-30 I had complete_days %in% c(1,2)
  filter(!ID %in% c(3715, 3712)) %>%
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  ungroup() %>%
  group_by(ID, SABLE) %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47"))

sable_hourly_data <- sable_hourly_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))


#### Options 1 and 2 relate to how I reattach (or just retain) columns while calculating average EE during each hr ####
    #Option 2 is the one that I chose to go with, but I am keeping Option 1 in the script because it shows
    # the logic for how I worked through this calculation 
#Option 1: 
#Average EE during each hour (0-23) for each mouse ID for each Sable time point
#EE_hourly_1 <- sable_hourly_data %>%
  #mutate(hour = hour(DateTime)) %>%
  #group_by(ID, SABLE, hour) %>%
  #summarise(EE_hour_avg = mean(value, na.rm = TRUE),.groups = "drop")

#Units are in kcal per hr --> want to check if kcal per day is reasonable
#Hourly_EE_per_day_1 <-EE_hourly_1 %>%
#mutate(kcal_day=EE_hour_avg*24) %>%

# reattach GROUP and DRUG
#EE_hourly_1<- EE_hourly_1 %>%
  #mutate(GROUP = case_when(ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      #ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    #DRUG = case_when(ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      #ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  #mutate(SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Option 2:**this is the one that I chose instead of option 1
#Alternative to reattaching GROUP and DRUG (also retains lights and hr from sable_hourly_data)

#This sanity check confirms that I have the correct data structure to use first() in the section below
sanity_check <- sable_hourly_data %>%
  mutate(hour = hour(DateTime)) %>%
  group_by(ID, SABLE, hour) %>%
  summarise(
    lights_n = n_distinct(lights),
    group_n  = n_distinct(GROUP),
    drug_n   = n_distinct(DRUG),
    .groups = "drop") %>%
  filter(lights_n > 1 | group_n > 1 | drug_n > 1)
#The df that results from the sanity check returns zero rows meaning that I am safe to use first()

EE_hourly <- sable_hourly_data %>%
  mutate(hour = hour(DateTime)) %>%
  group_by(ID, SABLE, hour) %>%
  summarise(
    EE_hour_avg = mean(value, na.rm = TRUE),
    lights      = first(lights),
    GROUP       = first(GROUP),
    DRUG        = first(DRUG),
    hr          =first(hr),
    .groups = "drop") %>%
  mutate(
    SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance","BW regain")))

lights_on_SABLE <- EE_hourly %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  filter(lights=="on")

#Format plot
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#3498DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

#####----- Barplot - NZO hourly energy expenditure (measured values from EE_hourly) ####
lights_on_SABLE_plot <- lights_on_SABLE %>%
  filter(GROUP=="Weight cycled") %>%
  ggplot(aes(x = hour, y = EE_hour_avg, fill = GROUP)) + 
  stat_summary( # mean bars
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  stat_summary( # error bars (mean ± SE)
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(  # individual data points
    aes(color = GROUP), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  labs(title = "Hourly energy expenditure during diurnal period", 
       y = "TEE (kcal/hour)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color="black", size=14),
    plot.title = element_text(hjust = 0.5, face = "bold")) +
  facet_wrap(~SABLE)
lights_on_SABLE_plot

#### Do hourly TEE adjusted for lean mass ####
#Need to reattach EchoMRI data
#echoMRI
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join body composition info (echoMRI_data) into hourly energy expenditure data frame (EE_hourly)
EE_hourly_adj <- EE_hourly %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean, Fat, Weight, adiposity_index),
    by = c("ID", "SABLE"))



#### *Haven't finished this/left off here on 12-2-25 at 5:45pm.  ####
      #How should I account for the hour element?
#### Create mixed model that adjusts for lean mass ####
  #### Build linear mixed model for RMR (adjusted by lean) ####
  model_EE_hr_lean <- lmer(EE_hour_avg ~ SABLE * GROUP * hr + Lean + (1 | ID), data = EE_hourly_adj)
  summary(model_EE_hr_lean)

####------------------------------------
#### EE during the Light and Dark period for each mouse during each Sable time point ####
#For each mouse at each Sable time point, average the minute by minute TEE (kcal_hr)
  #which are during the dark period. If this doesn't work I can average the hourly values
   
EE_light_dark <- sable_hourly_data %>%
    #mutate(hour = hour(DateTime)) %>%
    group_by(ID, SABLE, lights) %>%
    summarise(
      EE_lightcycle = mean(value, na.rm = TRUE),
      GROUP       = first(GROUP),
      DRUG        = first(DRUG),
      #complete_days = first(complete_days),
      .groups = "drop"
    ) %>%
    mutate(
      SABLE = factor(SABLE, 
                     levels = c("Baseline", 
                                "Peak obesity", 
                                "BW loss", 
                                "BW maintenance", 
                                "BW regain")))

#####----- Barplot - NZO hourly energy expenditure during light/dark (unadjusted) ####
  #EE_light_dark_plot <- EE_light_dark %>%
    #filter(GROUP=="Weight cycled") %>%
    ggplot(EE_light_dark, aes(x = lights, y = EE_lightcycle, fill = GROUP)) + 
    stat_summary( # mean bars
      fun = mean,
      geom = "col",
      position = position_dodge(width = 0.8),
      color = "black", width = 0.7, alpha = 0.7) +
    stat_summary( # error bars (mean ± SE)
      fun.data = mean_se,
      geom = "errorbar",
      position = position_dodge(width = 0.8),
      width = 0.3) +
    geom_point(  # individual data points
      aes(color = GROUP), 
      position = position_dodge(width = 0.8),
      alpha = 0.7, size = 2) + 
    scaleFill + 
    scaleColor +
    theme_minimal() +
    labs(title = "Hourly energy expenditure during diurnal period", 
         y = "TEE (kcal/hr)", 
         x= "Time point", 
         color = "Light cycle",
         fill = "Light cycle") +
    format.plot +
    theme(
      legend.position = "top",
      axis.ticks.y = element_line(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(color="black", size=14),
      plot.title = element_text(hjust = 0.5, face = "bold")) +
    facet_wrap(~SABLE)
  #EE_light_dark_plot
  
#### *Haven't finished this/left off here on 12-2-25  ####
  #How should I account for lights, SABLE, and GROUP...interaction between all...what do I want to compare?
#### Create mixed model that adjusts for lean mass ####
#### Build linear mixed model for EE during light/dark (adjusted by lean) ####
  model_EE_hr_lean <- lmer(EE_hour_avg ~ SABLE * GROUP * hr + Lean + (1 | ID), data = EE_hourly_adj)
  summary(model_EE_hr_lean)
