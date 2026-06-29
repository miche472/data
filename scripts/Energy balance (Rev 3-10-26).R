
#Energy balance
#Revised: 3-6-26

#This script is revised from: "New method energy balance (rev 3-6-26) version B.R"

#Use Carolina's method for demarcating complete day 1 and 2 (rather than my method for identifying 24hr period)
#Since this method includes gaps in recording for some mice, utilize hourly averages rather than direct summation

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
library(lme4)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds")

#---
#Calculate TEE ####
#Prepare data frame for step 1 (step 1 is section below this one)
sable_TEE_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

  
#---
#Step 1a: for each hour, for each ID, find the mean of the EE measurements (in kcal/hr) taken each minute
  mod_TEE_avg_hr <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_TEE_hr = mean(value))
  
#Step 1b: graph hourly TEE from step 1 (units are kcal/hr) --> using "complete days" rather than "LM_complete_day"
  ggplot(mod_TEE_avg_hr, aes(x = hr, y = avg_TEE_hr, fill = GROUP)) +
    geom_bar(stat = "summary", 
             fun = "mean", 
             position = position_dodge(width = 0.9)) +
    geom_errorbar(stat = "summary", 
                  fun.data = mean_se, 
                  position = position_dodge(width = 0.9), 
                  width = 0.3) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 12)) +
    labs(
      title="Average energy expenditure each hour (kcal/hr)",
      x = "Time (Hour)",
      y = "TEE (kcal/hour)",
      fill = "Treatment group") +
    facet_wrap(~SABLE)
  
#---
#Step 2a: for the 24hrs of measurement, sum the average EE from each hour for each mouse
mod_TEE_day <- mod_TEE_avg_hr %>%
    ungroup() %>%
    group_by(SABLE, ID, GROUP) %>%
    summarise(sum_avg_TEE_hr = sum(avg_TEE_hr)) #gives a daily average TEE per mouse
  
#Step 2b: Graph daily TEE from step 2 (raw/unadjusted values)
#use df = mod_TEE_day
ggplot(mod_TEE_day, aes(x = SABLE, y = sum_avg_TEE_hr, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Daily TEE in kcal/day (sum of hourly average)",
    x = "Time point",
    y = "TEE (kcal/day)",
    fill = "Treatment group")

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

#Combine echoMRI_data with mod_TEE_day 
Echo_mod_TEE_day <- mod_TEE_day    %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))

#MLR for TEE
#Build multiple linear regression model for TEE (use df = mod_TEE_day)
model_mod_TEE_day_lean <- lmer(sum_avg_TEE_hr ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_mod_TEE_day)
summary(model_mod_TEE_day_lean)

# Calculate estimated marginal means 
emm_mod_TEE_day_lean <- emmeans(model_mod_TEE_day_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_mod_TEE_day_lean_df <- as.data.frame(emm_mod_TEE_day_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_mod_TEE_day_lean <- contrast(emm_mod_TEE_day_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_mod_TEE_day_lean_df <- as.data.frame(contrasts_by_group_mod_TEE_day_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_mod_TEE_day_lean <- contrast(emm_mod_TEE_day_lean, method = "pairwise", by = "SABLE")
contrasts_SABLE_mod_TEE_day_lean_df <- as.data.frame(contrasts_by_SABLE_mod_TEE_day_lean)

#Conclusions from MLR for TEE adjusted for lean mass (df=Echo_mod_TEE_day)
#As expected: At BW maintenance and BW regain,control mice have significantly higher TEE compared to weight cycled mice
#As expected: for weight cycled mice, TEE at baseline and at peak obesity are both significantly higher than TEE at BW loss
#As expected: TEE is significantly lower at BW maintenance compared to baseline
#***This is generally a pretty good method for TEE in terms of supporting my hypothesis that weight cycling ####
    #leads to a sustained reduction in TEE

#Graph MLR for TEE adjusted for lean mass --> Bar plot
barplot_emm_mod_TEE_day_lean <- emm_mod_TEE_day_lean_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Daily TEE adjusted for lean mass (df=Echo_mod_TEE_day)",
    y = "Daily TEE adjusted for lean mass",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_mod_TEE_day_lean

#formatted differently (just aesthetic differences)
barplot_emm_mod_TEE_day_lean <- emm_mod_TEE_day_lean_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
        axis.title = element_text(face = "bold")) +
  labs(
    title = "Total energy expenditure (kcal/day)",
    y = "TEE (kcal/day) adjusted for lean mass",
    x = "Time point",
    fill = "Treatment group")
barplot_emm_mod_TEE_day_lean


#---
#Step 3a: Break down TEE by light and dark
mod_light_TEE_avg_hr <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, hr, lights) %>%
  arrange(DateTime) %>%
  summarise(avg_TEE_hr = mean(value))

#sum for light and dark period
mod_light_TEE_day <- mod_light_TEE_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(sum_avg_TEE_light = sum(avg_TEE_hr))

#Step 3b: graph TEE during lights on and off
#Light cycle on x-axis
ggplot(mod_light_TEE_day, aes(x = lights, y = sum_avg_TEE_light, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~SABLE) +
  labs(
    title="TEE during light and dark period (sum of hourly average)",
    x = "Light cycle",
    y = "TEE (kcal/period)",
    fill = "Treatment group")

#Step 3c: Break down TEE by light and dark, but use kcal/hr as the metric for comparing light vs dark
#this is necessary because the duration of light and dark are not equal
#This is from above --> get 
mod_light_period_hrly_avg <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  arrange(DateTime) %>%
  summarise(avg_hrly_TEE = mean(value))

#Graph 3c
ggplot(mod_light_period_hrly_avg, aes(x = lights, y = avg_hrly_TEE, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~SABLE) +
  labs(
    title="Average hourly EE during light and dark period",
    x = "Light cycle",
    y = "Hourly EE (kcal/hr)",
    fill = "Treatment group")


#---
#Calculate FI ####
#Approach using Sable measurements --> too many measurements seemed to be incorrect. I don't know that Sable measurements of FI are reliable
sable_FI_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("FoodA_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Find average intake during hr 0, 1, 2, 3, 4, etc. for each ID at each SABLE time point (next step is to 
#Add these values to get 1 day of intake for each mouse at each time point
#Remove Sable FI measurements for weight cycled mice during BW loss and BW maintenance
mod_sable_FI <- sable_FI_modCS %>%
  #filter(!(ID==3707 & SABLE== "BW loss")) %>%
  filter(!(ID==3707)) %>%
  select(-KCAL_PER_GR) %>% # this is wrong, so remove and replace with KCAL_G
  ungroup() %>%
  group_by(ID, SABLE, GROUP, complete_days) %>% #should I group by complete days here?
  arrange(DateTime) %>%       # make sure data is ordered
  mutate(intake = lag(value) - value) %>%
  drop_na() %>% # change in food from one minute to the next
  mutate(
    intake = if_else(intake < 0, 0, intake),  #If mouse doesn't eat between min. x and min. x+1, intake=0
    intake = if_else((intake >2), 0, intake)) %>% #removes values that are illogically high for 1 minute
  ungroup () %>%
  group_by(ID, SABLE, hr, GROUP, complete_days) %>%
  summarise(hr_eaten_gr = sum(intake)) %>%
  ungroup () %>% 
group_by(ID, SABLE, hr, GROUP) %>% #get average intake for each hour across complete days 1 and 2
  summarise(avg_hr_eaten_gr = mean(hr_eaten_gr)) %>%
  ungroup () %>% 
group_by(ID, SABLE, GROUP) %>% #calculate sum of average intake during each of the 24hrs
  summarise(Day_FI_gr = sum(avg_hr_eaten_gr)) %>%
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
  mutate(Day_FI_kcal = Day_FI_gr*KCAL_G) %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice



#Manual FI (first create a df for all mice, then process this to filter for only weight cycled mice at BW loss/maintenance)

#Create df "FI_manual" -> Process FI.csv file (i.e. manual measurements of FI)
FI_manual <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  #filter(corrected_intake_gr < 20 & corrected_intake_gr > 0) %>% #removes 1-29-25 measurements 
  filter(INTAKE_GR <20) %>%
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  #filter(!is.na(INTAKE_GR)) %>% 
  #mutate(INTAKE_GR = replace_na(INTAKE_GR, 0),) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
      sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance", "BW regain"))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  ungroup() 

#Calculate manual FI for all mice
FI_manual_ALL <- FI_manual %>%
  mutate(ID = factor(ID)) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  mutate(INTAKE_kcal = (INTAKE_GR*3.82)) %>%
  #summarise(Day_FI_kcal = mean(INTAKE_kcal)) 
  summarise(Sable_day_FI_kcal= (sum(INTAKE_kcal))/4)
  
#Graph of manual FI --> all mice
ggplot(FI_manual_ALL, aes(x = SABLE, y = Sable_day_FI_kcal, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.9)) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.9),
    width = 0.3) +
  # individual mouse points
  geom_point(
    color = "black",
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
    size = 2)+
  # mouse ID labels
  geom_text_repel(
    aes(label = ID, color = "black"),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9),
    size = 3,
    show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Manually measured FI in Sable (all mice)",
    x = "Time point",
    y = "Daily FI (kcal)",
    fill = "Treatment group")


#FI manual just for weight cycled mice
FI_manual_cycled <- FI_manual %>%
  mutate(ID = factor(ID)) %>%
  filter(GROUP == "Weight cycled") %>%
  filter(SABLE %in% c("BW loss", "BW maintenance")) %>% 
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  #summarise(avg_corrected_intake_kcal = mean(corrected_intake_kcal)) %>% #In this case, this method is fine, but to be consistent I will use the other method
  mutate(INTAKE_kcal = (INTAKE_GR*3.82)) %>%
  summarise(Day_FI_kcal = mean(INTAKE_kcal)) 
  
#Join Sable measured FI and manually measured FI
Join_mod_sable_FI <- mod_sable_FI %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  select(SABLE, GROUP, ID, Day_FI_kcal)

FI_for_balance <- bind_rows(FI_manual_cycled, Join_mod_sable_FI)  %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP)


#---Energy balance (original/less accurate method) ####
#EE is calculated in this script. For EE, use df= mod_TEE_day where daily EE is sum_avg_TEE_hr
#FI is calculated in this script. For FI, use df = FI_for_balance where daily FI is Day_FI_kcal

#Join the FI and the EE data frames; create new variable called E_balance_kcal
E_balance_kcal <- mod_TEE_day  %>%
  left_join(
    FI_for_balance %>% 
      select(ID, SABLE, GROUP, Day_FI_kcal),
    by = c("ID", "SABLE", "GROUP")) %>%
  mutate(balance_kcal= Day_FI_kcal - sum_avg_TEE_hr)

#---Energy balance (alternative/more accurate method) ####
#EE is calculated in this script. For EE, use df= mod_TEE_day where daily EE is sum_avg_TEE_hr
#FI is calculated in this script. For FI, use df = FI_manual_ALL where daily FI is Sable_day_FI_kcal
#Alternative: join manual FI in Sable and EE data frames
Alt_E_balance_kcal <- mod_TEE_day  %>%
  left_join(
    FI_manual_ALL %>% 
      select(ID, SABLE, GROUP, Sable_day_FI_kcal),
    by = c("ID", "SABLE", "GROUP")) %>%
  mutate(balance_kcal= Sable_day_FI_kcal - sum_avg_TEE_hr)

#Graph energy balance (raw values) --> alternative method (i.e. manually measured FI)
ggplot(Alt_E_balance_kcal, aes(x = SABLE, y = balance_kcal, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  # individual mouse points
  geom_point(
    color = "black",
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
    size = 2)+
  # mouse ID labels
  geom_text_repel(
    aes(label = ID, color = "black"),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9),
    size = 3,
    show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Alternative method: Daily energy balance (Manual FI - TEE)",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group")

#Graph Energy balance (raw values) --> original method (i.e. FI recorded by sable)
ggplot(E_balance_kcal, aes(x = SABLE, y = balance_kcal, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Daily energy balance (FI - TEE)",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group")


#Conduct MLR on energy balance (original method/less accurate method)
  #Multiple linear regression to determine if daily energy balance is significantly different
  #between time points and control/weight cycled mice

#Build multiple linear regression model for energy balance (use df = E_balance_kcal )
model_balance_kcal <- lmer(balance_kcal ~ SABLE * GROUP + (1 | ID), data = E_balance_kcal)
summary(model_balance_kcal)

# Calculate estimated marginal means 
emm_balance_kcal <- emmeans(model_balance_kcal, ~ SABLE * GROUP, cov.reduce = mean)
emm_balance_kcal_df <- as.data.frame(emm_balance_kcal)

# Pairwise contrasts within each GROUP
contrasts_by_group_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "GROUP")
contrasts_by_group_balance_kcal_df <- as.data.frame(contrasts_by_group_balance_kcal)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "SABLE")
contrasts_SABLE_balance_kcal_df <- as.data.frame(contrasts_by_SABLE_balance_kcal)

#Graph estimated marginal means
barplot_emm_emm_balance_kcal <- emm_balance_kcal_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Predicted average daily energy balance (kcal)",
    y = "Energy balance (FI-EE, kcal)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_emm_balance_kcal

#Conduct MLR on energy balance (alternative method/more accurate method)
#Multiple linear regression to determine if daily energy balance is significantly different
#between time points and control/weight cycled mice

#Build multiple linear regression model for energy balance (use df = Alt_E_balance_kcal)
model_Alt_E_balance_kcal <- lmer(balance_kcal ~ SABLE * GROUP + (1 | ID), data = Alt_E_balance_kcal)
summary(model_Alt_E_balance_kcal)

# Calculate estimated marginal means 
emm_Alt_E_balance_kcal <- emmeans(model_Alt_E_balance_kcal, ~ SABLE * GROUP, cov.reduce = mean)
emm_Alt_E_balance_kcal_df <- as.data.frame(emm_Alt_E_balance_kcal)

# Pairwise contrasts within each GROUP
contrasts_by_group_Alt_E_balance_kcal <- contrast(emm_Alt_E_balance_kcal, method = "pairwise", by = "GROUP")
contrasts_by_group_Alt_E_balance_kcal_df <- as.data.frame(contrasts_by_group_Alt_E_balance_kcal)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_Alt_E_balance_kcal <- contrast(emm_Alt_E_balance_kcal, method = "pairwise", by = "SABLE")
contrasts_SABLE_Alt_E_balance_kcal_df <- as.data.frame(contrasts_by_SABLE_Alt_E_balance_kcal)

#Graph estimated marginal means
barplot_emm_Alt_E_balance_kcal <- emm_Alt_E_balance_kcal_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Alternative method: Predicted daily energy balance (kcal)",
    y = "Energy balance (FI-EE, kcal)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_Alt_E_balance_kcal




#Add change in BW over the course of Sable --> quality control

  #consensus paper says to do this when looking at energy balance for quality control
  # if mice have positive energy balance during recording they should gain weight.
  #If this isn't true, there is likely an issue with the recording for a cage and
  #that cage should be removed from analysis
  #consensus method advises to use first and last BW measurement of recording for each mouse

#Prepare data frame for step 1 (step 1 is section below this one)
sable_BW_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("BodyMass_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Check if there are any BWs that are zero or less than zero

#Find first and last BW measurement for each mouse during each Sable time point


Find_sable_BW_modCS <- sable_BW_modCS %>%
  filter(value > 0) %>%
  filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  mutate(
    start_bw = first(value),
    first_DateTime = first(DateTime),
    last_DateTime = last(DateTime))




Find_sable_BW_modCS <- sable_BW_modCS %>%
  #filter(value > 0) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  mutate(
    first_bw = first(value),
    last_bw  = last(value),
    bw_change = last_bw - first_bw)

Find_sable_BW_modCS <- sable_BW_modCS %>%
  #filter(value > 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop")

Find_sable_BW_modCS2 <- sable_BW_modCS %>%
  filter(value > 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop")

Find_sable_BW_modCS2 <- sable_BW_modCS %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  #filter(value == 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop") %>%
  filter(bw_change == 0)

#Create df "FI_manual" -> Process FI.csv file (i.e. manual measurements of FI)
BW_manual <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
    sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
    TRUE ~ NA_character_)) %>% 
  #filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance", "BW regain"))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  ungroup() 
