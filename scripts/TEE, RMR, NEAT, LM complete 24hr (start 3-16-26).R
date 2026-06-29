#TEE, NEAT, and RMR using 24hr perfect day and hourly averaging method. 
#both of these approaches are intended to reduce the impact of outliers and 
#prevent the presence of missing observations from skewing data 

# #Libraries####
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

#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))}

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#Step 1: Establish a 24hr day for each mouse at each time point ####
#Took code from New complete day method (Rev. 3-16-26)

##------- Get only the All_meters parameter --> value column has All_meters data
filter_loc1 <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) 

#Get complete day 1 and 2 by explicitly definining the zt_time and SABLE_DAY
filter_loc2 <- filter_loc1 %>%
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
    #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
    #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
    #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
    #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_loc3 <- filter_loc2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2))

#Check number of mice in filter_loc3
filter_loc3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#------- Get only the kcal_hr parameter --> value column has kcal_hr data
filter_EE1 <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
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
  group_by(ID) 

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
filter_EE2 <- filter_EE1 %>% 
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
    #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
    #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
    #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
    #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_EE3 <- filter_EE2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2)) #remove observations not from complete day 1 or 2

#Check number of mice in filter_EE3
filter_EE3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#---
#In df filter_loc3 use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE3 use mutate to make a column called kcal_hr_ using the value column data

filter_locom3 <- filter_loc3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_energy3 <- filter_EE3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_locom_energy3 <- filter_locom3 %>%
  left_join(
    filter_energy3 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

# Step 2: Use Teske method to delineate RMR & NEAT based on whether the mouse moved or not between minutes ####
#Took code from script: New complete day method (Rev. 3-16-26).R

#### Identify when TEE is RMR for each ID at each SABLE time point ####
#Approach: Filter for minutes when mouse didn't move. During these minutes TEE is entirely RMR
# (i.e.) physical activity is not contributing to TEE when mouse isn't moving

#RMR across one complete day (24hrs) for each mouse at each SABLE time point
ID_RMR3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) #only keep data from minutes when the mouse moved

#### Identify when TEE is NEAT for each ID at each SABLE time point ####
#Approach: Filter for minutes when mouse moved. During these minutes TEE is NEAT
# (i.e.) mouse us moving so TEE is due to NEAT (activity)

#NEAT across one complete day (24hrs) for each mouse at each SABLE time point
ID_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) #only keep data from minutes when the mouse moved
  
#Step 3: Calculate daily RMR: calculate avg hourly RMR (kcal/hr) and sum RMR during each hr to get daily RMR (kcal/day) ####
#Took code from script: Energy balance (Rev 3-13-26).R

# RMR_hr: For each mouse and SABLE time point, find average RMR for each hr during a 24hr period
mod_RMR_avg_hr <- ID_RMR3 %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_RMR_hr = mean(Kcal_Hr))

#Graph RMR_hr: RMR for each hour during each SABLE time point (units are kcal/hr)
ggplot(mod_RMR_avg_hr, aes(x = hr, y = avg_RMR_hr, fill = GROUP)) +
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
    title="Average RMR each hour (kcal/hr)",
    x = "Time (Hour)",
    y = "RMR (kcal/hour)",
    fill = "Treatment group") +
  facet_wrap(~SABLE)

# RMR_daily: For the 24hrs of measurement, sum the average RMR from each hour for each mouse
mod_RMR_day <- mod_RMR_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  summarise(sum_avg_RMR_hr = sum(avg_RMR_hr)) #gives a daily average RMR per mouse

#Graph RMR_daily: Daily RMR (kcal/day) --> raw values
#use df = mod_TEE_day
ggplot(mod_RMR_day, aes(x = SABLE, y = sum_avg_RMR_hr, fill = GROUP)) +
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
    title="Daily RMR in kcal/day (sum of hourly average)",
    x = "Time point",
    y = "RMR (kcal/day)",
    fill = "Treatment group")

# t-test comparing RMR for control and weight cycled at BW regain 
#(RMR for each mouse was calculated as the sum of hourly average RMR)
#method for complete day was LM's method for perfect 24hr day
mod_RMR_day_regain <- mod_RMR_day %>% 
  filter(SABLE=="BW regain") %>%
  ungroup() %>% 
  group_by(GROUP)

#Summary statistics
mod_RMR_day_regain %>%
  group_by(GROUP) %>%
  get_summary_stats(sum_avg_RMR_hr, type = "mean_sd")

#T test
res <- t.test(sum_avg_RMR_hr ~ GROUP, data = mod_RMR_day_regain)
res
#For t-test of control vs weight cycled at BW regain, p = 0.070 --> almost significant

#---
#Step 4: Calculate daily NEAT: calculate avg NEAT for each of 24hrs (kcal/hr); sum NEAT from each hr to get daily NEAT (kcal/day) for each ID ####
#Took code from script: Energy balance (Rev 3-13-26).R

# NEAT_hr: For each mouse and SABLE time point, find average NEAT for each hr during a 24hr period
mod_NEAT_avg_hr <- ID_NEAT3 %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_NEAT_hr = mean(Kcal_Hr))

#Graph NEAT_hr: NEAT for each hour during each SABLE time point (units are kcal/hr)
ggplot(mod_NEAT_avg_hr, aes(x = hr, y = avg_NEAT_hr, fill = GROUP)) +
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
    title="NEAT each hour (kcal/hr)",
    x = "Time (Hour)",
    y = "NEAT (kcal/hour)",
    fill = "Treatment group") +
  facet_wrap(~SABLE)

# NEAT_daily: For the 24hrs of measurement, sum the average NEAT from each hour for each mouse
mod_NEAT_day <- mod_NEAT_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  summarise(sum_avg_NEAT_hr = sum(avg_NEAT_hr)) #gives daily NEAT per mouse

#Graph NEAT_daily: Daily NEAT (kcal/day) --> raw values
#use df = mod_NEAT_day
ggplot(mod_NEAT_day, aes(x = SABLE, y = sum_avg_NEAT_hr, fill = GROUP)) +
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
    title="Daily NEAT in kcal/day (sum of average from each hr)",
    x = "Time point",
    y = "NEAT (kcal/day)",
    fill = "Treatment group")

# T-test comparing NEAT for control and weight cycled at BW regain 
  #NEAT for each mouse was calculated as the sum of hourly average NEAT
  #method for complete day was LM's method for perfect 24hr day
mod_NEAT_day_regain <- mod_NEAT_day %>% #Filter for just BW regain timepoint
  filter(SABLE=="BW regain") %>%
  ungroup() %>% 
  group_by(GROUP)

mod_NEAT_day_regain %>% #Summary statistics
  group_by(GROUP) %>%
  get_summary_stats(sum_avg_NEAT_hr, type = "mean_sd")

#T test
res <- t.test(sum_avg_NEAT_hr ~ GROUP, data = mod_NEAT_day_regain)
res
# NEAT control vs NEAT weight cycled at BW regain, p = 0.0681 --> almost significant

#---
#### Step 5: Calculate TEE ####
#For TEE, skip Step 2 (TEE encompasses both when the mouse moved and when it didn't move...NEAT+RMR)
#So, use df= filter_locom_energy3

# TEE_hr: For each mouse and SABLE time point, find average TEE for each hr during a 24hr period
mod_TEE_avg_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_TEE_hr = mean(Kcal_Hr))

#Graph TEE_hr: TEE for each hour during each SABLE time point (units are kcal/hr)
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
    title="TEE each hour (kcal/hr)",
    x = "Time (Hour)",
    y = "TEE (kcal/hour)",
    fill = "Treatment group") +
  facet_wrap(~SABLE)

# TEE_daily: For the 24hrs of measurement, sum the avg TEE from each hour for each ID (avg the 60 observations taken each hr)
mod_TEE_day <- mod_TEE_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  summarise(sum_avg_TEE_hr = sum(avg_TEE_hr)) #gives daily TEE per mouse

#Graph TEE_daily: Daily TEE (kcal/day) --> raw values
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
    title="Daily TEE in kcal/day (sum of average from each hr)",
    x = "Time point",
    y = "TEE (kcal/day)",
    fill = "Treatment group")

# T-test comparing TEE for control and weight cycled at BW regain 
#TEE for each mouse was calculated as the sum of hourly average TEE
#method for complete day was LM's method for perfect 24hr day
mod_TEE_day_regain <- mod_TEE_day %>% #Filter for just BW regain timepoint
  filter(SABLE=="BW regain") %>%
  ungroup() %>% 
  group_by(GROUP)

mod_TEE_day_regain %>% #Summary statistics
  group_by(GROUP) %>%
  get_summary_stats(sum_avg_TEE_hr, type = "mean_sd")

#T test
res <- t.test(sum_avg_TEE_hr ~ GROUP, data = mod_TEE_day_regain)
res
# TEE control vs TEE weight cycled at BW regain, p = 0.1053 --> 
    #closer to significant than linear mixed models using different combinations of 
    #complete days & calculation methods for TEE have been

#--
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

#Combine echoMRI_data, mod_RMR_day, mod_NEAT_day, and mod_TEE_day
Echo_mod_EE_day <- mod_TEE_day %>%
  left_join(echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE")) %>%
  left_join(mod_NEAT_day %>% 
      select(ID, SABLE, sum_avg_NEAT_hr),
    by = c("ID", "SABLE")) %>%
  left_join(mod_RMR_day %>% 
              select(ID, SABLE, sum_avg_RMR_hr),
            by = c("ID", "SABLE")) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  mutate(RMR_NEAT = sum_avg_RMR_hr + sum_avg_NEAT_hr,
        diff= RMR_NEAT - sum_avg_TEE_hr)
#Issue: the averaging method is not compatible with Teske's method for distinguishing between NEAT vs RMR

#Trying to solve issue identified in step 5
#Step 6: look at what percent of minutes mice moved vs didn't move
#Graph 
#RMR across one complete day (24hrs) for each mouse at each SABLE time point
ID_RMR3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) #only keep data from minutes when the mouse moved

#### Identify when TEE is NEAT for each ID at each SABLE time point ####
#Approach: Filter for minutes when mouse moved. During these minutes TEE is NEAT
# (i.e.) mouse us moving so TEE is due to NEAT (activity)

#NEAT across one complete day (24hrs) for each mouse at each SABLE time point
ID_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) #only keep data from minutes when the mouse moved
