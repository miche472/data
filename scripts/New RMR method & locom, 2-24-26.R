#RMR --> calculated using additive method
#Locomotion --> calculated using additive method

#RMR and locomotion calculated in this script use NEAT_LM, 2-23-26 as starting point

#Started:2-24-26
#Revised: 2-24-26

#Left off at bottom of script. Need to figure out why TEE ~ RMR and NEAT > TEE using this method ####

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

#### Calculate RMR for All SABLE time points for all NZO mice ####
#Use AllMeters to determine the minutes during which a mouse was moving vs staying still

filter_EE <-sable_dwn %>%
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
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>%
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#------- Get only the AllMeters parameter --> value column has AllMeters data
filter_loc <-sable_dwn %>%
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
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>%
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  #For NZO_Figure7b-RRM_correctedbyLean I used both complete days 1 and 2...complete_days %in% c(1,2)
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725, 3711)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#---
#In df filter_loc use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE use mutate to make a column called kcal_hr_ using the value column data

filter_locom <- filter_loc %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_energy <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_locom_energy <- filter_locom %>%
  left_join(
    filter_energy %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Calculate daily RMR for each mouse ID at each SABLE time point
#Filter for times when mouse didn't move. During these moments TEE is entirely resting metabolic rate
  # (i.e.) physical activity is not contributing to TEE when mouse isn't moving
ID_RMR <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse moved
  mutate(RMR_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop")

#Calculate average daily RMR within each GROUP at each SABLE time point
GROUP_RMR <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>%
  mutate(RMR_per_min= (Kcal_Hr/60)) %>%
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_RMR = mean(RMR_teske), 
            SD_GROUP_RMR = sd(RMR_teske),
            SE_GROUP_RMR = sd(RMR_teske)/sqrt(n()))

#### Extrapolate RMR from EE during the minutes when mice aren't moving ####
# Issue for both methods: RMR seems too high given what TOTAL energy expenditure is (using other methods to calculate TEE)


#Approach B: Find avg hourly RMR for Control and Weight cycled mice during dark and light period
  #Similar method used by Banks et al. 2025 in "Consensus" paper
  #They calculated average hourly TEE during dark and light period
RMR_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(RMR_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, GROUP, lights) %>% n 
  summarise(RMR_GROUP = mean(RMR_kcal_hr), .groups="drop") %>& 

#My original version --> Extrapolate RMR from EE during the minutes when mice aren't moving
  #Calculate avg EE when each mouse isn't moving during dark period and during light period. 
  #Extrapolate to get EE when mouse isn't moving for the entire day by multiplying by the number of hrs
  #of light and dark. Then add together RMR during light and dark period
  RMR_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(RMR_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
  mutate(RMR_per_photo = photo_length_hrs*RMR_kcal_hr) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
#To get total daily RMR for each ID
  summarise(Daily_RMR_ID = sum(RMR_per_photo), .groups="drop") %>%
#To get total daily RMR for each GROUP
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(Daily_RMR_GROUP = mean(Daily_RMR_ID), .groups="drop")
  

#Try using my original method for calculating RMR (above) to calculate TEE
  #Is TEE calculated with this method unreasonably high?
TEE_extrap <- filter_locom_energy %>%
    ungroup() %>%
    arrange(DateTime) %>%     # make sure rows are in time order
    group_by(SABLE, ID, GROUP) %>%
    mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
    #filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
    ungroup() %>%
    group_by(SABLE, ID, GROUP, lights) %>%
    summarise(TEE_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
    group_by(SABLE, ID, GROUP, lights) %>%
    mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
    mutate(TEE_per_photo = photo_length_hrs*TEE_kcal_hr) %>%
    ungroup() %>%
    group_by(ID, SABLE, GROUP) %>%
    #To get total daily RMR for each ID
    summarise(Daily_TEE_ID = sum(TEE_per_photo), .groups="drop") %>%
    #To get total daily RMR for each GROUP
    ungroup() %>%
    group_by(SABLE, GROUP) %>%
    summarise(Daily_TEE_GROUP = mean(Daily_TEE_ID), .groups="drop")
#Daily TEE is essentially the same as daily RMR...this isn't plausible
#Try using this method to calculate NEAT

#Use this method to calculate NEAT
NEAT_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>% #only keep data from minutes when the mouse moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(NEAT_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
  mutate(NEAT_per_photo = photo_length_hrs*NEAT_kcal_hr) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  #To get total daily RMR for each ID
  summarise(Daily_NEAT_ID = sum(NEAT_per_photo), .groups="drop") %>%
  #To get total daily RMR for each GROUP
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(Daily_NEAT_GROUP = mean(Daily_NEAT_ID), .groups="drop")
#Values for NEAT are larger than values for TEE...there is an issue with this method of calculation I think
#Need to think more about this

