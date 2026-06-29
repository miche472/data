
# Began with a section of code from New RMR method & locom, (2-26-26 issues).R
# Started: 2-26-26
#Revised: 3-17-26
#On 3-17-26 I removed the sections with FI and E. balance that were in previous version

#Goals: calculate TEE, RMR, and NEAT by summing one complete day for each ID at each SABLE point
  #criteria for selecting a complete day was minimizing the number of missing observations 
  #(i.e. minimizing the number of minutes without data)
#Compared TEE and RMR calculated using the NEAT+RMR summation method and two gating methods: 
    #new method is perfect 24hr day (LM developed) and CS's method for gating (method used for prelim)

#Method & reasoning used to select complete days method (my new method) is in script: 
    #"New RMR method & locom, (2-26-26 issues).R"

#Additional goals: calculate meters moved and time spent moving during a complete day
#I used AllMeters for this, but I actually should have used PedMeters to get strictly ambulation


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
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))}

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

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

# Calculate daily RMR for each mouse and for each GROUP at each SABLE time point ####
  #Use minute by minute summation method
  #Approach: Filter for minutes when mouse didn't move. During these minutes TEE is entirely RMR
  # (i.e.) physical activity is not contributing to TEE when mouse isn't moving

#RMR across one complete day for each mouse at each SABLE time point
ID_RMR3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse moved
  mutate(RMR_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop")

#Average daily RMR within each GROUP at each SABLE time point
GROUP_RMR3 <- filter_locom_energy3 %>%
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

#### Calc. daily NEAT for each mouse and for each GROUP at each SABLE time point ####
#use summation method (sum minute values for NEAT to get daily NEAT)

#NEAT across one complete day for each mouse at each SABLE time point
ID_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>% #only keep data from minutes when the mouse moved
  mutate(NEAT_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop")

#Average daily NEAT within each GROUP at each SABLE time point
GROUP_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  mutate(NEAT_per_min= (Kcal_Hr/60)) %>%
  summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_NEAT = mean(NEAT_teske), 
            SD_GROUP_NEAT = sd(NEAT_teske),
            SE_GROUP_NEAT = sd(NEAT_teske)/sqrt(n()))

#Calculate daily TEE for each mouse and for each GROUP at each SABLE time point ####
  #Use summation method (sum minute values for EE to get daily EE)

#TEE across one complete day for each mouse at each SABLE time point
ID_TEE3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(TEE_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(TEE_teske = sum(TEE_per_min), .groups="drop")

#Average daily TEE within each GROUP at each SABLE time point
GROUP_TEE3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(TEE_per_min= (Kcal_Hr/60)) %>%
  summarise(TEE_teske = sum(TEE_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_TEE = mean(TEE_teske), 
            SD_GROUP_TEE = sd(TEE_teske),
            SE_GROUP_TEE = sd(TEE_teske)/sqrt(n()))


#Compare the summation method for calculating TEE & RMR (developed here) to original method ####
  #which used a different method for determining what a complete day was (still used summation for TEE)
  #(Scripts: "NZO_figure7 - TEE corrected by BW (rev. LM)" and NZO_Figure7b-RMR_correctedbyLean_(LM, 12-5))
  #The original script for RMR calculated RMR as the 30min interval with the lowest avg TEE
#Both methods use summation, the difference is in how complete days were determined

#Original method for TEE: calculated avg daily TEE for each ID between two complete days (gated with CS's method) 
  #calculated in script called NZO_Figure7b-RMR_correctedbyLean_(LM, 12-5)
  #code to make df=sable_TEE_adj_RMR which has tee and RMR ...pasted into bottom of this script
#Average daily TEE within each GROUP at each SABLE time point (using tee from CS's original method)
old_GROUP_TEE <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #sable cage tech issues
  #filter(!(ID==3711)) %>% #visually looked like possible outlier
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_TEE = mean(tee), 
            SD_GROUP_TEE = sd(tee),
            SE_GROUP_TEE = sd(tee)/sqrt(n()))

#Original method for RMR: look at avg TEE across every 30min window over approx 48hr period.
  #multiply out the EE during the 30min window with the lowest avg TEE measurement to get RMR for 24hrs
#Average daily RMR within each GROUP at each SABLE time point (original method for tee and RMR) 
#RMR_kcal_day is per ID and is a variable located in df called sable_TEE_adj_RMR
#created in script: NZO_Figure7b-RMR_correctedbyLean_(LM, 12-5)...pasted into bottom of this script
old_GROUP_RMR <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  #filter(!(ID==3711)) %>% #visually looked like possible outlier
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_RMR = mean(RMR_kcal_day), 
            SD_GROUP_RMR = sd(RMR_kcal_day),
            SE_GROUP_RMR = sd(RMR_kcal_day)/sqrt(n()))

##---------------------------------------------####
#Multiple linear regression and bar plots for EE Summation method

#Process/format echoMRI data ####
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
  filter(!(ID == 3711 & Date == as.Date("2025-01-27"))) %>% #repeated echo
  filter(!(ID == 3727 & Date == as.Date("2025-02-07"))) %>% #repeated echo
  mutate(STATUS = factor(STATUS,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  rename(SABLE = STATUS) 

#Check number of mice in echoMRI_data
echoMRI_data %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#Combine echoMRI_data with ID_TEE3
Echo_TEE3 <- ID_TEE3  %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))

#Combine df with TEE (sum) and echo, with RMR (sum) -> df= Echo_TEE_RMR3
Echo_TEE_RMR3 <- Echo_TEE3  %>%
  left_join(
    ID_RMR3 %>% 
      select(ID, SABLE, RMR_teske),
    by = c("ID", "SABLE"))

#Combine echo, TEE, RMR df with NEAT --> df=Echo_TEE_RMR_NEAT3
Echo_TEE_RMR_NEAT3 <- Echo_TEE_RMR3  %>%
  left_join(
    ID_NEAT3 %>% 
      select(ID, SABLE, NEAT_teske),
    by = c("ID", "SABLE"))

#Confirm that data frames have correct number of mice (good, we have n=16 for all SABLEs)
Echo_TEE_RMR3 %>% group_by(SABLE) %>% summarise(n_ID = n_distinct(ID))
Echo_TEE_RMR_NEAT3 %>% group_by(SABLE) %>% summarise(n_ID = n_distinct(ID))

#--Graph raw values for RMR, NEAT, and TEE
#use df = Echo_TEE_RMR_NEAT3
#Graph raw RMR (kcal/day) ####
ggplot(Echo_TEE_RMR_NEAT3, aes(x = SABLE, y = RMR_teske, fill = GROUP)) +
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
    title="RMR (kcal/day)--Raw",
    x = "Time point",
    y = "RMR (kcal/day)",
    fill = "Treatment group")

#Graph raw NEAT (kcal/day) ####
ggplot(Echo_TEE_RMR_NEAT3, aes(x = SABLE, y = NEAT_teske, fill = GROUP)) +
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
    title="NEAT (kcal/day)--Raw",
    x = "Time point",
    y = "NEAT (kcal/day)",
    fill = "Treatment group")

#Graph raw TEE (kcal/day) ####
ggplot(Echo_TEE_RMR_NEAT3, aes(x = SABLE, y = TEE_teske, fill = GROUP)) +
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
    title="TEE (kcal/day)--Raw",
    x = "Time point",
    y = "TEE (kcal/day)",
    fill = "Treatment group")

#--MLR: TEE, NEAT, RMR (direct minute summation method with LM's method to identify 24hr period)-----####
    #Use df = Echo_TEE_RMR3

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

###-----RMR adj. adiposity index (fat/lean) (minute summation method)-----####
#Build multiple linear regression model
model_RMR_AI <- lmer(RMR_teske ~ SABLE * GROUP + adiposity_index + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_AI)

# Calculate estimated marginal means #
emm_RMR_AI <- emmeans(model_RMR_AI, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_AI_df <- as.data.frame(emm_RMR_AI)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_AI <- contrast(emm_RMR_AI, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_AI_df <- as.data.frame(contrasts_by_group_RMR_AI)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_AI <- contrast(emm_RMR_AI, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_AI_df <- as.data.frame(contrasts_by_SABLE_RMR_AI)

# Bar plot - Graph predicted RMR adjusted for Lean mass #
barplot_emm_RMR_AI <- emm_RMR_AI_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar( 
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "RMR adjusted for adiposity index (fat/lean)",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_AI

#--
#For reference: adiposity index graph
ggplot(Echo_TEE_RMR3, aes(x = SABLE, y = adiposity_index, fill = GROUP)) +
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
  labs(
    title="Adiposity index (fat/lean) at iCAL time points",
    x = "Time point",
    y = "Adiposity index (fat/lean)",
    fill = "Treatment group")

#---
#For reference: BW graph (BW taken from EchoMRI, not manual measurements)
ggplot(Echo_TEE_RMR3, aes(x = SABLE, y = Weight, fill = GROUP)) +
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
  labs(
    title="Body weight at iCAL time points",
    x = "Time point",
    y = "Weight (grams)",
    fill = "Treatment group")

#---
#For reference: Lean mass graph
ggplot(Echo_TEE_RMR3, aes(x = SABLE, y = Lean, fill = GROUP)) +
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
  labs(
    title="Lean mass at iCAL time points",
    x = "Time point",
    y = "Lean mass (grams)",
    fill = "Treatment group")

#---
#For reference: Fat mass graph
ggplot(Echo_TEE_RMR3, aes(x = SABLE, y = Fat, fill = GROUP)) +
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
  labs(
    title="Fat mass at iCAL time points",
    x = "Time point",
    y = "Fat mass (grams)",
    fill = "Treatment group")

#---
###-----NEAT adj. for lean mass (minute summation method)-----####
#Build multiple linear regression model
model_NEAT_Lean <- lmer(NEAT_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR_NEAT3)
summary(model_NEAT_Lean)

# Calculate estimated marginal means #
emm_NEAT_Lean<- emmeans(model_NEAT_Lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_Lean_df <- as.data.frame(emm_NEAT_Lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_Lean <- contrast(emm_NEAT_Lean, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_Lean_df <- as.data.frame(contrasts_by_group_NEAT_Lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_Lean <- contrast(emm_NEAT_Lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_NEAT_Lean_df <- as.data.frame(contrasts_by_SABLE_NEAT_Lean)

# Bar plot - Graph predicted NEAT adjusted for Lean mass #
barplot_emm_NEAT_Lean <- emm_NEAT_Lean_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "NEAT (kcal/day) adjusted for lean mass",
    y = "Adjusted NEAT (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_NEAT_Lean
#Results: at BW regain, control vs. weight cycled --> p=0.097
#For weight cycled, NEAT is significantly higher at BW loss compared to BW regain


###-----NEAT adj. for BW (minute summation method)-----####
#Build multiple linear regression model
model_NEAT_BW <- lmer(NEAT_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR_NEAT3)
summary(model_NEAT_BW)

# Calculate estimated marginal means #
emm_NEAT_BW<- emmeans(model_NEAT_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_BW_df <- as.data.frame(emm_NEAT_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_BW <- contrast(emm_NEAT_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_BW_df <- as.data.frame(contrasts_by_group_NEAT_BW)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_BW <- contrast(emm_NEAT_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_NEAT_BW_df <- as.data.frame(contrasts_by_SABLE_NEAT_BW)

#Conclusions: Non-significant difference at BW regain (p=0.228)

# Bar plot - Graph predicted NEAT adjusted for BW #
barplot_emm_NEAT_BW <- emm_NEAT_BW_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "NEAT adjusted for body weight",
    y = "Adjusted NEAT (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_NEAT_BW

###-----RMR adj. for BW (minute summation method)-----#####
#Build multiple linear regression model
model_RMR_BW <- lmer(RMR_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_BW)

# Calculate estimated marginal means #
emm_RMR_BW <- emmeans(model_RMR_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_BW_df <- as.data.frame(emm_RMR_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_BW <- contrast(emm_RMR_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_BW_df <- as.data.frame(contrasts_by_group_RMR_BW)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_BW <- contrast(emm_RMR_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_BW_df <- as.data.frame(contrasts_by_SABLE_RMR_BW)

# Bar plot - Graph predicted RMR adjusted for Lean mass #
barplot_emm_RMR_BW <- emm_RMR_BW_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "RMR adjusted for body weight",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_BW
#Conclusion: Completely counterintuitive result

###-----RMR adj. for lean mass (minute summation method)-----####
#Build multiple linear regression model
model_RMR_lean <- lmer(RMR_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_lean)

# Calculate estimated marginal means #
emm_RMR_lean <- emmeans(model_RMR_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_lean_df <- as.data.frame(emm_RMR_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_lean <- contrast(emm_RMR_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_lean_df <- as.data.frame(contrasts_by_group_RMR_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_lean <- contrast(emm_RMR_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_lean_df <- as.data.frame(contrasts_by_SABLE_RMR_lean)

# Bar plot - Graph predicted RMR adjusted for Lean mass #
barplot_emm_RMR_lean <- emm_RMR_lean_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Resting metabolic rate (RMR) adjusted for lean mass",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_lean
#Conclusion: completely counter intuitive result

###-----TEE adj. lean mass (minute summation method)-----####
#Build multiple linear regression model
model_TEE_lean <- lmer(TEE_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR3)
summary(model_TEE_lean)

# Calculate estimated marginal means #
emm_TEE_lean <- emmeans(model_TEE_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_lean_df <- as.data.frame(emm_TEE_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_lean <- contrast(emm_TEE_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_lean_df <- as.data.frame(contrasts_by_group_TEE_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_lean <- contrast(emm_TEE_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_TEE_lean_df <- as.data.frame(contrasts_by_SABLE_TEE_lean)

# Bar plot - Graph predicted TEE adjusted for Lean mass #
barplot_emm_TEE_lean <- emm_TEE_lean_df %>%
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
    title = "Total energy expenditure adjusted for lean mass",
    y = "Adjusted TEE (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_TEE_lean
#Conclusion: Significant difference at BW loss, maintenance, and regain --> makes sense

###-----TEE adj. BW (minute summation method)-----####
#Build multiple linear regression model
    model_TEE_BW <- lmer(TEE_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR3)
    summary(model_TEE_BW)

# Calculate estimated marginal means #
emm_TEE_BW <- emmeans(model_TEE_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_BW_df <- as.data.frame(emm_TEE_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_BW <- contrast(emm_TEE_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_BW_df <- as.data.frame(contrasts_by_group_TEE_BW)
      
# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_BW <- contrast(emm_TEE_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
      contrasts_SABLE_TEE_BW_df <- as.data.frame(contrasts_by_SABLE_TEE_BW)

# Bar plot - Graph predicted TEE adjusted for Body weight #
  barplot_emm_TEE_BW <- emm_TEE_BW_df %>%
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
          title = "Total energy expenditure adjusted for body weight",
          y = "Adjusted TEE (kcal/day)",
          x = "Time point",
          fill = "Treatment group") +
        format.plot +
        theme(
          legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
      barplot_emm_TEE_BW
#Conclusions: at BW regain, TEE is not significantly different between weight cycled and control mice (p=0.153)
      

#----EE without adjustment for lean or BW (minute summation method)-----#
### NEAT not adj. (minute summation method) ####
  #Build multiple linear regression model for NEAT not adjusted for BW or lean #
      model_NEAT <- lmer(NEAT_teske ~ SABLE * GROUP + (1 | ID), data = Echo_TEE_RMR_NEAT3)
      summary(model_NEAT)
      
      #Calculate estimated marginal means #
      emm_NEAT <- emmeans(model_NEAT, ~ SABLE * GROUP, cov.reduce = mean)
      emm_NEAT_df <- as.data.frame(emm_NEAT)
      
      # Pairwise contrasts within each GROUP
      contrasts_by_group_NEAT <- contrast(emm_NEAT, method = "pairwise", by = "GROUP")
      contrasts_by_group_NEAT_df <- as.data.frame(contrasts_by_group_NEAT)
      
      # Pairwise contrasts within each SABLE (time point)
      contrasts_by_SABLE_NEAT <- contrast(emm_NEAT, method = "pairwise", by = "SABLE")
      contrasts_SABLE_NEAT_df <- as.data.frame(contrasts_by_SABLE_NEAT)
      
  # Bar plot - Graph predicted NEAT not adjusted for BW or lean #
      barplot_emm_NEAT <- emm_NEAT_df %>%
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
          title = "NEAT (kcal/day) unadjusted",
          y = "NEAT (kcal/day)",
          x = "Time point",
          fill = "Treatment group") +
        format.plot +
        theme(legend.position = "top",
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1))
      barplot_emm_NEAT
  #Results: at BW regain, Control vs Weight cycled --> p=0.042 --> expected!
      #Weight cycled is also significantly lower than control at BW loss and maintenance
      
###-----RMR not adj. (minute summation method) -----####
#unadjusted RMR: Build multiple linear regression model for RMR not adjusted for BW or lean #
      #Build multiple linear regression model
      model_RMR <- lmer(RMR_teske ~ SABLE * GROUP + (1 | ID), data = Echo_TEE_RMR3)
      summary(model_RMR)
      
      # Calculate estimated marginal means #
      emm_RMR <- emmeans(model_RMR, ~ SABLE * GROUP, cov.reduce = mean)
      emm_RMR_df <- as.data.frame(emm_RMR)
      
      # Pairwise contrasts within each GROUP
      contrasts_by_group_RMR <- contrast(emm_RMR, method = "pairwise", by = "GROUP")
      contrasts_by_group_RMR_df <- as.data.frame(contrasts_by_group_RMR)
      
      # Pairwise contrasts within each SABLE (time point)
      contrasts_by_SABLE_RMR <- contrast(emm_RMR, method = "pairwise", by = "SABLE")
      # Convert to a data frame
      contrasts_SABLE_RMR_df <- as.data.frame(contrasts_by_SABLE_RMR)
      
      # Bar plot - Graph predicted RMR unadjusted for lean or BW #
      barplot_emm_RMR <- emm_RMR_df %>%
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
          title = "RMR (kcal/day) unadjusted",
          y = "RMR (kcal/day)",
          x = "Time point",
          fill = "Treatment group") +
        format.plot +
        theme(legend.position = "top",
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1))
      barplot_emm_RMR
  #Conclusion: Control and weight cycled are not significantly different at any time point
      #Both control and weight cycled seem to decrease over time--almost like RMR declines with age
      
###-----TEE not adj. (minute summation method)-----####
#unadjusted TEE: Build multiple linear regression model for TEE not adjusted for BW or lean #
      #Build multiple linear regression model
      model_TEE <- lmer(TEE_teske ~ SABLE * GROUP + (1 | ID), data = Echo_TEE_RMR3)
      summary(model_TEE)
      
      # Calculate estimated marginal means #
      emm_TEE <- emmeans(model_TEE, ~ SABLE * GROUP, cov.reduce = mean)
      emm_TEE_df <- as.data.frame(emm_TEE)
      
      # Pairwise contrasts within each GROUP
      contrasts_by_group_TEE <- contrast(emm_TEE, method = "pairwise", by = "GROUP")
      contrasts_by_group_TEE_df <- as.data.frame(contrasts_by_group_TEE)
      
      # Pairwise contrasts within each SABLE (time point)
      contrasts_by_SABLE_TEE <- contrast(emm_TEE, method = "pairwise", by = "SABLE")
      contrasts_SABLE_TEE_df <- as.data.frame(contrasts_by_SABLE_TEE)
      
  # Bar plot - Graph predicted TEE unadjusted for lean or BW #
      barplot_emm_TEE <- emm_TEE_df %>%
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
          title = "Total energy expenditure (kcal/day) unadjusted",
          y = "TEE (kcal/day)",
          x = "Time point",
          fill = "Treatment group") +
        format.plot +
        theme(
          legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
      barplot_emm_TEE
  #Conclusion: at BW regain, control is significantly higher than weight cycled (p=0.0049)
      #At BW loss and BW maintenance control is also significantly higher than weight cycled
      
#---  
# Calculate hourly RMR for each ID at each SABLE time point ####
#Use minute by minute summation method, but calculate hourly and photoperiod sums rather than daily
#Approach: Filter for minutes when mouse didn't move. During these minutes TEE is entirely RMR
# (i.e.) physical activity is not contributing to TEE when mouse isn't moving
      
#By hour: average TEE (kcal/hr) each hour during one complete day for each ID at each SABLE time point
      ID_TEE3_hr <- filter_locom_energy3 %>%
        ungroup() %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        group_by(SABLE, ID, GROUP) %>%
        mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
        #filter(move==0) %>% #only keep data from minutes when the mouse moved
        ungroup() %>% 
        group_by(SABLE, ID, hr) %>%  #want kcals each hour
        #summarise(TEE_teske_hr = sum(Kcal_Hr)*(1/60), .groups="drop") 
        summarise(n_obs = n(),  # count observations per group
                TEE_teske_hr = sum(Kcal_Hr) / n_obs, .groups = "drop")
      #summarise line above adds the EE measured each minute for each minute of observation. 
      #This gives a value for EE in kcal_hr, then divide by the number of observations--> gives the 
      #avg hourly TEE for that hour (for each ID at each Sable) in kcal/hr. This is 
      #arguably better than directly summing because missing observations don't artificially reduce EE
      #Note: in the rest of the script I calculated the Kcal_min for each minute and then added those 
      #together for each hour --> that is NOT what I did in this section
      
#By photo period: avg. TEE (kcal/hr) in light/dark during one complete day for each ID at each SABLE time point
      ID_TEE3_photo <- filter_locom_energy3 %>%
        ungroup() %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        group_by(SABLE, ID, GROUP) %>%
        mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
        #filter(move==0) %>% #only keep data from minutes when the mouse moved
        ungroup() %>% 
        group_by(SABLE, ID, lights) %>%  #want kcals each hour
        summarise(n_obs = n(),  # count observations per group
                  TEE_teske_photo = sum(Kcal_Hr) / n_obs, .groups = "drop")

ggplot(mod_RMR_day, aes(x = SABLE, y = TEE_teske_photo, fill = GROUP)) +
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
     
      
      #RMR by photoperiod during one complete day for each ID at each SABLE time point
      ID_RMR3_photo <- ID_RMR3_hr %>%
        ungroup() %>%
        group_by(SABLE, ID, GROUP, lights) %>%
        summarise(RMR_teske_photo = mean(RMR_teske_hr))

      
      #### Calc. daily NEAT for each mouse and for each GROUP at each SABLE time point ####
      #use summation method (sum minute values for NEAT to get daily NEAT)
      
      #NEAT across one complete day for each mouse at each SABLE time point
      ID_NEAT3 <- filter_locom_energy3 %>%
        ungroup() %>%
        group_by(SABLE, ID, GROUP) %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
        filter(move==1) %>% #only keep data from minutes when the mouse moved
        mutate(NEAT_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
        summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop")
      
      #Average daily NEAT within each GROUP at each SABLE time point
      GROUP_NEAT3 <- filter_locom_energy3 %>%
        ungroup() %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        group_by(SABLE, ID, GROUP) %>%
        mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
        filter(move==1) %>%
        mutate(NEAT_per_min= (Kcal_Hr/60)) %>%
        summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop") %>%
        ungroup() %>%
        group_by(SABLE, GROUP) %>%
        summarise(GROUP_NEAT = mean(NEAT_teske), 
                  SD_GROUP_NEAT = sd(NEAT_teske),
                  SE_GROUP_NEAT = sd(NEAT_teske)/sqrt(n()))
      
      #Calculate daily TEE for each mouse and for each GROUP at each SABLE time point ####
      #Use summation method (sum minute values for EE to get daily EE)
      
      #TEE across one complete day for each mouse at each SABLE time point
      ID_TEE3 <- filter_locom_energy3 %>%
        ungroup() %>%
        group_by(SABLE, ID, GROUP) %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        mutate(TEE_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
        summarise(TEE_teske = sum(TEE_per_min), .groups="drop")
      
      #Average daily TEE within each GROUP at each SABLE time point
      GROUP_TEE3 <- filter_locom_energy3 %>%
        ungroup() %>%
        arrange(DateTime) %>%     # make sure rows are in time order
        group_by(SABLE, ID, GROUP) %>%
        mutate(TEE_per_min= (Kcal_Hr/60)) %>%
        summarise(TEE_teske = sum(TEE_per_min), .groups="drop") %>%
        ungroup() %>%
        group_by(SABLE, GROUP) %>%
        summarise(GROUP_TEE = mean(TEE_teske), 
                  SD_GROUP_TEE = sd(TEE_teske),
                  SE_GROUP_TEE = sd(TEE_teske)/sqrt(n()))
#---------------------------------------------#
#LOCOMOTION...should do pedmeters for locomotion actually ####
#Look at distance traveled during each time period for each group of mice. Also look at this during day and night time
#AllMeters is not zero for the first measurement of a complete day. 
Distance1 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP, LM_complete_day) %>%
  arrange(DateTime) %>%
  mutate(distance_m = All_meters-lag(All_meters),
         distance_m = if_else(distance_m <0, 0, distance_m),
         moving_min = if_else(distance_m >0, 1, 0)) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP, LM_complete_day, hr) %>%
  drop_na() %>%
  summarise(total_distance = sum(distance_m), 
            total_moving_min = sum(moving_min), .groups="drop")

ggplot(Distance1, aes(x = hr, y = total_moving_min, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(size = 2, alpha = 0.3, position = position_dodge(width = 0.8)) +
  labs(
    title= "Total time spent moving in 24hrs",
    x = "Time point",
    y = "Time (minutes)",
    fill = "Treatment group") +
  theme_minimal() +
  format.plot +
  facet_wrap(~SABLE) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))




#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")


#---------------------------------------------------------#
#---------------------------------------------------------#
#Code to make df called sable_TEE_adj_RMR (from script: NZO_Figure7b-RMR_correctedbyLean_(LM, 12-5)) ####
#### Resting metabolic rate (RMR): Identify 30min with lowest avg. TEE 
#Use the df created in this chunk (sable_RMR_data) for the code that identifies the 30min window
#This code is basically creating sable_TEE_data, but the steps that calculate avg tee are deleted

sable_RMR_data <- sable_dwn %>% 
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
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  group_by(ID, SABLE) %>% 
  # reattach GROUP and DRUG
  mutate(GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47"))

sable_RMR_data <- sable_RMR_data %>%
  mutate(SABLE = factor(SABLE, 
  levels = c("Baseline", "Peak obesity","BW loss","BW maintenance","BW regain")))

#### Code including all IDs at all SABLE time points #
#This step requires a df that has minute data for TEE rather than pre-calculated average daily TEE

# Compute sliding 30-minute averages for each mouse *and* period
lowest_windows_all <- sable_RMR_data %>%
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(avg_30min_value = slide_dbl(
    .x = value,
    .f = mean,
    .before = 29,       # previous 29 rows + current = 30 minutes
    .complete = TRUE),
    window_end_time   = DateTime,
    window_start_time = DateTime - minutes(29)) %>%
  filter(!is.na(avg_30min_value)) %>%
  ungroup()

# For each ID and SABLE, find the lowest 30-minute average
lowest_windows_summary <- lowest_windows_all %>%
  group_by(ID, SABLE) %>%
  slice_min(avg_30min_value, n = 1) %>%
  ungroup() %>%
  select(ID, SABLE, window_start_time, window_end_time, avg_30min_value)

# View summary
lowest_windows_summary

#### Change RMR units from kcal_hr to kcal_day to match tee units#
lowest_windows_summary <- lowest_windows_summary %>%
  rename(RMR_kcal_hr = avg_30min_value) %>%
  mutate(RMR_kcal_day = RMR_kcal_hr*24) %>%
  group_by(ID, SABLE)


##### Process sable_dwn into sable_TEE_data to get Avg daily TEE (tee) for each mouse at each time point#
# build the summarized dataset 
#version with creation of tee 
#Join it with echo data to create sable_tee_adj 
#Then join sable_tee_adj with lowest_windows_summary (i.e. df with RMR) to create sable_TEE_adj_RMR
#use this compiled code to do tee-BMR and for linear regression models and graphing)
sable_TEE_data <- sable_dwn %>% 
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
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>%
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47"))

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))


#### Attach echoMRI_data to sable_TEE_data --> name new df as sable_TEE_adj 
#echo info is from NZO_Figure7 - TEE_correctedbyLean (rev. LM).R on (LM accessed on 10-16-25)

#Process echoMRI info for NZO mice
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
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join Lean, Fat, and Weight info into TEE data set
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat),
    by = c("ID", "SABLE"))

##### Combine lowest_windows_summary with sable_TEE_adj #
sable_TEE_adj_RMR <- sable_TEE_adj %>%
  left_join(
    lowest_windows_summary %>% 
      select(ID, SABLE, window_start_time, window_end_time, RMR_kcal_day),
    by = c("ID", "SABLE")) %>%
  group_by(ID, SABLE) %>%
  mutate(TEE_minus_RMR = tee - RMR_kcal_day)



