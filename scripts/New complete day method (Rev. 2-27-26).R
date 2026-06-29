# Began with a secion of code from New RMR method & locom, (2-26-26 issues).R
# Started: 2-26-26

#Goals: calculate TEE, RMR, and NEAT by summing one complete day for each ID at each SABLE point
#criteria for selecting a complete day was minimizing the number of missing observations 
#(i.e. minimizing the number of minutes without data)
#compared TEE and RMR calculated using the new method to TEE and RMR calculated using the old method
#(old method was used for my prelim exam)

#Method & reasoning used to select complete days is in script: 
    #New RMR method & locom, (2-26-26 issues).R)

#Additional goals: calculate meters moved and time spent moving during a complete day

#Left off ### 
#To do --> multiple linear regression and bar plot for TEE, RMR, NEAT, locomotion distance, locomotion minutes

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

#Get complete day 1 and 2 by explicitly definining the zt_time and SABLE_DAY
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

#### Calculate daily RMR for each mouse and for each GROUP at each SABLE time point ####
  #Use summation method
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

#### Calculate daily NEAT for each mouse and for each GROUP at each SABLE time point ####
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

#### Calculate daily TEE for each mouse and for each GROUP at each SABLE time point ####
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
#which used averages (Scripts: "NZO_figure7 - TEE corrected by BW (rev. LM)" and NZO_Figure7b-RMR_correctedbyLean_(LM, 12-5))

#Original method for TEE: calculated avg TEE for each ID across all min of ~48hrs. 
  #df called sable_TEE_adj_BW df found in script "NZO_figure7 - TEE corrected by BW (rev. LM).r"
  #Extrapolate this out to 24hrs to get "tee" --> tee is the extrapolation for daily TEE for each mouse ID
#Average daily TEE within each GROUP at each SABLE time point (original method)
old_GROUP_TEE <- sable_TEE_adj_BW %>%
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
#Average daily RMR within each GROUP at each SABLE time point (original method) 
#RMR_kcal_day is per ID and is a variable located in df called sable_TEE_adj_RMR
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

#First, need to process/format echoMRI data
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
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28"))) %>% #repeated
  mutate(STATUS = factor(STATUS,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  rename(SABLE = STATUS) 

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

#Confirm that data frame has correct number of mice
n_distinct(Echo_TEE_RMR3$ID) #good we have 16 animals
n_distinct(Echo_TEE_RMR_NEAT3$ID) #good we have 16 animals


##-----MLR: EE (summation) adjusted for BW or Lean mass-----####
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
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

###-----RMR (summation method) and adiposity index (fat/lean)-----##
#Build multiple linear regression model
model_RMR_AI <- lmer(RMR_teske ~ SABLE * GROUP + adiposity_index + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_AI)

#### Calculate estimated marginal means ####
emm_RMR_AI <- emmeans(model_RMR_AI, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_AI_df <- as.data.frame(emm_RMR_AI)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_AI <- contrast(emm_RMR_AI, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_AI_df <- as.data.frame(contrasts_by_group_RMR_AI)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_AI <- contrast(emm_RMR_AI, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_AI_df <- as.data.frame(contrasts_by_SABLE_RMR_AI)

#### Bar plot - Graph predicted RMR adjusted for Lean mass ####
barplot_emm_RMR_AI <- emm_RMR_AI_df %>%
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
    title = "Resting metabolic rate (RMR) adjusted for adiposity index (fat mass/lean mass)",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_AI

#For reference: adiposity index graph
#start delete
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
    title="Adiposity index (fat mass/lean mass) over weight cycle",
    x = "Time point",
    y = "Adiposity index (fat/lean)",
    fill = "Treatment group")
#end delete

###-----NEAT (summation method) and Lean mass-----##
#Build multiple linear regression model
model_NEAT_Lean <- lmer(NEAT_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR_NEAT3)
summary(model_NEAT_Lean)

#### Calculate estimated marginal means ####
emm_NEAT_Lean<- emmeans(model_NEAT_Lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_Lean_df <- as.data.frame(emm_NEAT_Lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_Lean <- contrast(emm_NEAT_Lean, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_Lean_df <- as.data.frame(contrasts_by_group_NEAT_Lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_Lean <- contrast(emm_NEAT_Lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_NEAT_Lean_df <- as.data.frame(contrasts_by_SABLE_NEAT_Lean)

#### Bar plot - Graph predicted NEAT adjusted for Lean mass ####
barplot_emm_NEAT_Lean <- emm_NEAT_Lean_df %>%
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
    title = "Non-exercise activity thermogenesis (NEAT) adjusted for lean mass",
    y = "Adjusted NEAT (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_NEAT_Lean

###-----NEAT (summation method) and BW-----##
#Build multiple linear regression model
model_NEAT_BW <- lmer(NEAT_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR_NEAT3)
summary(model_NEAT_BW)

#### Calculate estimated marginal means ####
emm_NEAT_BW<- emmeans(model_NEAT_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_BW_df <- as.data.frame(emm_NEAT_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_BW <- contrast(emm_NEAT_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_BW_df <- as.data.frame(contrasts_by_group_NEAT_BW)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_BW <- contrast(emm_NEAT_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_NEAT_BW_df <- as.data.frame(contrasts_by_SABLE_NEAT_BW)

#### Bar plot - Graph predicted NEAT adjusted for BW ####
barplot_emm_NEAT_BW <- emm_NEAT_BW_df %>%
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
    title = "Non-exercise activity thermogenesis (NEAT) adjusted for body weight",
    y = "Adjusted NEAT (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_NEAT_BW

###-----RMR (summation method) and BW-----##
#Build multiple linear regression model
model_RMR_BW <- lmer(RMR_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_BW)

#### Calculate estimated marginal means ####
emm_RMR_BW <- emmeans(model_RMR_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_BW_df <- as.data.frame(emm_RMR_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_BW <- contrast(emm_RMR_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_BW_df <- as.data.frame(contrasts_by_group_RMR_BW)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_BW <- contrast(emm_RMR_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_BW_df <- as.data.frame(contrasts_by_SABLE_RMR_BW)

#### Bar plot - Graph predicted RMR adjusted for Lean mass ####
barplot_emm_RMR_BW <- emm_RMR_BW_df %>%
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
    title = "Resting metabolic rate (RMR) adjusted for body weight",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_BW

###-----RMR (summation method) and lean mass-----##
#Build multiple linear regression model
model_RMR_lean <- lmer(RMR_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR3)
summary(model_RMR_lean)

#### Calculate estimated marginal means ####
emm_RMR_lean <- emmeans(model_RMR_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_lean_df <- as.data.frame(emm_RMR_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_lean <- contrast(emm_RMR_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_lean_df <- as.data.frame(contrasts_by_group_RMR_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_lean <- contrast(emm_RMR_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_lean_df <- as.data.frame(contrasts_by_SABLE_RMR_lean)

#### Bar plot - Graph predicted RMR adjusted for Lean mass ####
barplot_emm_RMR_lean <- emm_RMR_lean_df %>%
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
    title = "Resting metabolic rate (RMR) adjusted for lean mass",
    y = "Adjusted RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_lean

###-----TEE (summation method) and Lean mass-----##
#Build multiple linear regression model
model_TEE_lean <- lmer(TEE_teske ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_RMR3)
summary(model_TEE_lean)

#### Calculate estimated marginal means ####
emm_TEE_lean <- emmeans(model_TEE_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_lean_df <- as.data.frame(emm_TEE_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_lean <- contrast(emm_TEE_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_lean_df <- as.data.frame(contrasts_by_group_TEE_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_lean <- contrast(emm_TEE_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_TEE_lean_df <- as.data.frame(contrasts_by_SABLE_TEE_lean)


#### Bar plot - Graph predicted TEE adjusted for Lean mass ####
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

###-----TEE (summation method) and BW-----##
#Build multiple linear regression model
    model_TEE_BW <- lmer(TEE_teske ~ SABLE * GROUP + Weight + (1 | ID), data = Echo_TEE_RMR3)
    summary(model_TEE_BW)

# Calculate estimated marginal means ####
emm_TEE_BW <- emmeans(model_TEE_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_BW_df <- as.data.frame(emm_TEE_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_BW <- contrast(emm_TEE_BW, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_BW_df <- as.data.frame(contrasts_by_group_TEE_BW)
      
# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_BW <- contrast(emm_TEE_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
      contrasts_SABLE_TEE_BW_df <- as.data.frame(contrasts_by_SABLE_TEE_BW)

# Bar plot - Graph predicted TEE adjusted for Body weight
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

#---------------------------------------------
------------------####LOCOMOTION...should do pedmeters for locomotion actually ####--------------------
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
  scaleFill+
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
#---------------------------------------------

#FI over LM_complete_day =1 or 2
#This chunk of code applies my new method for defining a complete day to FI during Sable recording


#Combine code for EE and locomotion in the top part of this script with FI in sable from NZO_FI_during_Sable (2-23-26)

#------- Get only the FoodA parameter --> value column has FoodA data
filter_FI1 <-sable_dwn %>%
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
  filter(grepl("FoodA_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) 

#Get complete day 1 and 2 by explicitly definining the zt_time and SABLE_DAY
filter_FI2 <- filter_FI1 %>% 
  mutate(LM_complete_day = case_when( 
    #Baseline
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
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_FI3 <- filter_FI2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2)) #remove observations not from complete day 1 or 2

#Check number of mice in filter_FI3
filter_FI3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #we have 16 NZO in all SABLE periods

#---Make df called filter_loc_EE_FI3 with: sable FI, locomotion, and EE 

#In df filter_food3 use mutate to make a column called Food_g
filter_food3 <- filter_FI3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Food_g = value) %>%
  rename(parameter_FoodA = parameter) %>%
  rename(fix_value_FoodA = fix_value)

#Join filter_food3 and filter_locom_energy3 into a df called filter_loc_EE_FI3
  #Add Food_g to filter_locom_energy3 (by ID, DateTime, sable_idx)
filter_loc_EE_FI3 <- filter_locom_energy3 %>%
  left_join(
    filter_food3 %>% 
      select(Food_g, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Do FI calculations. Then after these are done attach the manual FI calculations ####
#Sum FI each minute to get total food FI on complete day 1 and complete day 2
#(Note: the source df for this section already removed incomplete days)
filter_sable_FI3 <- filter_loc_EE_FI3 %>%
  group_by(ID, LM_complete_day, DRUG, GROUP, SABLE) %>%
  arrange(ID, DateTime) %>%       # make sure data is ordered
  mutate(
    intake = lag(Food_g) - Food_g,       # change in food from one minute to the next
    intake = if_else(intake < 0, 0, intake)) %>%  #If mouse doesn't eat between min. x and min. x+1, intake=0
  mutate (intake = if_else(intake >8, 0, intake)) %>% #removes values that are illogically high for 1 minute
  drop_na() %>% 
  select(-KCAL_PER_GR) %>% # this is wrong, so remove and replace with KCAL_G
  summarise(
    total_eaten_gr = sum(intake),      # total FI per day in grams
    .groups = "drop") %>%
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
  mutate(total_eaten_kcal = total_eaten_gr*KCAL_G) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>% 
  summarise (avg_corrected_intake_kcal = mean(total_eaten_kcal)) %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice

#### MANUAL FI for weight cycled MICE at BW loss and BW maintenance ####
#Create df "FI_manual_cycled" -> Process FI.csv file (i.e. manual measurements of FI)
FI_manual_cycled <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(sable_idx = case_when(
      #BW loss 
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      ID %in% c(3710) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      ID %in% c(3708, 3714) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      #BW maintenance (SABLE_DAY_18 for all mice --> complete day 2) 
      ID %in% c(3708, 3710) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      ID %in% c(3714, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      ID %in% c(3722, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
      sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  ungroup() 

Manual_FI_cycled <- FI_manual_cycled %>%
  mutate(ID = factor(ID)) %>%
  filter(GROUP == "Weight cycled") %>%
  filter(SABLE %in% c("BW loss", "BW maintenance")) %>% 
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  mutate(INTAKE_kcal = (INTAKE_GR*3.82)) %>%
  summarise(avg_corrected_intake_kcal = mean(INTAKE_kcal))

#Combine FI measured by cages with manual measurements from food restriction
FI_in_Sable3 <- bind_rows(Manual_FI_cycled, filter_sable_FI3)

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
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

#Bar graph of FI calculated using avg FI during all 4 days of Sable 
#df=FI_in_Sable is a compilation of manual FI for weight cycled mice 
#during BW loss/maintenance (Manual_FI) and FI directly recorded by Sable
ggplot(FI_in_Sable3, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
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
  geom_point(size = 2, alpha = 0.3, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = custom_colors) +
  labs(
    title= "Avg daily (LM_complete_day) food intake",
    x = "Time point",
    y = "Average daily food intake (kcal/day)",
    fill = "Treatment group") +
  format.plot +
  theme_bw(base_size = 10) +
  theme(
    legend.position ="top",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))

#Energy balance ####
#join df=FI_in_Sable3 and df=ID_TEE3
Energy_balance3 <- FI_in_Sable3 %>%
  left_join(
    ID_TEE3 %>% 
      select(TEE_teske, SABLE, ID, GROUP),
    by = c("SABLE", "ID", "GROUP")) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  mutate (EE_balance = avg_corrected_intake_kcal-TEE_teske) %>%
  filter(!(EE_balance >20))

ggplot(Energy_balance3, aes(x = SABLE, y = EE_balance, fill = GROUP)) +
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
  geom_point(size = 2, alpha = 0.3, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = custom_colors) +
  labs(
    title= "Energy balance = FI - TEE (LM_complete_day)",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group") +
  format.plot +
  theme_bw(base_size = 10) +
  theme(
    legend.position ="top",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))
