#Quality control to ensure that gas analyzers, mass monitors, and beam breaks 
#are working after the installation of the GA that was repaired second, after
#readdressing mass monitors so that non-functional MMs are address #2 (water),
#and after calibrating the functional MMs

#
#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

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
library(hms)

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#Read in Sable data ####
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
##Filter Sable data to include only mice that are cannulated
sable_dwn_19 <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()

#Test run Sable data (~4:30pm 8/27 to ~10:45am 8/28)

#---------------------------------------------------.#
#---------------------------------------------------.#
# Part 1: Test EE and locomotion ####
#---------------------------------------------------.#

##Test kcal_hr (aka EE, kcal/hr) ####
Test_EE_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

##Test AllMeters (aka all movement) ####
Test_loc_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

#Test PedMeters (aka movement in XY plane) ####
Test_ped_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("PedMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

#---
#Join EE and locomotion data 
Test_loc_EE <- Test_loc_1 %>%
  left_join(
    Test_EE_1 %>% 
      select(Kcal_Hr, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  left_join(
    Test_ped_1 %>% 
      select(Ped_meters, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(ID=as.factor(ID)) 

#Filter for data collected during the test period (8/27 and 8/28)
#Also remove data that will be recorded in the afternoon on 8/28/2026
Test_loc_EE_2 <- Test_loc_EE %>%
  ungroup()%>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  filter(date %in% c("2026-08-27", "2026-08-28")) %>%
  #filter(!(DateTime > "2026-08-28 10:30:00")) %>%
  filter(DateTime > "2026-08-28 10:30:00")
  

#Look at data starting at 5pm since mice had about 30min to normalize behavior after recording started
Test_loc_EE_60min <- Test_loc_EE_2 %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(
    start_recording = min(DateTime),
    minutes_post_recording = as.numeric(difftime(DateTime, start_recording, units = "mins"))) %>%
  ungroup() %>%
  group_by(ID, DateTime) %>% 
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) 

Test_sum_EE_hr_bins <- Test_loc_EE_60min %>%
  ungroup()%>%
  #group_by(ID, recording_bin) %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  summarise(sum_EE_60min = sum(TEE_per_min)) 

#Calculate total EE for the afternoon run on 8/28/26 (after gas calibration)
post_calib_EE_sum <- Test_loc_EE_2 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  summarise(sum_EE_60min = sum(TEE_per_min)) 

#---------------------------------------------------.#
#---------------------------------------------------.#
# Part 2: Test Mass Monitors ####
#---------------------------------------------------.#

## Test BodyMass ####
## Test FoodA ####
## Test Water ####

#---------------------------------------------------.#
#---------------------------------------------------.#
# Part 3: Compare set up 1 and set up 2 ####
#compare same mice in different cages: 
#Complete days created: incomplete (Fri), Fri->Sat, Sat->Sun, Sun->Mon, incomplete (Mon->Tues)
#---------------------------------------------------.#
# Set up 1: Use  6pm on Friday to 6pm on Saturday 
# Set up 2: Use 6pm on Sunday to 6pm on Monday 

## TEE from sable_dwn
TEE_quality <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date %in% c("2026-08-28", "2026-08-29", "2026-08-30", "2026-08-31")) %>% 
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
         lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
         SET_UP= case_when(
           date %in% c("2026-08-28", "2026-08-29") ~ "one",
           date %in% c("2026-08-30", "2026-08-31") ~ "two")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, SET_UP, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, complete_days, SET_UP) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),min(DateTime), units="hours")), zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#All movement (AllMeters) from sable_dwn
#(Once I run TEE from above and confirm that it works then copy and paste for locomotion)
ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

#Ambulation (PedMeters)
loc_quality <-sable_dwn %>%

  
#Join AllMeters and PedMeters
Ped_quality <-sable_dwn %>%

  
#Join EE and locomotion data 
Quality_loc_EE <- loc_quality %>%
  left_join(
    TEE_quality %>% 
      select(Kcal_Hr, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  left_join(
    Ped_quality %>% 
      select(Ped_meters, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(ID=as.factor(ID)) 



#Calculate EE (daily sum)...could also look at mean or median hourly EE





