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
  filter(date %in% c("2026-08-27", "2026-08-28")) %>%
  filter(!(DateTime > "2026-08-28 10:30:00"))

#Look at data starting at 5pm since mice had about 30min to normalize behavior after recording started
Test_loc_EE_60min <- Test_loc_EE %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(
    start_recording = min(DateTime),
    minutes_post_recording = as.numeric(difftime(DateTime, start_recording, units = "mins")))
  ungroup() %>%
  group_by(ID, DateTime) %>% 
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) 
  #drop_na(minutes_post_recording)

#---------------------------------------------------.#
#---------------------------------------------------.#
# Part 2: Test Mass Monitors ####
#---------------------------------------------------.#

## Test BodyMass ####
## Test FoodA ####
## Test Water ####