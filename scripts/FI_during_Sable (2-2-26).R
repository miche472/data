#NZO food intake and energy balance while recording in sable cages 

  #This script started as a cleaned up version of NZO_FI (LM Rev. 1-30-26)
#Script started: 1-30-26
#Date of last revision: 1-30-26
  #In this script I look at the concurrent FI and EE rather than FI during time periods between Sable measurements
# Bottom of this script:
  #Used mouse 3706 to establish a script for extracting FI data from sable recording ("Establishing method using 3706")

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
#### Functions and source data #####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}
#------- sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 


#### Process FI recorded by Sable (all IDs and time points) ####
    #code written on 1-28-26
#This chunk of code processes sable_dwn to sable_FI_data 
#sable_FI_data only has FoodA under the value column and only includes complete days 1 and 2
sable_FI_data <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light"),
    SABLE = case_when(
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                               "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                               "SABLE_DAY_7") ~ "Baseline",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("FoodA_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  group_by(ID, DRUG,GROUP,DIET_FORMULA,SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,DIET_FORMULA,SABLE) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter( is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain")))

#Sum the mass of food eaten each minute to get total food eaten on complete days 1 and 2
  #sable_FI_minutes excludes mice that were recorded in Sable cages that were problematic
sable_FI_minutes <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, GROUP, SABLE) %>%
  #filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725=cage5 issues and 3724=cage 6
  #filter(complete_days == "2") %>% #Use this line to calculate just one of the two complete days
  arrange(ID, DateTime, hr, complete_days) %>%       # make sure data is ordered
  group_by(ID, complete_days, DRUG, GROUP, SABLE) %>% 
  mutate(
    intake = lag(value) - value,       # change in food from one minute to the next
    intake = if_else(intake < 0, 0, intake)) %>%  #If mouse doesn't eat between min. x and min. x+1, intake=0
  drop_na() %>% 
  summarise(
    total_eaten_gr = sum(intake),      # total FI per day in grams
    .groups = "drop") %>%
  mutate(total_eaten_kcal = (total_eaten_gr)*3.82) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise (avg_corrected_intake_gr = mean(total_eaten_gr), .groups = "drop") %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice

#Average daily FI (g) recorded by Sable for each GROUP
  #Avg_Sable_FI doesn't include accurate values for Weight cycled mice at BW loss/maintenance
Avg_Sable_FI <- sable_FI_minutes %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_avg_corrected_intake_gr= mean(avg_corrected_intake_gr), .groups="drop")

#### Use manual FI measurements from Sable recording at BW loss for Weight cycled mice ####
  #For BW loss: Sable days 12-15 (13 and 14 are the "middle days")
  #For BW maintenance: Sable days 16-19 (17 and 18 are the "middle days")

#Create df "FI_manual_SABLE" -> Process FI.csv file (manually measured FI)
    #Extract accurate FI for weight cycled mice during BW loss/maintenance
FI_manual_SABLE <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>%
  #filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
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
  mutate(
    SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
      sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) %>% 
  ungroup() 

#Calculate manual FI during BW loss -> Two possible methods to calclate avg FI (daily)
  #during Sable recording for BW loss/maintenance in Weight cycled mice

#Method 1:Average FI (g) for the two middle days of recording for each weight cycled mouse
  #at BW loss and BW maintenance time points
  #This is similar to the concept of including only complete days 1 and 2 for FI measured by Sable
Manual_FI_complete_days <- FI_manual_SABLE %>%
  mutate(ID = factor(ID)) %>%
  #filter(GROUP == "Weight cycled") %>%
  #filter(sable_idx %in% c("SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_17", "SABLE_DAY_18")) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise(avg_corrected_intake_gr = mean(corrected_intake_gr))

#Method 2: Average FI (g) across all 4 days of Sable recording for each weight cycled mouse
  #at BW loss and BW maintenance time points
Manual_FI <- FI_manual_SABLE %>%
  mutate(ID = factor(ID)) %>%
  #filter(GROUP == "Weight cycled") %>%
  #filter(SABLE %in% c("BW loss", "BW maintenance")) %>% #includes all 4 Sable days for BW loss and BW maintenance
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise(avg_corrected_intake_gr = mean(corrected_intake_gr))

#Is FI significantly different when using the 2 days vs 4 day approach? 
#Approximate the answer to this by graphing and comparing to when the complete days (i.e. 2 days of data is used)
#Make a data frame with 4 day data for weight cycled at BW loss/maintenance and sable_FI_minutes
#FI_in_Sable_gr_kcal
FI_in_Sable_4 <- bind_rows(sable_FI_minutes, Manual_FI)
#Add column with avg daily KCAL
FI_in_Sable_4_gr_kcal <- FI_in_Sable_4 %>%
  group_by(ID, SABLE, GROUP) %>%
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82),
         avg_corrected_intake_kcal = (avg_corrected_intake_gr)*(KCAL_G))

#Graph FI_in_Sable_4_gr_kcal (can do kcal or grams by changing y variable)
ggplot(FI_in_Sable_4_gr_kcal, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
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
  labs(
    title= "Four days: Average daily food intake during Sable recording",
    x = "Time point",
    y = "Food intake (kcal/day)",
    fill = "Treatment group") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))
#Using 4 days or 2 days of FI data for weight cycled mice at BW loss/maintenance looks very similar
#Is there a mathematical outlier using either approach?

#### End of 1-29-26 coding ####

# ------------------------------------------- #

##### Compiled df with FI for control and weight cycled mice during all Sable time points ####
#Wrote this section on 1-30-26
    #df1= Manual_FI_complete_days (BW loss and maintenance for weight cycled mice. Avg of 2 days)
    #df2= sable_FI_minutes (All time points for Control and baseline/peak obesity/BW regain for weight cycled mice)
#Excludes sable data recorded by cages that were problematic (excluded from sable_FI_minutes)

#FI_in_Sable_gr_kcal
  FI_in_Sable <- bind_rows(sable_FI_minutes, Manual_FI_complete_days)

  #Add column with avg daily KCAL
  FI_in_Sable_gr_kcal <- FI_in_Sable %>%
    group_by(ID, SABLE, GROUP) %>%
    mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82),
          avg_corrected_intake_kcal = (avg_corrected_intake_gr)*(KCAL_G))

#Graph FI_in_Sable_gr_kcal (can do kcal or grams by changing y variable)
ggplot(FI_in_Sable_gr_kcal, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(size = 2, alpha = 0.3, position = position_dodge(width = 0.8)) +
  labs(
    title= "Two days: Average daily food intake during Sable recording",
    x = "Time point",
    y = "Food intake (kcal/day)",
    fill = "Treatment group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))

#2-2-26 Check for outliers in FI measurements df= FI_in_Sable_gr_kcal
#Look at outliers in Excel and using statistical approaches 
#Note the Excel sheet "considering possible outliers for FI" is more complete than this section of code

# 1. Check for outliers in FI during Sable recording
#Graph boxplots
library("ggpubr")
ggboxplot(FI_in_Sable_gr_kcal, x = "SABLE", y = "avg_corrected_intake_kcal", color = "GROUP",
          palette = c("#00AFBB", "#E7B800"))

ggboxplot(FI_in_Sable_4_gr_kcal, x = "SABLE", y = "avg_corrected_intake_kcal", color = "GROUP",
          palette = c("#00AFBB", "#E7B800"))

#Statistical tests for significant differences

#Method 1: Fixed effects ANOVA
ANOVA_sable_FI_kcal <- aov(
  avg_corrected_intake_kcal ~ GROUP * SABLE,
  data = FI_in_Sable_gr_kcal)

    # Extract studentized residuals
    FI_in_Sable_gr_kcal$stud_resid <- rstudent(ANOVA_sable_FI_kcal)
    
    # Flag potential outliers
    FI_in_Sable_gr_kcal$outlier_model <- abs(FI_in_Sable_gr_kcal$stud_resid) > 3
    FI_in_Sable_gr_kcal %>% filter(outlier_model)
    
    #View flagged results
    FI_in_Sable_gr_kcal %>% filter(outlier_model)
    #revealed two outliers: 3713 at Baseline (kcal=72.2) and 3720 at BW regain (kcal=53.6)

#Method 2: Mixed effects model
    library(lme4)
    
    lmm_sable_FI_kcal <- lmer(
      avg_corrected_intake_kcal ~ GROUP * SABLE + (1 | ID),
      data = FI_in_Sable_gr_kcal)
    # UPDATE: need to addadd statistical test for linear mixed model ####

#####Graph the manual FI measurements during Sable recording ####

#Add column with FI in kcal
    #Manual_FI_gr_kcal
    #Add column with avg daily KCAL
    Manual_FI_gr_kcal <- Manual_FI %>%
      ungroup() %>%
      group_by(ID, SABLE, GROUP) %>%
      mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82),
             avg_corrected_intake_kcal = (avg_corrected_intake_gr)*(KCAL_G))

#Manual_FI
    #Box plot
ggboxplot(Manual_FI_gr_kcal, x = "SABLE", y = "avg_corrected_intake_kcal", color = "GROUP",
          palette = c("#FAAC41","#5392DB"))

#Bar bar plot with data points (using manual food intake for all time points and mice)
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
ggplot(Manual_FI_gr_kcal, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
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
    title= "Manually measured average daily food intake during Sable recording",
    x = "Time point",
    y = "Food intake (kcal/day)",
    fill = "Treatment group") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))


#### Next step --> energy balance ####
# 2.  Energy balance
      # Extract code that calculated raw EE at each time point
      # Compare the raw FI and EE (i.e. unadjusted for BW/lean mass) for the groups at the different time points

#Use the section of this code where I was more accurate with the Sable dates for each mouse
      #to revise the code in NZO_Figure7b-RMR_correctedbyLean to also reflect these dates
      #then filter Manual_FI_gr_kcal and FI_in_Sable_gr_kcal to remove mice with problematic Sable recordings
      #Use left join to combine TEE/RMR data frame and each of the two FI data frames 
      #Calculate energy balance using each of the FI data frames and determine which is more logical


# ------------------------------------------- #
#This section of code was written on 1/27/26, before any of the code above
#### Establish method using 3706 ####
#I retained the notes which describe my train of logic while writing the code
#This section of code was written on 1/27/26

#Before running this you need to load in zt_time function and sable_dwn, located at the top of the script
#Start of 1/27attempt using section from locomotion_LM_12-2.R to do initial processing of sable_dwn to select just FoodA
sable_FI_data <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light"),
    SABLE = case_when(
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                               "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                               "SABLE_DAY_7") ~ "Baseline",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("FoodA_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  group_by(ID, DRUG,GROUP,DIET_FORMULA,SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,DIET_FORMULA,SABLE) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter( is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain")))
#This created a df from sable_dwn with only FoodA under the "value" column

#####This method worked (manually measured intake =10.42g and Sable measured 10.89g)
#Filter for just 3706 on the first complete day of SABLE=Peak obesity
sable_FI_3706 <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  filter(complete_days == "1") %>%
  arrange(hr) %>%
  ungroup() %>%
  group_by(hr) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706 <- sable_FI_3706 %>%
  group_by(ID, GROUP, DRUG, SABLE, hr) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706 <- hr_FI_3706 %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
            .groups = "drop")

#Filter for just 3706 on the second complete day of SABLE=Peak obesity
sable_FI_3706_2 <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  filter(complete_days == "2") %>%
  arrange(hr) %>% #should I arrange by DateTime instead since hrs 20-23 actually occurred before hrs 0-19?
  ungroup() %>%
  group_by(hr) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_2 <- sable_FI_3706_2 %>%
  group_by(ID, GROUP, DRUG, SABLE, hr) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706_2 <- hr_FI_3706_2 %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
            .groups = "drop")

#To get the sum of FI during day 1 and day 2 using one chunk of code I grouped by complete_days 
#this approach worked (sum of days 1 and 2 = 10.89)
#Filter for 3706 and sable day 9
sable_FI_3706_bothdays <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == "1") %>% # commenting this out leaves both complete days in
  arrange(hr, date) %>% 
  ungroup() %>%
  group_by(hr, date) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_bothdays <- sable_FI_3706_bothdays %>%
  group_by(ID, GROUP, DRUG, SABLE, hr, complete_days) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706_bothdays <- hr_FI_3706_bothdays %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
            .groups = "drop")

#Here I used the same approach as immediately above, but I arranged by DateTime. I think this will be better
#for when I want to expand this code to be used for multiple time points

sable_FI_3706_cum <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == 1) %>% # commenting this out leaves both complete days in
  arrange(DateTime) %>%
  ungroup() %>%
  group_by(hr, date) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_cum <- sable_FI_3706_cum %>%
  group_by(ID, GROUP, DRUG, SABLE, hr, complete_days) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),.groups = "drop") 

sum_3706_cum <- hr_FI_3706_cum %>%
  summarise(day_cum_FI = sum(FI_cum_hr),.groups = "drop")

# Alternative approach: Sum minute by minute FI to obtain FI for complete day 1 and 2
  #this code works
sable_FI_3706_minutes <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == "2") %>% #Use this line to calculate just one of the two complete days
  arrange(ID, DateTime, hr, complete_days) %>%       # make sure data is ordered
  group_by(ID, complete_days, SABLE, GROUP, DRUG) %>% 
  mutate(
    intake = lag(value) - value,       # change in meters
    intake = if_else(intake < 0, 0, intake)) %>%  
  #the line above says that if the calculated intake is less than zero then put 0 in the cell. 
  #Otherwise, put the intake (g) which was calculated into the cell.
  drop_na() %>% 
  summarise(
    total_eaten_gr = sum(intake),      # total FI per day in grams
    .groups = "drop") 
# End of establish method using 3706 ####
