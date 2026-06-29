#NZO food intake and energy balance while recording in sable cages 

  #This script started as a cleaned up version of NZO_FI (LM Rev. 1-30-26)
#Script started: 1-30-26
#Date of last revision: 2-13-26
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
library(lme4)
#### Functions and source data #####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}
#------- sable data 
#Set working directory: session->Set working directory->Choose directory->Documents->GitHub->data->data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 


#### SABLE RECORDED FI recorded for all IDs at all time points ####
    #processes sable_dwn to sable_FI_data 
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
  filter(is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain")))

#Sum FI each minute to get total food FI on complete day 1 and complete day 2
  #(Note: the source df for this section already removed incomplete days)
sable_FI_minutes <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, GROUP, SABLE) %>%
  arrange(ID, DateTime, hr, complete_days) %>%       # make sure data is ordered
  group_by(ID, complete_days, DRUG, GROUP, SABLE) %>% 
  mutate(
    intake = lag(value) - value,       # change in food from one minute to the next
    intake = if_else(intake < 0, 0, intake)) %>%  #If mouse doesn't eat between min. x and min. x+1, intake=0
  #filter(intake > 1) #delete this line after using it to see if there are any times when there is a huge jump between minute to minute measurements
  mutate (intake = if_else(intake >8, 0, intake)) %>% #removes values that are illogically high for 1 minute
drop_na() %>% 
  select(-KCAL_PER_GR) %>%
  summarise(
    total_eaten_gr = sum(intake),      # total FI per day in grams
    .groups = "drop") %>%
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
  mutate(total_eaten_kcal = total_eaten_gr*KCAL_G) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>% #UPDATE: need to get average across the two complete_days...how to do this?
  summarise (avg_corrected_intake_kcal = mean(total_eaten_kcal)) %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice

#UPDATE: move this to the section after you attach the FI data for weight cycled mice at BW loss/maintenance ####
#Average daily FI (kcal) recorded by Sable for each GROUP 
    #(Note: currently excludes Weight cycled BW loss/maintenance)
Avg_Sable_FI <- sable_FI_minutes %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  #filter(!(SABLE == "BW regain" & ID == "3720")) %>% #remove visually apparent outlier, 3720 at BW regain (kcal=53.6)
  #filter(!(ID == "3713" & SABLE == "Baseline")) %>% #remove visually apparent outlier, 3713 at Baseline (kcal=72.2)
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_avg_corrected_intake_gr= mean(avg_corrected_intake_kcal), .groups="drop")

#### MANUAL FI for ALL MICE at all time points ####
    #(Note: this df will be further processed to isolate exclusively th FI during BW loss/maintenance 
    #Sable recording for Weight cycled mice

#Create df "FI_manual_SABLE" -> Process FI.csv file (i.e. manual measurements of FI)
    #(Note: in a later step I can extract FI for just weight cycled mice during BW loss/maintenance)
FI_manual_SABLE <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
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

#UPDATE: rename this df to highlight that it is just weight cycled mice and manual FI ####
#Filter FI_manual_SABLE to extract only values for weight cycled mice during BW loss/maintenance ####
  # Method: Average the daily FI (g) across all manual FI measurements taken during Sable recording for 
  #each weight cycled mouse at both the BW loss and then the BW maintenance time points. 
    #Note: Must use all 4 days since there were only 2 manual FI measurements during Sable recording (i.e. FI not measured on consecutive days)
  Manual_FI <- FI_manual_SABLE %>%
  mutate(ID = factor(ID)) %>%
  filter(GROUP == "Weight cycled") %>%
  filter(SABLE %in% c("BW loss", "BW maintenance")) %>% #includes all 4 Sable days for BW loss and BW maintenance
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise(avg_corrected_intake_kcal = mean(corrected_intake_kcal))

#COMPILED Sable recorded and manual FI for control and weight cycled mice during all time points ####

# Method: Combine Avg. Manual FI during Sable at BW loss and BW maintenance for weight cycled mice 
        #and average sable recorded FI on complete days 1 and 2 for other SABLE/GROUP combinations
FI_in_Sable <- bind_rows(Manual_FI, sable_FI_minutes)

# For bar graph -> Define custom colors
  custom_colors <- c("Control" = "#FAAC41", "Weight cycled" = "#3498DB")
    
#Bar graph of FI calculated using avg FI during all 4 days of Sable 
  #df=FI_in_Sable is a compilation of manual FI for weight cycled mice 
  #during BW loss/maintenance (Manual_FI) and FI directly recorded by Sable
ggplot(FI_in_Sable, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
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
    title= "Average daily food intake during Sable recording (compilation)",
    x = "Time point",
    y = "Average daily food intake (kcal/day)",
    fill = "Treatment group") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))

#2-2-26 Check for outliers in FI measurements df= FI_in_Sable ####
#Look at outliers in Excel and using statistical approaches 
#Note: Excel sheet "considering possible outliers for FI" is more complete than this section of code

# Check for outliers in FI during Sable recording
    #Graph boxplots
    ggboxplot(FI_in_Sable, x = "SABLE", y = "avg_corrected_intake_kcal", color = "GROUP",
              palette = c("#FAAC41", "#3498DB"))

    #UPDATE: Statistical tests for significant differences ####
            # Fixed effects ANOVA
            # Mixed effects model
    
# Graph the manual FI measurements during Sable recording for all mice ####
    #Recall: df=FI_manual_SABLE has MANUAL FI during recording for all mice at all time points
    
  #Method: Manual Average daily FI (g) across duration of Sable for ALL MICE
    Manual_FI_all_mice <- FI_manual_SABLE %>%
      mutate(ID = factor(ID)) %>%
      ungroup() %>%
      group_by(ID, SABLE, GROUP) %>%
      summarise(avg_corrected_intake_kcal = mean(corrected_intake_kcal)) %>%
          #Optional: remove mice with sable recording issues (to prepare for energy balance calculation)
          group_by(ID) %>%
          filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
          ungroup()

#Graph Manual_FI_all_mice
    #Box plot
    ggboxplot(Manual_FI_all_mice, x = "SABLE", y = "avg_corrected_intake_kcal", color = "GROUP",
          palette = c("#FAAC41","#5392DB"))

#Bar bar plot 
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
ggplot(Manual_FI_all_mice, aes(x = SABLE, y = avg_corrected_intake_kcal, fill = GROUP)) +
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
    y = "Avg daily food intake (kcal/day)",
    fill = "Treatment group") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5))


#### Next step --> Energy balance and Food efficiency ####
# 2-5-26 #### 1st: energy expenditure 
# 1.  Energy balance (left join data frames below to facilitate)
      #1a.  EE (sable_TEE_adj_RMR) - FI, sable recorded (FI_in_Sable)
      #1b.  EE (sable_TEE_adj_RMR) - FI, manual for all mice (Manual_FI_all_mice) 
  #Note: remove mice that had problems with Sable recording
    #Conclusion: FI_in_Sable yields more logical results than manual measurements

# Original code to generate tee and RMR data frame with echoMRI (taken from Energy balance (LM, 12-8-25).R) ####
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
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"))

sable_RMR_data <- sable_RMR_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))

#### Code including all IDs at all SABLE time points ####
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
  dplyr::select(ID, SABLE, window_start_time, window_end_time, avg_30min_value)

# View summary
lowest_windows_summary

#### Change RMR units from kcal_hr to kcal_day to match tee units####
lowest_windows_summary <- lowest_windows_summary %>%
  rename(RMR_kcal_hr = avg_30min_value) %>%
  mutate(RMR_kcal_day = RMR_kcal_hr*24) %>%
  group_by(ID, SABLE)

##### Process sable_dwn into sable_TEE_data to get Avg daily TEE (tee) for each mouse at each time point####
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"))

sable_TEE_data <- sable_TEE_data %>%
  mutate(SABLE = factor(SABLE, 
            levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))


#### Attach echoMRI_data to sable_TEE_data --> name new df as sable_TEE_adj ####
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  dplyr::select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
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

# Left join Lean, Fat, and Weight info into TEE dataset
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% 
      dplyr::select(ID, SABLE, Lean, Weight, Fat),
    by = c("ID", "SABLE"))

##### Combine lowest_windows_summary with sable_TEE_adj ####
sable_TEE_adj_RMR <- sable_TEE_adj %>%
  left_join(
    lowest_windows_summary %>% 
      dplyr::select(ID, SABLE, window_start_time, window_end_time, RMR_kcal_day),
    by = c("ID", "SABLE")) %>%
  group_by(ID, SABLE) %>%
  mutate(TEE_minus_RMR = tee - RMR_kcal_day)

# End of original code to get RMR and tee data frame ####

#Energy balance calculations ####

# 1a. Combine EE data and FI during Sable (compiled FI, not just manually measured FI) ####
    #sable_TEE_adj_RMR and FI_in_SABLE. Calculate energy balance for each mouse.
    #Calculate average energy balance for each GROUP at each SABLE time point

#Remove mice with problematic cages from FI_in_SABLE
  #These mice were already removed from EE data frame when sable_TEE_adj_RMR was created
FI_in_Sable_balance <- FI_in_Sable %>%
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725))
  
#Left join to combine data frames
#Calculate energy balance for each mouse
EE_FI_in_Sable <- sable_TEE_adj_RMR %>%
  left_join(
    FI_in_Sable_balance %>% 
      select(ID, SABLE, GROUP, avg_corrected_intake_kcal),
    by = c("ID", "GROUP", "SABLE")) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  filter(!(SABLE == "BW regain" & ID == "3720")) %>% #remove visually apparent FI outliers. 3713 at Baseline (kcal=72.2) and 3720 at BW regain (kcal=53.6)
  filter(!(ID == "3713" & SABLE == "Baseline")) %>%
  mutate(tee_balance = avg_corrected_intake_kcal - tee) %>%
  ungroup()

#calculate energy balance for each GROUP at each SABLE time point
GROUP_balance <- EE_FI_in_Sable %>%
  group_by(SABLE, GROUP) %>%
  summarise (Energy_balance_kcal = mean(tee_balance), .groups = "drop") 

#graph energy balance using compiled FI
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

#Barplot of energy balance using COMPILED FI (combo of manual and Sable recorded)
ggplot(EE_FI_in_Sable, aes(x = SABLE, y = tee_balance, fill = GROUP)) +
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
    title= "Energy balance: Daily food intake (kcal) vs total energy expenditure (kcal) during Sable recording",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group") +
  theme_minimal() +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))

# FI_in_Sable statistical analyses: energy balance (Compiled Sable measured/manual FI) ####
#ANOVA

#Linear mixed model
#model_balance <- lmer(tee_balance ~ SABLE * GROUP + Weight + (1 | ID), data = EE_FI_in_Sable)
model_balance <- lmer(tee_balance ~ SABLE * GROUP + (1 | ID), data = EE_FI_in_Sable)
summary(model_balance)

n_distinct(EE_FI_in_Sable$ID) #good we have 16 animals

#### calculate estimated marginal means (emmeans) ####
emm_balance <- emmeans(model_balance, ~ SABLE * GROUP, cov.reduce = mean)
emm_balance_df <- as.data.frame(emm_balance)

# Pairwise contrasts within each GROUP
contrasts_by_group_balance <- contrast(emm_balance, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_balance_df <- as.data.frame(contrasts_by_group_balance)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_balance <- contrast(emm_balance, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_balance_df <- as.data.frame(contrasts_by_SABLE_balance)

#### Bar plot - Graph predicted energy balance ####
# Define custom colors
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

barplot_emm_balance <- emm_balance_df %>%
  #filter(SABLE %in% c("Baseline", "BW regain")) %>% #only includes EE gap time points
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
    title = "Energy balance (linear mixed model)",
    y = "Energy balance (kcal)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_balance

#---------------------------------#
# 1b. Combine EE data and FI during Sable (just manually measured FI) ####
    #sable_TEE_adj_RMR and Manual_FI_all_mice
#Remove mice with problematic cages from FI_in_SABLE
#These mice were already removed from EE data frame when sable_TEE_adj_RMR was created
Manual_FI_all_mice_balance <- Manual_FI_all_mice %>%
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725))

#Left join to combine data frames
#Calculate energy balance for each mouse
EE_Manual_FI_all_mice <- sable_TEE_adj_RMR %>%
  left_join(
    Manual_FI_all_mice_balance %>% 
      select(ID, SABLE, GROUP, avg_corrected_intake_kcal),
    by = c("ID", "GROUP", "SABLE")) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  mutate(tee_balance = avg_corrected_intake_kcal - tee) %>%
  ungroup()

#calculate energy balance (using manual FI) for each GROUP at each SABLE time point
Manual_GROUP_balance <- EE_Manual_FI_all_mice %>%
  group_by(SABLE, GROUP) %>%
  summarise (Energy_balance_kcal = mean(tee_balance), .groups = "drop") 

#Barplot of energy balance using Manual FI for all mice
ggplot(EE_Manual_FI_all_mice, aes(x = SABLE, y = tee_balance, fill = GROUP)) +
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
    title= "Energy balance: Manual daily FI (kcal) vs daily TEE (kcal) during Sable recording",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group") +
  theme_minimal() +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))

# FI_in_Sable statistical analyses: energy balance (Manual FI for all mice) ####

#UPDATE: Need to do ANOVA too ####

#Linear mixed model (multiple linear regression)
#model_balance <- lmer(tee_balance ~ SABLE * GROUP + Weight + (1 | ID), data = EE_FI_in_Sable)
model_balance_manual <- lmer(tee_balance ~ SABLE * GROUP + (1 | ID), data = EE_Manual_FI_all_mice)
summary(model_balance_manual)

n_distinct(EE_Manual_FI_all_mice$ID) #good we have 16 animals

#### calculate estimated marginal means (emmeans) ####
emm_balance_manual <- emmeans(model_balance_manual, ~ SABLE * GROUP, cov.reduce = mean)
emm_balance_manual_df <- as.data.frame(emm_balance_manual)

# Pairwise contrasts within each GROUP
contrasts_by_group_balance_manual <- contrast(emm_balance_manual, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_balance_manual_df <- as.data.frame(contrasts_by_group_balance_manual)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_balance_manual <- contrast(emm_balance_manual, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_balance_manual_df <- as.data.frame(contrasts_by_SABLE_balance_manual)

#### Bar plot - Graph predicted energy balance ####
# Define custom colors
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

barplot_emm_balance_manual <- emm_balance_manual_df %>%
  #filter(SABLE %in% c("Baseline", "BW regain")) %>% #only includes EE gap time points
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
    title = "Energy balance (linear mixed model) manual FI",
    y = "Energy balance (kcal)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_balance_manual

# Left off here on 2-5-26 #### (file is complete for the time being)
# Viewed again on 2-13-26. I should check for outliers in energy balance



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
