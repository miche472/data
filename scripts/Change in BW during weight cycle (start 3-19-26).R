#Script: Change in BW during weight cycle ####

#Script was started with code from Energy efficiency & %BW (Rev. 3-10-26)
#Started: 3-19-26
#Last revised: 3-19-26

  
#Aim 1: Compare changes in BW between Sable time points (i.e. during stages of weight cycle)
        #Graph raw, create linear mixed model and calculate 
        #estimated marginal means, Graph predictions from model
#Aim 2: Compare rate of change over time for BW between Sable time points

#Objectives ####
  # % BW change from starting BW
  # All measurements of BW over the course of the study
  # STAGES of weight cycle
  # Calculate BW & duration (& FI) for each stage of weight cycle
  #Stat analysis: ΔBW by stage


#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(patchwork)

BW <- read_csv("../data/BW.csv")
FI <- read_csv("../data/FI.csv")
#-------------------------------------------------------#
# % BW change from starting BW ####
#percent change in BW (relative to starting BW...BW on Sable day 7) at each subsequent BW measurement
  # Calculated for each mouse ID at every time point, not just SABLE sign posts

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(COMMENTS_BW = COMMENTS) %>%
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"),
    day_rel = DATE - first(DATE)) %>%
  mutate(STATE = case_when(
      ##Baseline: First day of LFD/first day of obesity development
      #Date is the first date after start of LFD for which there was a BW and FI measurement
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
      ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
      ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
      
      ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
      #First day of calorie restriction for Weight Cycled mice. 
      #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
      ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
      
      ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
      #Sable day 15
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
      
      ##BW maintenance: End of Sable recording for BW maintenance
      #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
      ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
      ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
      
      ##BW regain: Final day of experiment --> day of sac (End of regain)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!(ID == 3726 & DATE == as.Date("2025-04-28"))) %>%  #repeated
  mutate(STATE = factor(STATE, # Make STATE an ordered factor
                        levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#-----------------------------------------------------------------------------#
# All measurements of BW and FI over the course of the study ####

#Prepare a df to look at every measurement of FI and BW (df_all_BW_FI)
#The data frame called df_all_BW_FI in this script is based on a df called df1 in "FI & BW (started 2-16-26).R"
#To create df_all_BW_FI, first run the chunk of code below which combines FI.csv and BW.csv to create FI_BW_joined
#df_all_BW_FI retains measurements from all days 
#df1 only retains measurements from sign post days
df_all_BW_FI <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW)) %>% 
  mutate(day_rel = DATE - first(DATE),
         FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
         STATE = case_when(
           ##Baseline: First day of LFD/first day of obesity development
           #Date is the first date after start of LFD for which there was a BW and FI measurement
           ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
           ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
           ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
           ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
           
           ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
           #First day of calorie restriction for Weight Cycled mice. 
           #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
           ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                     3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                     3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
           
           ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
           #Sable day 15
           ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
           ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           
           ##BW maintenance: End of Sable recording for BW maintenance
           #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           
           ##BW regain: Final day of experiment --> day of sac (End of regain)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  #in this script we want to retain data from all dates, not just phase "sign post" dates
  mutate(STATE = factor(STATE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 

#-------------------------------------------------------#
# STAGES of weight cycle ####
#Prepare data frames to calculate daily FI and delta BW 
##Create df1 
#Create df1 directly in this script (originally created in "FI & BW (Started 2-16-26).R")
#Prepare BW.csv
BW_to_join <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>%
  rename(COMMENTS_BW = COMMENTS)

#Prepare FI.csv
FI_to_join <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% #remove time when fed chow
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>%
  rename(COMMENTS_FI = COMMENTS)

##Create df = FI_BW_joined 
#Join prepared BW and FI data frames
FI_BW_joined <- FI_to_join %>%
  left_join(
    BW_to_join %>% 
      select(ID, DATE, BW, COMMENTS_BW),
    by = c("ID", "DATE"))

# Create df1 
#Adds variables: GROUP, DRUG, STATE, day_rel, FI_rel, FI_cum to joined BW & FI data
#Removes measurements that aren't from sign post dates (i.e. STATE)
df1 <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(day_rel = DATE - first(DATE),
         FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
         STATE = case_when(
           ##Baseline: First day of LFD/first day of obesity development
           #Date is the first date after start of LFD for which there was a BW and FI measurement
           ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
           ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
           ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
           ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
           
           ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
           #First day of calorie restriction for Weight Cycled mice. 
           #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
           ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                     3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                     3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
           
           ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
           #Sable day 15
           ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
           ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           
           ##BW maintenance: End of Sable recording for BW maintenance
           #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           
           ##BW regain: Final day of experiment --> day of sac (End of regain)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(STATE)) %>%   #Remove measurements that aren't from sign post dates (i.e. STATE)
  mutate(STATE = factor(STATE,levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 

#---#
#Calculate FI, BW, & duration ####
## FI by stage ####
#Use df1 to calculate cumulative FI for each transition period 
# Summarize cumulative FI per ID and STATE
FI_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(FI_cum_end = max(FI_cum_INTAKE_kcal, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = FI_cum_end) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `Peak obesity` - Baseline,
    kcal_peak_to_loss = `BW loss` - `Peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
FI_stage_long <- FI_stage_summary %>%
  select(ID, DRUG, GROUP,
         kcal_baseline_to_peak,
         kcal_peak_to_loss,
         kcal_loss_to_maint,
         kcal_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("kcal_"),
    names_to = "Transition",
    values_to = "kcal") %>%
  mutate(
    Transition = factor(
      Transition,
      levels = c("kcal_baseline_to_peak",
                 "kcal_peak_to_loss",
                 "kcal_loss_to_maint",
                 "kcal_maint_to_regain"),
      labels = c("Obesity development",         #formerly Baseline → Peak obesity
                 "Weight loss",                 #formerly Peak obesity → BW loss
                 "Weight maintenance",          #formerly BW loss → BW maintenance
                 "Weight regain")))             #formerly BW maintenance → BW regain

##Duration by stage ####
#use df1 to calculate duration (days) for each stage 
# Summarize duration (days) per ID and between STATE
Days_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(max_day_rel = max(day_rel, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = max_day_rel) %>%
  # Calculate days between stages
  mutate(
    Days_baseline_to_peak = `Peak obesity` - Baseline,
    Days_peak_to_loss = `BW loss` - `Peak obesity`,
    Days_loss_to_maint = `BW maintenance` - `BW loss`,
    Days_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
Days_stage_long <- Days_stage_summary %>%
  select(ID, DRUG, GROUP,
         Days_baseline_to_peak,
         Days_peak_to_loss,
         Days_loss_to_maint,
         Days_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("Days_"),
    names_to = "Transition",
    values_to = "Duration_day") %>%
  mutate(
    Transition = factor(
      Transition,
      levels = c("Days_baseline_to_peak",
                 "Days_peak_to_loss",
                 "Days_loss_to_maint",
                 "Days_maint_to_regain"),
      labels = c("Obesity development",    #formerly Baseline → Peak obesity
                 "Weight loss",            #formerly Peak obesity → BW loss
                 "Weight maintenance",     #formerly BW loss → BW maintenance
                 "Weight regain")))        #formerly BW maintenance → BW regain

# df2: Join duration (days) and cummulative FI during each phase
FI_duration_joined <- FI_stage_long %>%
  left_join(
    Days_stage_long %>% 
      select(Transition, ID, Duration_day),
    by = c("ID", "Transition"))

df2 <- FI_duration_joined %>%
  group_by(GROUP, Transition, ID) %>%
  mutate(Daily_kcal = kcal/as.numeric(Duration_day))

#Verify that there are 22 mice in each of the 4 transition periods 
df2 %>% 
  group_by(Transition) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO in each transition period

##ΔBW (g) by stage ####
#use df1 to calculate delta_BW (g) for each stage
# Summarize Δ BW (g) per ID and STATE
BW_stage_summary <- df1 %>% #df1 has 5 values for each mouse (one per time point)
  group_by(ID, GROUP, STATE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,         #time period (5 total)
    values_from = BW_end) %>%   #BW at end of each time period
  # Calculate change in BW between time points for each ID
  mutate(
    BWg_baseline_to_peak = `Peak obesity` - Baseline,
    BWg_peak_to_loss = `BW loss` - `Peak obesity`,
    BWg_loss_to_maint = `BW maintenance` - `BW loss`,
    BWg_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
BW_stage_long <- BW_stage_summary %>%
  select(ID, GROUP,
         BWg_baseline_to_peak,
         BWg_peak_to_loss,
         BWg_loss_to_maint,
         BWg_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("BWg_"),
    names_to = "Transition",       #name of new column which lists all transition periods for all mice
    values_to = "delta_BW_g") %>%  #name of new column with delta BW values for each transition period
  mutate(
    Transition = factor(
      Transition,
      levels = c("BWg_baseline_to_peak",
                 "BWg_peak_to_loss",
                 "BWg_loss_to_maint",
                 "BWg_maint_to_regain"),
      labels = c("Obesity development",    #formerly Baseline → Peak obesity
                 "Weight loss",            #formerly Peak obesity → BW loss
                 "Weight maintenance",     #formerly BW loss → BW maintenance
                 "Weight regain")))        #formerly BW maintenance → BW regain

#Combine cumulative FI, duration, ΔBW for each phase -> Add BW_stage_long to df2 
  #df2 has cumulative FI and # of days in each stage
FI_duration_BW_joined <- df2 %>%
  left_join(
    BW_stage_long %>% 
      select(Transition, ID, delta_BW_g),
    by = c("ID", "Transition"))

# Change in BW per kcal (phase_BW_FI) for each ID during each stage
# Change in BW per kcal is variable in df4 called --> phase_BW_FI
df4 <- FI_duration_BW_joined %>%
  group_by(GROUP, Transition, ID) %>%
  mutate(phase_BW_FI = delta_BW_g/kcal)

df4 %>% #Verify that there are 22 mice in each of the 4 transition periods 
  group_by(Transition) %>%
  summarise(n_ID = n_distinct(ID)) #we have 22 NZO in each transition period

## (use) Graph ΔBW (g) by stage ####
#df4 has delta_BW_grams which is the change in BW for each ID during each transition period
# (use) Graph ΔBW (g) for Control and Weight cycled mice
ggplot(df4, aes(x = Transition, y = delta_BW_g, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
        axis.title = element_text(face = "bold"))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  labs(
    title="Change in body weight (g)",
    x = "Stage of weight cycle",
    y = "Change in BW (g) during stage",
    fill = "Treatment group")
#This shows change in BW in grams, but percent change in BW may be more useful since there was a wide range in BW

##Percent change in BW during each stage ####
#Need to group by ID and state and then do--> bw_rel = 100 * (BW - first(BW)) / first(BW)
#Then filter out measurements that don't have a STATE?...would this give me % change 
  #between start and end of a transition period?

#Method for percent change in BW
  #Modify the code that creates BW_stage_summary.(safer approach than using lag(BW))
  #Modification: calculate % change rather than change in grams between phases
  #ΔBW during each stage: use df1 to calculate delta_BW (g) for each stage

# Summarize Δ BW (g) per ID and STATE
BW_stage_summary_percent <- df1 %>% #df1 has 5 values for each mouse (one per time point)
  group_by(ID, GROUP, STATE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,         #time period (5 total)
    values_from = BW_end) %>%   #BW at end of each time period
  # Calculate change in BW between time points for each ID
  group_by(ID) %>%
  mutate(
    BWp_baseline_to_peak = 100*((`Peak obesity` - Baseline)/Baseline), 
    BWp_peak_to_loss = 100*((`BW loss` - `Peak obesity`)/`Peak obesity`),
    BWp_loss_to_maint = 100*((`BW maintenance` - `BW loss`)/`BW loss`),
    BWp_maint_to_regain = 100*((`BW regain` - `BW maintenance`)/`BW maintenance`))

# Convert to long format
BW_stage_long_percent <- BW_stage_summary_percent %>%
  select(ID, GROUP,
         BWp_baseline_to_peak,
         BWp_peak_to_loss,
         BWp_loss_to_maint,
         BWp_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("BWp_"),
    names_to = "Transition",       #name of new column which lists all transition periods for all mice
    values_to = "delta_BW_perce") %>%  #name of new column with delta BW values for each transition period
  mutate(
    Transition = factor(
      Transition,
      levels = c("BWp_baseline_to_peak",
                 "BWp_peak_to_loss",
                 "BWp_loss_to_maint",
                 "BWp_maint_to_regain"),
      labels = c("Obesity development",    #formerly Baseline → Peak obesity
                 "Weight loss",            #formerly Peak obesity → BW loss
                 "Weight maintenance",     #formerly BW loss → BW maintenance
                 "Weight regain")))        #formerly BW maintenance → BW regain

# df5: Add BW_stage_long to df4. df4 has cumulative FI, # of days, and ΔBW (g) in each stage
df5 <- df4 %>%
  left_join(
    BW_stage_long_percent %>% 
      select(Transition, ID, delta_BW_perce),
    by = c("ID", "Transition"))

## (use) Graph % change in BW for Control and Weight cycled mice ####
ggplot(df5, aes(x = Transition, y = delta_BW_perce, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  #format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
        axis.title = element_text(face = "bold"))+
  labs(
    title="Change in body weight (%)",
    x = "Stage of weight cycle", 
    y = "Change in BW (%) during stage",
    fill = "Treatment group")

#Stat analysis: ΔBW by stage ####
# df5 has change in BW during each stage for each mouse in terms of grams and %

#Analysis for change in BW in grams
###-----LMM: BW change (g) during stages of weight cycle -----####
#Build multiple linear regression model
model_change_BW_g <- lmer(delta_BW_g ~ Transition * GROUP + (1 | ID), data = df5)
summary(model_change_BW_g)

# Calculate estimated marginal means #
emm_change_BW_g <- emmeans(model_change_BW_g, ~ Transition * GROUP, cov.reduce = mean)
emm_change_BW_g_df <- as.data.frame(emm_change_BW_g)

# Pairwise contrasts within each GROUP
contrasts_by_group_change_BW_g <- contrast(emm_change_BW_g, method = "pairwise", by = "GROUP")
contrasts_by_group_change_BW_g_df <- as.data.frame(contrasts_by_group_change_BW_g)

# Pairwise contrasts within each stage of weight cycling
contrasts_by_SABLE_change_BW_g <- contrast(emm_change_BW_g, method = "pairwise", by = "Transition")
contrasts_SABLE_change_BW_g_df <- as.data.frame(contrasts_by_SABLE_change_BW_g)

# Bar plot - Graph predicted changed in BW (grams) #
barplot_emm_change_BW_g <- emm_change_BW_g_df %>%
  ggplot(aes(x = Transition, y = emmean, fill = GROUP)) +
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
    title = "MLR: Change in BW (grams) unadjusted",
    y = "Change in BW (grams)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold", size=15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11))
barplot_emm_change_BW_g
#Conclusion: The change in BW (grams) during BW regain stage is significantly larger 
    #for weight cycled vs control (linear mixed model, estimated marginal means, p=3.86e-9)

#Analysis for percent change in BW
###-----LMM: BW change (%) during stages of weight cycle -----####
#unadjusted change in BW (%) 
#Build multiple linear regression model
model_delta_BW_perce <- lmer(delta_BW_perce ~ Transition * GROUP + (1 | ID), data = df5)
summary(model_delta_BW_perce)

# Calculate estimated marginal means 
emm_delta_BW_perce <- emmeans(model_delta_BW_perce, ~ Transition * GROUP, cov.reduce = mean)
emm_delta_BW_perce_df <- as.data.frame(emm_delta_BW_perce)

# Pairwise contrasts within each GROUP
contrasts_by_group_delta_BW_perce <- contrast(emm_delta_BW_perce, method = "pairwise", by = "GROUP")
contrasts_by_group_delta_BW_perce_df <- as.data.frame(contrasts_by_group_delta_BW_perce)

# Pairwise contrasts within each stage of weight cycling
contrasts_by_SABLE_delta_BW_perce <- contrast(emm_delta_BW_perce, method = "pairwise", by = "Transition")
contrasts_SABLE_delta_BW_perce_df <- as.data.frame(contrasts_by_SABLE_delta_BW_perce)

# Bar plot - Graph predicted percent change in BW (emmeans for MLR) 
barplot_emm_delta_BW_perce <- emm_delta_BW_perce_df %>%
  ggplot(aes(x = Transition, y = emmean, fill = GROUP)) +
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
    title = "MLR: Percent change in BW, unadjusted",
    y = "Change in BW (%)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold", size=15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11))
barplot_emm_delta_BW_perce
#Conclusion: The % change in BW during BW regain stage is significantly greater for weight cycled  
    #compared to control mice (linear mixed model, estimated marginal means, p=1.2e-8)

#------#
# Aim 2: Change over time (slope) within stages ####
#Use iCAL time points as sign posts. synthesize individual measurements to assess how BW changes over time
#compare between stages for the same GROUP and for the same stage, but between GROUPS

#Label all observations with what stage they are part of by using fill()
df_stage <- BW_data %>%
  arrange(ID, DATE) %>%   # make sure time is ordered within mouse
  group_by(ID) %>%        # VERY important: do this per mouse
  fill(STATE, .direction = "up") %>%  # carry forward last SABLE
  ungroup() %>%
  rename(stage = STATE)

#Check if there are any NAs in the "stage" column (both methods below are useful)
df_stage %>% #see which rows (if any) have NAs
  filter(is.na(stage))

df_stage %>% #check per mouse
  group_by(ID) %>%
  summarise(n_na = sum(is.na(stage))) %>%
  print(n=22)

#Create a linear mixed model to calculate slopes
#First, create a time variable from DATE...kind of like CS's day_rel approach
df_stage <- df_stage %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(time = as.numeric(DATE - min(DATE)), time_sc = time / 7) #converts to weeks...change in BW weekly

#Build linear mixed model
model_df_stage <- lmer(BW ~ time_sc * stage * GROUP + (1 | ID), data = df_stage)
#summary(model_df_stage)

#Calculate estimated marginal slopes (analog of calculating estimated mariginal mean)
    #Use emtrends rather than emmeans because I am comparing slopes. emmeans is 
    #for single measurements (such as just BW at sign posts)
emm_trends <- emtrends(model_df_stage, ~ stage * GROUP,var = "time_sc")
emm_trends_df <- as.data.frame(emm_trends)

# Pairwise contrasts within each GROUP
contrasts_by_group <- contrast(emm_trends,method = "pairwise",by = "GROUP")
contrasts_by_group_df <- as.data.frame(contrasts_by_group)

# Pairwise contrasts within each stage
contrasts_by_stage <- contrast(emm_trends,method = "pairwise", by = "stage")
contrasts_by_stage_df <- as.data.frame(contrasts_by_stage)

#Left off --> graph this to make sure it is correct and that I understand 

#-------------------------------------------------#
#Food intake
#We did not measure food intake daily, so 
#if I estimate daily intake by dividing the amount of food eaten by the number 
#of days between measurements and then take the average of that extrapolation to get
#an average daily intake over the entire time period, I will have taken the average twice. 
#Additionally, mice ate significantly less on Sable acclimation days, days when they had 
#echoMRIs, and days when their cage was changed. 
#Another complication is that different mice were in each stage for a different number of days.
#If I calculate cumulative total daily intake during the duration of a stage and divide by the number 
#of days during that stage, that gives an estimate of daily intake which has not been manipulated 
#prior to conducting the calculation. 
#I could look at weekly food intake, but that would get very messy and be difficult to verify that the code is working correctly

#I could estimate by comparing the rate of increase in cumulative food intake, but that requires 
#an arbitrary choice of time unti over which the change occurs. I could try using weeks, but this could be 
#difficult if there was incomplete weeks

#I will use the daily intake calculated using my method of total cumulative intake/number of days in phase. 
#Then I will try using an ANOVA that accounts fo repeated measures, two treatment groups, multiple stages,
#and interaction between treatment group and stage
#I will also try using linear mixed model: 
#model <- lmer(daily FI ~ stage * GROUP + (1 | ID), data = ___)

