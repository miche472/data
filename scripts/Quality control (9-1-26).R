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

##Test PedMeters (aka movement in XY plane) ####
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

#Read in manually measured BW and FI
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19)
FI_LM_COHORT19 <- read_csv("~/Documents/GitHub/data/data/FI_LM.csv") %>%
  filter(COHORT == 19)

## Test BodyMass ####
BW_sable <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("BodyMass_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Body_Mass = value) %>%
  rename(parameter_Body_Mass = parameter) %>%
  rename(fix_value_Body_Mass = fix_value) 
  
BW_sable_summary <- BW_sable %>%
  filter(date > "2026-08-25")  %>%
  group_by(ID, date) %>%
  arrange(DateTime) %>%
  filter(Body_Mass >0) %>%
  summarise(BW_avg = mean(Body_Mass),
            BW_median= median(Body_Mass)) %>%
  rename("DATE"="date")  %>%
  mutate(ID = as.factor(ID))


#Manual BW
BW_manual <- BW_COHORT19 %>%
  filter(DATE > "2026-08-25") %>%
  mutate(ID = as.factor(ID))

#combine manual and Sable BW
BW_combined <- BW_manual %>%
  left_join(
    BW_sable_summary %>% 
      select(BW_avg, BW_median, DATE, ID), 
    by = c("ID", "DATE"))
  
#Combine sable and manual BW measurements
BW_combined <-BW_sable_summary  %>%
  left_join(
    BW_manual%>% 
      select(ID, BW, DATE, COMMENTS), 
    by = c("ID", "DATE"))
  
## Test FoodA ####
FI_sable <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("FoodA_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Food_G = value) %>%
  rename(parameter_Food_G = parameter) %>%
  rename(fix_value_Food_G = fix_value) %>%
  mutate(ID = as.factor(ID))

FI_sable2 <- FI_sable %>%
  filter(DATE > "2026-08-25") %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(FI = lag(Food_G)-Food_G)

#Manual FI
FI_manual <- FI_LM_COHORT19 %>%
  filter(DATE > "2026-08-25") %>%
  mutate(ID = as.factor(ID))

#---------------------------------------------------.#
#---------------------------------------------------.#
# Part 3: Compare set up 1 and set up 2 ####
#compare same mice in different cages: 
#Complete days created: incomplete (Fri), Fri->Sat, Sat->Sun, Sun->Mon, incomplete (Mon->Tues)
#Mice changed cages on sunday 8-30-2026 at 4:50pm. META.csv indicates the cage
  #they were in from 4:50pm to midnight on 8-30-26. Cage # and ID will be inaccurate
  #for 12:01am to 4:50pm on 8-30-26, so don't use these data.

    #3731: Cage 1 -> 5
    #3732: Cage 2 -> 6
    #3733: Cage 3 -> 7
    #3735: Cage 4 -> 8
    #3737: Cage 5 -> 1
    #3738: Cage 6 -> 2
    #3739: Cage 7 -> 3
    #3741: Cage 8 -> 4

#---------------------------------------------------.#
# Set up 1: Use  6pm on Friday (2026-08-28) to 6pm on Saturday (2026-08-29)
# Set up 2: Use 6pm on Sunday (2026-08-30) to 6pm on Monday (2026-08-31)

zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

## TEE from sable_dwn
TEE_quality <-sable_dwn_19 %>%
  filter(date %in% c("2026-08-28", "2026-08-29", "2026-08-30", "2026-08-31")) %>% 
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
         lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  # Identify which complete_days contain the dates of interest
  group_by(ID, complete_days) %>%
  mutate(
    has_aug28 = any(date == as.Date("2026-08-28")),
    has_aug29 = any(date == as.Date("2026-08-29")),
    has_aug30 = any(date == as.Date("2026-08-30")),
    has_aug31 = any(date == as.Date("2026-08-31"))) %>%
  mutate(
    SET_UP = case_when( # Assign setup based on which two dates occur within the complete day
      has_aug28 & has_aug29 ~ "one",
      has_aug30 & has_aug31 ~ "two",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(SET_UP)) %>%  # Keep ONLY the two complete days of interest
  #group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),min(DateTime), units="hours")), zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)


#---
#All movement (AllMeters) from sable_dwn 
loc_quality <-sable_dwn_19 %>%
 filter(date %in% c("2026-08-28", "2026-08-29", "2026-08-30", "2026-08-31")) %>% 
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
         lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  # Identify which complete_days contain the dates of interest
  group_by(ID, complete_days) %>%
  mutate(
    has_aug28 = any(date == as.Date("2026-08-28")),
    has_aug29 = any(date == as.Date("2026-08-29")),
    has_aug30 = any(date == as.Date("2026-08-30")),
    has_aug31 = any(date == as.Date("2026-08-31"))) %>%
  mutate(
    SET_UP = case_when( # Assign setup based on which two dates occur within the complete day
      has_aug28 & has_aug29 ~ "one",
      has_aug30 & has_aug31 ~ "two",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(SET_UP)) %>%  # Keep ONLY the two complete days of interest
  #group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),min(DateTime), units="hours")), zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 


#---
#Ambulation (PedMeters)
ped_quality <-sable_dwn_19 %>%
  filter(date %in% c("2026-08-28", "2026-08-29", "2026-08-30", "2026-08-31")) %>% 
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
         lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("PedMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  # Identify which complete_days contain the dates of interest
  group_by(ID, complete_days) %>%
  mutate(
    has_aug28 = any(date == as.Date("2026-08-28")),
    has_aug29 = any(date == as.Date("2026-08-29")),
    has_aug30 = any(date == as.Date("2026-08-30")),
    has_aug31 = any(date == as.Date("2026-08-31"))) %>%
  mutate(
    SET_UP = case_when( # Assign setup based on which two dates occur within the complete day
      has_aug28 & has_aug29 ~ "one",
      has_aug30 & has_aug31 ~ "two",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(SET_UP)) %>%  # Keep ONLY the two complete days of interest
  #group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),min(DateTime), units="hours")), zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

#Join EE and locomotion data 
Quality_loc_EE <- loc_quality %>%
  left_join(
    TEE_quality %>% 
      select(Kcal_Hr, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  left_join(
    ped_quality %>% 
      select(Ped_meters, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(ID=as.factor(ID)) %>%
  filter(!(ID==3740)) %>%
  mutate(GA_group = if_else(ID %in% c(3731, 3732, 3733, 3735), "group_1", "group_2"))

##EE Calculations #### 

Quality_EE_60min <- Quality_loc_EE %>%
  group_by(ID, DateTime) %>% 
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60)

###Total daily EE ####
EE_sum <- Quality_EE_60min %>%
  ungroup() %>%
  group_by(ID, complete_days) %>%
  arrange(DateTime) %>%
  summarise(n_obs = n(), #confirm that the same number of observations are being summed for each of the set ups
            sum_day_EE = sum(TEE_per_min)) %>%
  mutate(GA_group = if_else(ID %in% c(3731, 3732, 3733, 3735), "group_1", "group_2")) %>%
  mutate(complete_days=as.factor(complete_days)) 

####Plot x= GA_group, y= sum_day_EE, GROUP = set_up, color=set_up ####
ggplot(
  EE_sum,
  aes(x = complete_days,
    y = sum_day_EE,
    group = ID,
    color = ID)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_discrete(
  labels = c(
    "1" = "Setup 1",
    "3" = "Setup 2"
  )
) +
  labs(x = "Cage assignment",
    y = "Total energy expenditure (kcal/day)",
    color = "Mouse ID") +
  theme_classic() 


### EE: Photo period ####
#Avg EE and median EE during each photoperiod and 
Hour_avg_EE <- Quality_loc_EE %>%
  ungroup() %>%
  group_by(ID, SET_UP, lights) %>%
  summarise(avg_EE_lights = mean(Kcal_Hr),
            median_EE_lights = median(Kcal_Hr))

#---------------------------------------------------------------------------.
#---------------------------------------------------------------------------.
  
# Creat manual BW and FI dataframes ####
cohort_csv_files <- tibble(
  filepath = list.files("../data", full.names = TRUE)) %>% 
  filter(
    grepl("COHORT_[0-9]+[0-9]*.csv", filepath)) #now we can used cohort > 10
cohort_csv_files

cohort_open_files <- cohort_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map_dfr(
    ., function(X){
      read_csv(X$filepath) %>% 
        select(ID, FOOD_WEIGHT_START_G, FOOD_WEIGHT_END_G, DATE, DIET, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS) %>% 
        mutate(
          INTAKE_GR = (FOOD_WEIGHT_START_G - FOOD_WEIGHT_END_G),
          DATE = lubridate::mdy(DATE)
        ) %>% 
        select(ID, INTAKE_GR, DATE, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS) %>% 
        rename(
          BW = BODY_WEIGHT_G
        ) %>% 
        mutate(BW=as.numeric(BW), ID=as.factor(ID))})

# load food description
food_desc <- read_csv("../data/food_description.csv")

# load metadata
metadata <- read_csv("../data/META.csv") %>% 
  select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# output food-intake file
FI_LM <- cohort_open_files %>%
  select(ID, DIET_FORMULA, INTAKE_GR, DATE, COMMENTS) %>%
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    delta_alt = {
      intake_idx <- !is.na(INTAKE_GR) #creates a logical vector (TRUE/FALSE) where rows that have INTAKE_GR=NA --> FALSE and rows with a value for INTAKE_GR -->TRUE
      intake_dates <- DATE[intake_idx] #Unconfirmed: only keeps rows for which intake_idx is TRUE
      
      # compute differences only on valid intake rows
      diffs <- c(NA, as.numeric(diff(intake_dates)))
      
      # create full-length vector and fill only intake rows
      out <- rep(NA_real_, n())
      out[intake_idx] <- diffs
      out
    }
  ) %>%
  mutate(delta_measurement = DATE - lag(DATE)) %>% #just use to remove first observation for each mouse
  drop_na(delta_measurement) %>% #just use to remove first observation for each mouse
  mutate(corrected_intake_gr = INTAKE_GR / as.numeric(delta_alt)) %>%
  left_join(., food_desc, by = "DIET_FORMULA") %>%
  mutate(corrected_intake_kcal = corrected_intake_gr * KCAL_G) %>%
  left_join(., metadata, by = "ID")  %>%
  select(-delta_measurement)


# output bodyweight file
BW <- cohort_open_files %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, BW, DATE,COMMENTS) %>% 
  drop_na(BW) %>% 
  left_join(., metadata, by = "ID")

write_csv(x = FI_LM, "../data/FI_LM.csv")
write_csv(x = BW, "../data/BW.csv")

#Read in BW and FI
#Read in BW and filter for cohort 19 (Spring 2026 NZO mice)
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19)

#Read in FI and filter for cohort 19 (Spring 2026 NZO mice)
FI_LM_COHORT19 <- read_csv("~/Documents/GitHub/data/data/FI_LM.csv") %>%
  filter(COHORT == 19)

#Create df with BW and FI
BW_FI_19 <- BW_COHORT19 %>% #Join FI and BW
  left_join(
    FI_LM_COHORT19 %>% 
      select(ID, INTAKE_GR, DATE, delta_alt, corrected_intake_gr, corrected_intake_kcal, KCAL_G),
    by = c("ID", "DATE")) %>%
  mutate(ID = as.factor(ID)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(day_rel = DATE - first(DATE),
         day_rel = as.numeric(day_rel))
replace_na(list(#INTAKE_GR=0, 
  #delta_alt=0, 
  #corrected_intake_gr=0,
  #corrected_intake_kcal=0, 
  KCAL_G=3.82))



