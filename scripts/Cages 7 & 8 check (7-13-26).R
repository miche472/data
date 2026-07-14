# Verify that cages 7 and 8 gas analyzers are working properly
#Recorded data for 3733 and 3741 in cages 7 and 8, respectively from 
#7/8 through 7/10

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
library(hms)


#-__
#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

#Identify when mice moved using AllMeters
filter_loc1_test<-sable_dwn %>%
  filter(COHORT ==19) %>%
  filter(ID %in% c(3733, 3741)) %>%
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
    is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
    complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                             min(DateTime),
                                             units="hours")),
    zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 23, #complete day must have a total of at least 23 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
ungroup()

  
#Repeat for EE data
filter_TEE1_test <-sable_dwn %>%
  filter(COHORT ==19) %>%
  filter(ID %in% c(3733, 3741)) %>%
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
    is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
    complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, omplete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                             min(DateTime),
                                             units="hours")),
    zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 23, #complete day must have a total of at least 23 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()


#---
#In df filter_loc2 use mutate to make a column called AllMeters_ using the value column data
#In df filter_TEE2 use mutate to make a column called kcal_hr_ using the value column data

filter_loc2_test <- filter_loc1_test %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_TEE2_test <- filter_TEE1_test %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_loc2 and filter_TEE2 into a df called filter_loc_TEE2
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_loc_TEE2_test <- filter_loc2_test %>%
  left_join(
    filter_TEE2_test %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx? ####
    by = c("ID", "DateTime", "sable_idx")) 

#---------------------------------------------------------------------
#To look at just one complete day at a time
filter_loc_TEE3_test <- filter_loc_TEE2_test %>%
  filter(complete_days ==1)

# RMR calculated using percentile Percentile for RMR --> 
Minute_sum_EE_hr_exp2b <- filter_loc_TEE3_exp2b %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr) %>%
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    # Total Energy Expenditure (TEE) --> minute by minute summation
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    # Resting metabolic rate for each hour calculated using 10th percentile method for RMR:
    #Note: this is the rate of RMR rate within each hr not for the entire 24hrs
    #Calculating this will allow for calculation of NEAT
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    # Total RMR energy across observed time in the hour
    RMR_kcal = RMR_rate * (n_obs / 60),
    # NEAT= TEE-RMR when the mouse is moving
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE), .groups = "drop") %>%
  
  # Verify that TEE = RMR + NEAT --> TEEvsRMR_NEAT should be close to zero
  mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal)) %>%
  #Indicate when mice received vehicle and when they received 
  mutate(Treatment_stage = case_when(
    SABLE=="Peak obesity"~ "No injections",
    SABLE== "BW loss_ctrl" ~ "Vehicle",
    SABLE=="BW loss_TZP" ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss_ctrl", "BW loss_TZP")))

#Calculate daily sum (add all hours together)
Daily_EE_exp2b <- Minute_sum_EE_hr_exp2b %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))
