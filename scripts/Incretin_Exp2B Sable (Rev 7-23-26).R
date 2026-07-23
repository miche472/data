# Energy expenditure in GLP-1 experiment 2B (10 -> 15nmol/kg)

#Started: 7-9-26
#Updated:7-23-26


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


#Functions####
##USE -- Updated: zt_time ####
#This method gives integer values for zt_time and accounts for the 6:30 to 18:30 cycle
#It is more explicit in its format/logic compared to original zt_time function
zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

#Identify when mice moved using AllMeters
filter_loc1_exp2b <-sable_dwn %>%
  filter(COHORT ==19) %>%
  filter(ID %in% c(3744,3745,3746,3748,3752)) %>%
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      
      sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss_ctrl",
      
      sable_idx %in% c("SABLE_DAY_9", "SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12",
                       "SABLE_DAY_13","SABLE_DAY_14")~"BW loss_TZP")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(zt_time = zt_time(DateTime),
    is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
    complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, SABLE, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                             min(DateTime),
                                             units="hours")),
    zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
ungroup()

  
#Repeat for EE data
filter_TEE1_exp2b <-sable_dwn %>%
  filter(COHORT ==19) %>%
  filter(ID %in% c(3744,3745,3746,3748,3752)) %>%
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      
      sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss_ctrl",
      
      sable_idx %in% c("SABLE_DAY_9", "SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12",
                       "SABLE_DAY_13","SABLE_DAY_14")~"BW loss_TZP")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(zt_time = zt_time(DateTime),
    is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
    complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, SABLE, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                             min(DateTime),
                                             units="hours")),
    zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()


#---
#In df filter_loc2 use mutate to make a column called AllMeters_ using the value column data
#In df filter_TEE2 use mutate to make a column called kcal_hr_ using the value column data

filter_loc2_exp2b <- filter_loc1_exp2b %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_TEE2_exp2b <- filter_TEE1_exp2b %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_loc2 and filter_TEE2 into a df called filter_loc_TEE2
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_loc_TEE2_exp2b <- filter_loc2_exp2b %>%
  left_join(
    filter_TEE2_exp2b %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx? ####
    by = c("ID", "DateTime", "sable_idx")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss_ctrl","BW loss_TZP"))) 
 #filter(!ID=="3748")
  #filter(!ID %in% c(3748, 3749, 3751)) #issues with recording

#On some days mouse A was recording in cage I in the morning and mouse B starting recording in cage I
#in the evening. META.csv only indicates dates not hours, so I need to extricate the data from these two mice
#in the code

new_filter_loc_TEE2_exp2b <- filter_loc_TEE2_exp2b %>%
  group_by(ID, SABLE) %>%
  filter(max(complete_days) <= 2 |complete_days != max(complete_days)) %>% 
  ungroup()
#Now I should have accurate data for days when two mice record in the same cage number


#---------------------------------------------------------------------
#To look at just one complete day at a time
filter_loc_TEE3_exp2b <- new_filter_loc_TEE2_exp2b %>%
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

#---------------------------------------------------------------------
#To calculate TEE for each complete day rather than just the first complete day. 
#Then, for each mouse, find the average daily EE across the last two sable days

Minute_sum_EE_hr_all_days_exp2b <- new_filter_loc_TEE2_exp2b %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, complete_days) %>%
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

#For each complete day, Calculate daily sum (add all hours together)
Daily_EE_all_days_exp2b <- Minute_sum_EE_hr_all_days_exp2b %>%
  ungroup() %>%
  group_by(ID, SABLE, complete_days) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))

#Was there a big difference in EE across days within one sable time point?
#Manually look at this using mean and SD --> SD is relatively small
Avg_daily_EE_exp2B <- Daily_EE_all_days_exp2b %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(Avg_TEE = mean(TEE_kcal_day),
            SD_TEE = sd(TEE_kcal_day))

#Should start with a linear mixed model to see if there is significant day to day variation
#during any of the phases. If there isn't then averaging complete days within a phase would be fine

model_complete_days_exp2b <- lmer(TEE_kcal_day ~ SABLE*complete_days + (1 | ID), data = Daily_EE_all_days_exp2b)
summary(model_complete_days_exp2b)

model_complete_days2_exp2b <- lmer(TEE_kcal_day ~ SABLE + complete_days + (1 | ID), data = Daily_EE_all_days_exp2b)
summary(model_complete_days2_exp2b)
