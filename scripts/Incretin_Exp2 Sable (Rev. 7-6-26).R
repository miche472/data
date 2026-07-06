# Energy expenditure in GLP-1 experiment 2

#Started: 6-5-26
#Updated: 7-2-26

#Next week: do linear mixed model comparing every complete day. also graph the raw data including
#each individual measurement. For the stat analysis I need to remove mice that were measured in  
#non-functional cages at peak obesity

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

##Old zt_time ####
#zt_time <- function(hr){
  #return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))}

##Alternative updated zt_time: zt_time ####
#This method gives integer values for zt_time and accounts for the 6:30 to 18:30 cycle
#It mirrors the format that zt_time had prior to the change in the light cycle
#zt_time <- function(DateTime){
  #time <- hour(DateTime) + minute(DateTime)/60
  #return(if_else(time >= 18.5,floor(time - 18.5),floor(time + 5.5)))}

##zt_time for 7pm to 6am light cycle ####
#zt_time <- function(hr){
  #return(if_else(hr >= 19 & hr <= 23, hr-19, hr+5))}

# Start KEEP : replacement for above version of filter_loc1 ####
filter_loc1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss")) %>% 
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
  filter(recording_duration >= 23, #complete day must have a total of at least 23 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()
 
#figure out why some days have 47hrs duration
#filter_loc1_revised %>%
  #filter(ID == 3742, #issue for 3742. 3743, 3849, 3750, 3751, 3753
        # SABLE == "BW loss",
         #complete_days == 3) %>%
  #select(DateTime, zt_time, is_zt_init) %>%
  #View()

#Figure out if there is missing data
filter_loc1_revised %>% 
  #group_by(SABLE, ID, date) %>%
  #filter(SABLE=="Peak obesity")%>% 
  summarise(n_each_ID = n_distinct(hr)) %>%
  print(n=57) %>%
  filter(n_each_ID== 24) %>%
  print(n=42)
  

#Now create code to calculate locomotion on each complete_days (see below)

# end "KEEP this updated version" ####

# Start "revise this to match df filter_loc1_revised" ####

#Checks ####
#Figure out if there is missing data
filter_loc1_revised %>% 
  group_by(SABLE, ID, complete_days) %>%
  filter(SABLE=="BW loss")%>% 
  summarise(n_each_ID = n_distinct(hr)) %>%

#THis check could still be useful, but should be revised now that I am thinking about complete_days slightly differently
filter_loc1_revised %>% 
  group_by(SABLE, ID, hr, complete_days) %>%
  filter(complete_days ==1) %>%
  filter(SABLE=="BW loss") %>%
  filter(ID==3743) %>%
  summarise(n_each_ID = n_distinct(DateTime)) %>%
  print(n=24)

# Get just TEE 
filter_TEE1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
         lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
         SABLE= case_when(
           sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                            "SABLE_DAY_4") ~ "Peak obesity",
           sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                            "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss")) %>% 
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
  filter(recording_duration >= 23, #complete day must have a total of at least 23 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()


#---
#In df filter_loc2 use mutate to make a column called AllMeters_ using the value column data
#In df filter_TEE2 use mutate to make a column called kcal_hr_ using the value column data

filter_loc2 <- filter_loc1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_TEE2 <- filter_TEE1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_loc2 and filter_TEE2 into a df called filter_loc_TEE2
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_loc_TEE2 <- filter_loc2 %>%
  left_join(
    filter_TEE2 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx? ####
    by = c("ID", "DateTime", "sable_idx")) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss"))) %>%
  #filter(!ID %in% c(3748, 3749, 3751)) #issues with recording
  filter(!cage_number %in% c("5","6","8"))
  #try removing the cages that weren't working during Peak obesity rather than IDs that were in those cages

#OPTIONS ####
#To look at just one complete day at a time
filter_loc_TEE3 <- filter_loc_TEE2 %>%
  filter(complete_days ==1)

#To look at multiple complete days, distinguishing between SABLE phases
#(Group by mouse, SABLE phase, and complete day for downstream analyses)
filter_loc_TEE3_b <- filter_loc_TEE2 %>%
  group_by(ID, SABLE, complete_days)

#Taken from 4/16 script and modified
  # In previoius scripts I used whether or not the mouse moved between two minutes to classify
  #the EE during that minute as entirely NEAT or entirely RMR. This results in an 
  #overestimate of NEATsince even when I mouse is moving some of its TEE is RMR.
  #A way to get around this issue is to calculate the 10th percentile of TEE for a given 
  #observation period (hr, lights on/off, or daily [aka global])--> for now, do by hr
  #For min. when mouse didn't move: RMR=TEE and NEAT=0. 
  #For min. when mouse moved, NEAT = TEE - (RMR=10th percentile of RMR values across that hr for the ID)

#---------------------------------------------------------------------
# RMR calculated using percentile Percentile for RMR --> 
Minute_sum_EE_hr <- filter_loc_TEE3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, DRUG) %>%
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
  #Re-attach DRUG
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss")))

#Calculate daily sum (add all hours together)
Daily_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))


#EE, RMR, and NEAT broken down by light and dark cycle
#cumulative EE using the 10th percentile method for RMR 

#Photo periods: EE (cumulative)-> TEE, NEAT, and RMR ####
#distinguishing between light and dark photo periods
Photo_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG, lights) %>%
  summarise(TEE_kcal_photo = sum(TEE_kcal), 
            NEAT_kcal_photo = sum(NEAT_kcal),
            RMR_kcal_photo = sum(RMR_kcal))

  
##Try using zt_time rather than hr to calculate the calories expended each hour
#Conclusion: extremely similar values for TEE, RMR, and NEAT. In terms of the outcome
#I don't think zt_time or hr matters. I still don't know what is technically 
#correct. I susepct zt_time is more correct, but hr is more intuitive

# RMR calculated using percentile Percentile for RMR --> 
Minute_sum_EE_hr_zttime <- filter_loc_TEE3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, DRUG) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, zt_time) %>%
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
  #Re-attach DRUG
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss")))

#Calculate daily sum (add all hours together)
Daily_EE_zttime <- Minute_sum_EE_hr_zttime %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))


  
  
  
#Try looking at more than just complete_days ==1...look at all of the complete days ####

Minute_sum_EE_hr_zttime_alldays <- filter_loc_TEE2 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, complete_days) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, complete_days, zt_time) %>%
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
  #Re-attach DRUG
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss")))

#Calculate daily sum (add all hours together)
Daily_EE_zttime_alldays <- Minute_sum_EE_hr_zttime_alldays %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG, complete_days) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))


#---
#---
# Linear mixed models ####

#How can I distinguish between complete days within a SABLE time point? ####
#For now perhaps use the average of the last two complete days for each time point ####

#----EE (minute summation with 10th percentile RMR method)-----#

#For stats I should remove the mice that had problematic peak obesity readings
Daily_EE_zt_stats <- Daily_EE_zttime_alldays %>%
  filter(!ID %in% c(3748, 3749, 3751))

### NEAT --> all complete days ####
#Build multiple linear regression model for NEAT not adjusted for BW or lean #
model_NEAT_zt <- lmer(NEAT_kcal_day ~ SABLE * DRUG + (1 | ID), data = Daily_EE_zt_stats)
summary(model_NEAT_zt)

#Calculate estimated marginal means #
emm_NEAT_zt <- emmeans(model_NEAT_zt, ~ SABLE * DRUG, cov.reduce = mean)
emm_NEAT_zt_df <- as.data.frame(emm_NEAT_zt)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_zt <- contrast(emm_NEAT_zt, method = "pairwise", by = "DRUG")
contrasts_by_group_NEAT_zt_df <- as.data.frame(contrasts_by_group_NEAT_zt)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_zt <- contrast(emm_NEAT_zt, method = "pairwise", by = "SABLE")
contrasts_SABLE_NEAT_zt_df <- as.data.frame(contrasts_by_SABLE_NEAT_zt)

  
