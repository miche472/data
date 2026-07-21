# Food intake during Sable for mice in Exp 2A (10--> 20nmol/kg TZP)

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

#
filter_FI1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  #filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
  filter(grepl("FoodA_*", parameter)) %>%
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

#Check for issues in the FI measurements
#Filter for negative values for FoodA
negative_foodA <-filter_FI1_revised %>%
  filter(value<0)

#Calculated diff= lag(value)-value --> any values which are negative are issues with recording 
#since I didn't give the mice food during recording
#Calculate

FI_check <- filter_FI1_revised %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  arrange(DateTime) %>%
  mutate(diff=lag(value)-value) %>%
  filter(diff<0)

# cage 1 -> reasonable numbers
# cage 2 -> reasonable numbers
# cage 3 -> no reasonable numbers
# cage 4 -> reasonable numbers
# cage 5 -> reasonable numbers
# cage 6 -> reasonable numbers
# cage 7 -> reasonable numbers
# cage 8 -> reasonable numbers

#Water ####
filter_water1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  #filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
  filter(grepl("Water_*", parameter)) %>%
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

# cage 1 -> all reasonable values
# cage 2 -> not working on any days (super high or super low values...in the thousands)
# cage 3 -> all reasonable values
# cage 4 -> not working on any days (super high or super low values...in the thousands)
# cage 5 -> all reasonable values
# cage 6 -> all reasonable values
# cage 7 -> not working on any days (super high or super low values...in the thousands)
# cage 8 -> not working on any days (super high or super low values...in the thousands)

#BodyMass ####
filter_BW1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  #filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
  filter(grepl("BodyMass_*", parameter)) %>%
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

# cage 1 -> Not working. Every value is zero
# cage 2 -> all reasonable values
# cage 3 -> all reasonable values
# cage 4 -> Not working. Every value is zero
# cage 5 -> all reasonable values
# cage 6 -> all reasonable values
# cage 7 -> all reasonable values
# cage 8 -> all reasonable values
