# Energy expenditure in GLP-1 experiment 2

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

#functions####
#zt_time <- function(hr){
  #return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))}

#For now say the light cycle is 7pm to 6am...need to extract actual time (not just hr) from DateTime to be able to do 6:30am and pm
#functions####
zt_time <- function(hr){
  return(if_else(hr >= 19 & hr <= 23, hr-19, hr+5))}

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#The light cycle is now 6:30am - 6:30pm. When the code was initially written it was 6am to 8pm

filter_loc1 <-sable_dwn %>%
  filter(COHORT >18) %>%   
  mutate(lights  = if_else(hr %in% c(19, 20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4") ~ "Peak obesity")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) %>%
mutate(LM_complete_day = case_when(
  #Peak obesity (Sable wasn't stopped so no missing data)
  sable_idx=="SABLE_DAY_2" ~ 1,
  sable_idx=="SABLE_DAY_3" ~ 2,
  #BW loss (since sable will be stopped for injections so some data will be lose)
  zt_time %in% c(0, 1, 2, 3, 4) & sable_idx=="SABLE_DAY_5" ~ 1,
  zt_time>4 & zt_time<24 & sable_idx=="SABLE_DAY_6" ~ 1,
  zt_time %in% c(0, 1, 2, 3, 4) & sable_idx=="SABLE_DAY_6" ~ 2,
  zt_time>4 & zt_time<24 & sable_idx=="SABLE_DAY_7" ~ 2))

#Figure out if there is missing data
filter_loc1 %>% 
  group_by(SABLE, ID, date) %>%
  filter(SABLE=="Peak obesity")%>% #issue with several mice
  #filter(SABLE=="BW maintenance")%>% #All time points other than BW loss look good from this relatively macro level
  summarise(n_each_ID = n_distinct(hr)) %>%
  #print(n=32) %>%
  filter(n_each_ID== 24)
#Good, all 8 mice have measurements during each hour for a 48hr period

filter_loc1 %>% 
  group_by(SABLE, ID, hr, LM_complete_day) %>%
  filter(LM_complete_day ==1) %>%
  filter(SABLE=="Peak obesity") %>%
  filter(ID==3743) %>%
  summarise(n_each_ID = n_distinct(DateTime)) %>%
  print(n=24)

# Get just TEE 
filter_TEE1 <-sable_dwn %>%
  filter(COHORT >18) %>%   
  mutate(lights  = if_else(hr %in% c(19, 20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4") ~ "Peak obesity")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) %>%
  mutate(LM_complete_day = case_when(
    #Peak obesity (Sable wasn't stopped so no missing data)
    sable_idx=="SABLE_DAY_2" ~ 1,
    sable_idx=="SABLE_DAY_3" ~ 2,
    #BW loss (since sable will be stopped for injections so some data will be lose)
    zt_time %in% c(0, 1, 2, 3, 4) & sable_idx=="SABLE_DAY_5" ~ 1,
    zt_time>4 & zt_time<24 & sable_idx=="SABLE_DAY_6" ~ 1,
    zt_time %in% c(0, 1, 2, 3, 4) & sable_idx=="SABLE_DAY_6" ~ 2,
    zt_time>4 & zt_time<24 & sable_idx=="SABLE_DAY_7" ~ 2))

##---
#In df filter_loc2 use mutate to make a column called AllMeters_ using the value column data
#In df filter_TEE2 use mutate to make a column called kcal_hr_ using the value column data

filter_loc2 <- filter_loc1 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_TEE2 <- filter_TEE1 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_loc_TEE2 <- filter_loc2 %>%
  left_join(
    filter_TEE2 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx")) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss")))

#OPTIONS ####
#I have almost 4 days of recording. For now filter so I am just using the first complete day
#LM_complete_day = 1
filter_loc_TEE3 <- filter_loc_TEE2 %>%
  ungroup() %>%
  group_by(ID) %>%
  filter(LM_complete_day ==1) %>%
  ungroup()

#Taken from 4/16 script and modified

# In previoius scripts I used whether or not the mouse moved between two minutes to classify
#the EE during that minute as entirely NEAT or entirely RMR. This results in an 
#overestimate of NEATsince even when I mouse is moving some of its TEE is RMR.
#A way to get around this issue is to calculate the 10th percentile of TEE for a given 
#observation period (hr, lights on/off, or daily [aka global])--> for now, do by hr
#For min. when mouse didn't move: RMR=TEE and NEAT=0. 
#For min. when mouse moved, NEAT = TEE - (RMR=10th percentile of RMR values across that hr for the ID)

#---------------------------------------------------------------------#
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
  mutate(lights= if_else(hr %in% c(19,20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% #add photo-period
  ungroup() %>%
  group_by(ID, SABLE, DRUG, lights) %>%
  summarise(TEE_kcal_photo = sum(TEE_kcal), 
            NEAT_kcal_photo = sum(NEAT_kcal),
            RMR_kcal_photo = sum(RMR_kcal))

  
