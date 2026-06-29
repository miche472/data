#RMR --> calculated using additive method
#Locomotion --> calculated using additive method

#RMR and locomotion calculated in this script use NEAT_LM, 2-23-26 as starting point

#Started:2-24-26
#Revised: 2-24-26

#Left off at bottom of script. Need to figure out why TEE ~ RMR and NEAT > TEE using this method ####

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
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#### Calculate RMR for All SABLE time points for all NZO mice ####
#Use AllMeters to determine the minutes during which a mouse was moving vs staying still

filter_EE1 <-sable_dwn %>%
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
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  #filter(!(ID==3711)) %>% #visually looks like an outlier
  filter(is_complete_day == 1) %>% #keept only days that are complete
  filter(complete_days==1) %>% #keep only the 1st complete day 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Verify that there are 15 mice in each of the Sable periods 
filter_EE1 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID))  #we have 15 NZO in all SABLE periods

#------- Get only the AllMeters parameter --> value column has AllMeters data

#To get complete day 1 and 2, explicitly state zt time and sable day to get compl
#taken from original
filter_loc1 <-sable_dwn %>%
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
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
group_by(ID) #in original version she also grouped by complete_days here
    #End of portion taken from original
    #My method to identify complete day 1 and 2 during each sable period
filter_loc2 <- filter_loc1 %>%
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
  #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
  #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
  #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
  #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  #filter(!(ID == 3711)) %>% #looks visually like an outlier
mutate(
  GROUP = case_when(
    ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
  DRUG = case_when(
    ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
    ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))
  
#Select complete day 1 or 2 for each ID and time point ####
#for filter_loc2 verify number of mice  
filter_loc2 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #we have 15 NZO in all SABLE periods

#Check time points --> look deeper at BW loss issues
filter_loc2 %>% 
  group_by(SABLE, ID, date) %>%
 filter(SABLE=="BW loss")%>% #issue with several mice
  #filter(SABLE=="BW maintenance")%>% #All time points other than BW loss look good from this relatively macro level
  summarise(n_each_ID = n_distinct(hr)) %>%
  print(n=27)

#Determine whether LM_complete_day 1 or 2 has fewer missing observations for each ID at each SABLE time point ####
    #Especially focus on BW loss since I identified issues when using the code above to look at time points
    #Confirm that each mouse has 60 observations during each time point (checking IDs and Sable time points manually)
filter_loc2 %>% 
  group_by(SABLE, ID, hr, LM_complete_day) %>%
  filter(LM_complete_day ==1) %>%
  filter(SABLE=="BW loss") %>%
  filter(ID==3711) %>%
  summarise(n_each_ID = n_distinct(DateTime)) %>%
  print(n=24)

#Decisions about which complete days to use:
#Baseline and Peak obesity: LM_complete_day = 1 --> Only a few missing observations for each mouse
 ###For Baseline and Peak obesity use LM_complete_day = 1 for all mice
#BW loss:
      #Complete day 1:
      #IDs with data for all 24 hrs & <45 missing observations at 8am and 6pm: 3706,3707,3713,3719,3720,3721,3722
      #3708,3714,3716,3726 --> only 12 hrs (hrs 8-19) --> use day 2
      #3727,3728,3729 --> Day 1 is fine, but could be better. However, Day 2 only has 12 hrs for these mice, so use day 1
      #Complete day 2: has most data points for 3708, 3714,3716,3726
  ###For BW loss: Use LM_complete_day = 2 for 3708, 3710, 3714, 3716, 3726. 
  ###For BW loss: Use LM_complete_day = 1 for 3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729 (i.e. all other mice)
#BW maintenance: LM_complete_day = 1 --> about 45 missing observations at hr 18 because of BW and feeding of restricted mice
#BW maintenance: LM_complete_day = 2 --> only about 12 missing observation because restricted mice were fed
 ####For BW maintenance, use LM_complete_day = 2 for all mice
#For BW regain: LM_complete_day = 1 --> about 45min of missing observations due to daily injections, but this is true for all mice
  ###For BW regain, use LM_complete_day = 1 for all mice

#--
#Get only the complete day needed (either day 1 or 2)
filter_loc3 <- filter_loc2 %>%
mutate(use_day = case_when(
  SABLE == "Baseline" & LM_complete_day==1 ~ 1,
  SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
  SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
  SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
  SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
  SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
mutate(
  GROUP = case_when(
    ID %in% c(3706, 3707, 3713, 3716, 3719, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
  DRUG = case_when(
    ID %in% c(3706, 3707, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
    ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain"))) %>%
  filter(use_day %in% c(1,2))

filter_loc3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID))
#---
#In df filter_loc use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE use mutate to make a column called kcal_hr_ using the value column data

filter_locom3 <- filter_loc3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_energy <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_locom_energy <- filter_locom %>%
  left_join(
    filter_energy %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Calculate daily RMR for each mouse ID at each SABLE time point
#Filter for times when mouse didn't move. During these moments TEE is entirely resting metabolic rate
  # (i.e.) physical activity is not contributing to TEE when mouse isn't moving
ID_RMR <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse moved
  mutate(RMR_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop")

#Calculate average daily RMR within each GROUP at each SABLE time point
GROUP_RMR <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>%
  mutate(RMR_per_min= (Kcal_Hr/60)) %>%
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_RMR = mean(RMR_teske), 
            SD_GROUP_RMR = sd(RMR_teske),
            SE_GROUP_RMR = sd(RMR_teske)/sqrt(n()))

#### Extrapolate RMR from EE during the minutes when mice aren't moving ####
# Issue for both methods: RMR seems too high given what TOTAL energy expenditure is (using other methods to calculate TEE)


#Approach B: Find avg hourly RMR for Control and Weight cycled mice during dark and light period
  #Similar method used by Banks et al. 2025 in "Consensus" paper
  #They calculated average hourly TEE during dark and light period
RMR_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(RMR_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, GROUP, lights) %>% n 
  summarise(RMR_GROUP = mean(RMR_kcal_hr), .groups="drop") %>& 

#My original version --> Extrapolate RMR from EE during the minutes when mice aren't moving
  #Calculate avg EE when each mouse isn't moving during dark period and during light period. 
  #Extrapolate to get EE when mouse isn't moving for the entire day by multiplying by the number of hrs
  #of light and dark. Then add together RMR during light and dark period
  RMR_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(RMR_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
  mutate(RMR_per_photo = photo_length_hrs*RMR_kcal_hr) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
#To get total daily RMR for each ID
  summarise(Daily_RMR_ID = sum(RMR_per_photo), .groups="drop") %>%
#To get total daily RMR for each GROUP
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(Daily_RMR_GROUP = mean(Daily_RMR_ID), .groups="drop")
  

#Try using my original method for calculating RMR (above) to calculate TEE
  #Is TEE calculated with this method unreasonably high?
TEE_extrap <- filter_locom_energy %>%
    ungroup() %>%
    arrange(DateTime) %>%     # make sure rows are in time order
    group_by(SABLE, ID, GROUP) %>%
    mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
    #filter(move==0) %>% #only keep data from minutes when the mouse didn't moved
    ungroup() %>%
    group_by(SABLE, ID, GROUP, lights) %>%
    summarise(TEE_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
    group_by(SABLE, ID, GROUP, lights) %>%
    mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
    mutate(TEE_per_photo = photo_length_hrs*TEE_kcal_hr) %>%
    ungroup() %>%
    group_by(ID, SABLE, GROUP) %>%
    #To get total daily RMR for each ID
    summarise(Daily_TEE_ID = sum(TEE_per_photo), .groups="drop") %>%
    #To get total daily RMR for each GROUP
    ungroup() %>%
    group_by(SABLE, GROUP) %>%
    summarise(Daily_TEE_GROUP = mean(Daily_TEE_ID), .groups="drop")
#Daily TEE is essentially the same as daily RMR...this isn't plausible
#Try using this method to calculate NEAT

#Use this method to calculate NEAT
NEAT_extrap <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>% #only keep data from minutes when the mouse moved
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(NEAT_kcal_hr = mean(Kcal_Hr), .groups="drop") %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  mutate(photo_length_hrs = ifelse(lights == "on", 14, 10)) %>%
  mutate(NEAT_per_photo = photo_length_hrs*NEAT_kcal_hr) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  #To get total daily RMR for each ID
  summarise(Daily_NEAT_ID = sum(NEAT_per_photo), .groups="drop") %>%
  #To get total daily RMR for each GROUP
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(Daily_NEAT_GROUP = mean(Daily_NEAT_ID), .groups="drop")
#Values for NEAT are larger than values for TEE...there is an issue with this method of calculation I think
#Need to think more about this

# 2-25-26 ####
#Back to basics: for each minute we get TEE in kcal/hr. For each minute, convert units to kcal/min










