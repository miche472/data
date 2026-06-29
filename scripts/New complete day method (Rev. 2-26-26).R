# Began with a secion of code from New RMR method & locom, (2-26-26 issues).R
# Started: 2-26-26

#Goals: calculate TEE, RMR, and NEAT by summing one complete day for each ID at each SABLE point
#criteria for selecting a complete day was minimizing the number of missing observations 
#(i.e. minimizing the number of minutes without data)
#compared TEE and RMR calculated using the new method to TEE and RMR calculated using the old method
#(old method was used for my prelim exam)

#Method & reasoning used to select complete days is in script: 
    #New RMR method & locom, (2-26-26 issues).R)

#Additional goals: calculate meters moved and time spent moving during a complete day

#Left off ### 
#To do --> multiple linear regression and bar plot for TEE, RMR, NEAT, locomotion distance, locomotion minutes

#------- Get only the kcal_hr parameter --> value column has kcal_hr data

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
  group_by(ID) #in original version she also grouped by complete_days here
#End of portion taken from original

filter_EE2 <- filter_EE1 %>% #My method to identify complete day 1 and 2 during each sable period
  mutate(LM_complete_day = case_when( 
    #Baseline
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
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725, 3711)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_EE3 <- filter_EE2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706,3707,3713,3719,3720,3721,3722,3727,3728,3729) & LM_complete_day==1 ~ 1,
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
  filter(use_day %in% c(1,2)) #remove observations not from complete day 1 or 2

#---
#In df filter_loc3 use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE3 use mutate to make a column called kcal_hr_ using the value column data

filter_locom3 <- filter_loc3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_energy3 <- filter_EE3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_locom_energy3 <- filter_locom3 %>%
  left_join(
    filter_energy3 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Calculate daily RMR for each mouse ID at each SABLE time point ####
#Filter for times when mouse didn't move. During these moments TEE is entirely resting metabolic rate
# (i.e.) physical activity is not contributing to TEE when mouse isn't moving

#Calculate average daily RMR for each mouse at each Sable time point
ID_RMR3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse moved
  mutate(RMR_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop")

#Calculate average daily RMR within each GROUP at each SABLE time point
GROUP_RMR3 <- filter_locom_energy3 %>%
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

#Calculate TEE using summation of Kcal_HR ####
#Calculate average daily TEE for each mouse at each Sable time point
ID_EE3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  #filter(move==0) %>% #only keep data from minutes when the mouse moved
  mutate(TEE_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(TEE_teske = sum(TEE_per_min), .groups="drop")

#Calculate average daily TEE within each GROUP at each SABLE time point
GROUP_TEE3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  #filter(move==0) %>%
  mutate(TEE_per_min= (Kcal_Hr/60)) %>%
  summarise(TEE_teske = sum(TEE_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_TEE = mean(TEE_teske), 
            SD_GROUP_TEE = sd(TEE_teske),
            SE_GROUP_TEE = sd(TEE_teske)/sqrt(n()))

#Calculate NEAT using summation of Kcal_HR ####
#Calculate average daily TEE for each mouse at each Sable time point
ID_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>% #only keep data from minutes when the mouse moved
  mutate(NEAT_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop")

#Calculate average daily TEE within each GROUP at each SABLE time point
GROUP_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  mutate(NEAT_per_min= (Kcal_Hr/60)) %>%
  summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_NEAT = mean(NEAT_teske), 
            SD_GROUP_NEAT = sd(NEAT_teske),
            SE_GROUP_NEAT = sd(NEAT_teske)/sqrt(n()))

#compare these to raw values using original method
#use df called sable_TEE_adj_BW from NZO_figure7 - TEE corrected by BW (rev. LM)

#Calculate average daily TEE within each GROUP at each SABLE time point
old_GROUP_TEE <- sable_TEE_adj_BW %>%
  ungroup() %>%
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725, 3711)) %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_TEE = mean(tee), 
            SD_GROUP_TEE = sd(tee),
            SE_GROUP_TEE = sd(tee)/sqrt(n()))

#Compare raw RMR calculated using new method and old method
#Get old method of calculating RMR from df called sable_TEE_adj_RMR
old_GROUP_RMR <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725, 3711)) %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_RMR = mean(RMR_kcal_day), 
            SD_GROUP_RMR = sd(RMR_kcal_day),
            SE_GROUP_RMR = sd(RMR_kcal_day)/sqrt(n()))



#---------------------------------------------
#Look at distance traveled during each time period for each group of mice. Also look at this during day and night time
#AllMeters is not zero for the first measurement of a complete day. 
Distance1 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP, LM_complete_day) %>%
  arrange(DateTime) %>%
  mutate(distance_m = All_meters-lag(All_meters),
         distance_m = if_else(distance_m <0, 0, distance_m),
         moving_min = if_else(distance_m >0, 1, 0)) %>%
  drop_na() %>%
  summarise(total_distance = sum(distance_m), 
            total_moving_min = sum(moving_min), .groups="drop")

#---------------------------------------------
#### MANUAL FI for weight cycled MICE at BW loss and BW maintenance ####

#Create df "FI_manual" -> Process FI.csv file (i.e. manual measurements of FI)
FI_manual_cycled <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
  #filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725, 3711)) %>%
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
      #BW loss 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3710) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3714) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
     
      #BW maintenance (SABLE_DAY_18 for all mice --> complete day 2) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
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

Manual_FI_cycled <- FI_manual_cycled %>%
  mutate(ID = factor(ID)) %>%
  filter(GROUP == "Weight cycled") %>%
  filter(SABLE %in% c("BW loss", "BW maintenance")) %>% #includes all 4 Sable days for BW loss and BW maintenance
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  #summarise(avg_corrected_intake_kcal = mean(corrected_intake_kcal)) %>% #In this case, this method is fine, but to be consistent I will use the other method
  mutate(INTAKE_kcal = (INTAKE_GR*3.82)) %>%
  summarise(avg_corrected_intake_kcal = mean(INTAKE_kcal))

#Combine code for EE and locomotion in the top part of this script with FI in sable from NZO_FI_during_Sable (2-23-26)
