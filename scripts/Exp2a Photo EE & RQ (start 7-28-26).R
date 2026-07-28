
#Experiment 2A analysis of Sable data ####

#Started: 7/28/26 -> used Incretin_Exp2 Sable (Rev. 7-23-26) as starting point
#Revised: 

#Objectives ####
##  1a) Light/Dark period TEE and RMR (separate EE by photoperiod) ####
##  1b) Hourly TEE and RMR ####
##  1c) Day to day variation in TEE and RMR -> is there more variation during BW loss than during peak obesity ####

##  2a) Light/Dark period RQ (separate RQ by photoperiod) ####
##  2b) Hourly RQ ####
##  2c) Is RQ different during BW loss compared to peak obesity? Could reflect that fat and carbs are being handled differently during BW loss ####

## 3) Energy balance -> During BW loss is it positive for vehicle and negative for TZP mice? Is it significantly different?

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
#This method gives integer values for zt_time and accounts for the 6:30 to 18:30 cycle
#It is more explicit in its format/logic compared to original zt_time function
zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

#Get just locomotion for sable_dwn
filter_loc1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
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


# Get just TEE from sable_dwn
filter_TEE1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>%
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
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx?
    by = c("ID", "DateTime", "sable_idx")) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss"))) %>%
  filter(!ID %in% c(3748, 3749, 3751)) #issues with recording

#On some days mouse A was recording in cage I in the morning and mouse B starting recording in cage I
#in the evening. META.csv only indicates dates not hours, so I need to extricate the data from these two mice
#in the code
new_filter_loc_TEE2 <- filter_loc_TEE2 %>%
  group_by(ID, SABLE) %>%
  filter(max(complete_days) <= 2 |complete_days != max(complete_days)) %>% 
  ungroup()
#Now I should have accurate data for days when two mice record in the same cage number


#Calculate EE for complete day 1 ####
filter_loc_TEE3 <- new_filter_loc_TEE2 %>%
  filter(complete_days ==1)

#Try refers to the fact that I added "lights" to the second grouping to allow for calculation of EE in the two photo periods
# RMR calculated using percentile Percentile for RMR --> 
Try_Minute_sum_EE_hr <- filter_loc_TEE3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, DRUG) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, lights) %>% # group by lights to retain photoperiod in df for later use
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
Try_Daily_EE <- Try_Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))

#Calculate sum for light and sum for dark period on each complete day for each mouse
Try_Photo_EE <- Try_Minute_sum_EE_hr %>%
  group_by(ID, lights, SABLE, DRUG) %>%
  summarise(TEE_kcal_photo = sum(TEE_kcal),
            NEAT_kcal_photo = sum(NEAT_kcal),
            RMR_kcal_photo = sum(RMR_kcal),
            diff= abs(TEE_kcal_photo - sum(NEAT_kcal_photo + RMR_kcal_photo)))

#---
#Calculate EE for EACH complete day #### 

Minute_sum_EE_hr_zttime_alldays <- new_filter_loc_TEE2 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, complete_days) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, complete_days, zt_time, lights) %>% #added lights as an additional grouping factor
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

#Calculate sum for light and for dark period on each complete day
Try_Photo_EE_zttime_alldays <- Minute_sum_EE_hr_zttime_alldays %>%
  group_by(ID, lights, SABLE, DRUG, complete_days) %>%
  summarise(TEE_kcal_photo = sum(TEE_kcal),
            NEAT_kcal_photo = sum(NEAT_kcal),
            RMR_kcal_photo = sum(RMR_kcal),
            diff= abs(TEE_kcal_photo - sum(NEAT_kcal_photo + RMR_kcal_photo)))

#---
#---
#---
#---

# RQ ####

filter_RQ1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
  filter(grepl("RQ_*", parameter)) %>%
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

#Calculate hourly average for each mouse
Format_RQ_df <- filter_RQ1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(RQ = value) %>%
  rename(parameter_RQ = parameter) %>%
  rename(fix_value_RQ = fix_value) %>%
  filter(!ID %in% c(3748, 3749, 3751)) #issues with recording of gas -> impacts RQ and EE readings

#On some days mouse A was recording in cage I in the morning and mouse B starting recording in cage I
#in the evening. META.csv only indicates dates not hours, so I need to extricate the data from these two mice
#in the code
New_Format_RQ_df <- Format_RQ_df %>%
  group_by(ID, SABLE) %>%
  filter(max(complete_days) <= 2 |complete_days != max(complete_days)) %>% 
  ungroup()
#Now I should have accurate data for days when two mice record in the same cage number
  
Calc_RQ <- New_Format_RQ_df %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss"))) %>%
  group_by(ID, SABLE, complete_days, hr, DRUG) %>%
  summarise(hr_RQ = mean(RQ)) # this gives average RQ during each hour of each complete day for each mouse at each sable time point
 
#In df = Calc_hr_RQ I have one value for each mouse during each hour during which sable data was recorded.
#there are multiple ways to get a summary statistic for RQ. 
#SABLE and DRUG are the two things that I definitely won't collapse
#Keep in mind that all mice at all time points don't necessarily have 4 complete days (minimum of 2)

#Method 1: Group by SABLE and DRUG ####
#--> no longer distinguishes between individual mice or specific complete days (ID & complete_days collapsed)

Calc_RQ_2 <- Calc_RQ %>%
  group_by(DRUG, SABLE, hr) %>%
  summarise(Avg_hr_RQ = mean(hr_RQ))

##Bar graph ####
Plot_averaged_hr_RQ_exp2a <- ggplot(Calc_RQ_2, aes(x=hr, y=Avg_hr_RQ, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") + #454441
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  facet_wrap(~SABLE) +
  labs(x="Hour",
       y= "RQ",
       title= "Averaged 24hr Respiratory quotient",
       color="Treatment", fill="Treatment", group="Treatment")
Plot_averaged_hr_RQ_exp2a

##LMM: RQ--method 1 ####



#Method 2: Group by SABLE, DRUG, and complete_days ####
#only use complete_days = 1 and 2

Calc_RQ_3 <- Calc_RQ %>%
  filter(complete_days %in% c(1,2)) %>%
  group_by(complete_days, DRUG, SABLE, hr) %>%
  summarise(Avg_hr_RQ = mean(hr_RQ))

# Graph: how can I graph both days on the same graph? Need a new variable which shows hrs as 0-48 rather than 0-24 repeated for each complete_days
#mutate(48_hr= if_else(complete_days =2, hr+24, hr)) 
plot_Calc_RQ_3 <- Calc_RQ_3 %>%
  mutate(hr_Two_days= if_else((complete_days >1), hr+24, hr)) 

##Bar graph ####
Plot_48hr_RQ_exp2a <- ggplot(plot_Calc_RQ_3, aes(x=hr_Two_days , y=Avg_hr_RQ, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") + #454441
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  facet_wrap(~SABLE) +
  labs(x="Hour",
       y= "RQ",
       title= "Respiratory quotient",
       color="Treatment", fill="Treatment", group="Treatment")
Plot_48hr_RQ_exp2a

## LMM: RQ--method 2 ####

