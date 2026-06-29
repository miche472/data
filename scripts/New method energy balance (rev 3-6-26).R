
#Energy balance
#Revised: 3-6-26

#Use Carolina's method for demarcating complete day 1 and 2 (rather than my method for identifying 24hr period)
#Since this method includes gaps in recording for some mice, utilize hourly averages rather than direct summation

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

#THis portion is at the top of the script, but it is not the crux of the script ---
#Determine whether the following two methods of summing TEE for each mouse during each 
#sable during each day gives the same results --> both methods do

#Method 1: Carolina's initial method
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # Method 1 for calculating TEE for each day 
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop")


# Method 2: the method that makes logical sense to me
sable_TEE_data2 <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # Method 2 for calculating TEE for each day 
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee_two = sum(value/60), .groups="drop")

#combine dfs for methods 1 and 2 
sable_TEE_joined <- sable_TEE_data  %>%
  left_join(
    sable_TEE_data2 %>% 
      select(ID, SABLE, tee_two, complete_days),
    by = c("ID", "SABLE", "complete_days")) %>%
  mutate(tee_teetwo= tee - tee_two) %>%
# keep both complete days
filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  left_join(
    ID_TEE3 %>% 
      select(ID, SABLE, , LM_complete_day, TEE_teske),
    by = c("ID", "SABLE")) %>%
  group_by(ID, SABLE) %>%
  mutate(tee_TEE_teske= abs(tee - TEE_teske))

#What if I compare my 24hr method to the average between carolina's day 1 and 2
sable_TEE_data3 <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop")

sable_joined3 <- sable_TEE_data3 %>% 
  left_join(
    ID_TEE3 %>% 
      select(ID, SABLE, TEE_teske),
    by = c("ID", "SABLE")) %>%
  mutate(avgtee_TEE_teske = tee - TEE_teske)

#End of exploratory portion ----

# above is the summation method that I used with my one LM_complete_day compared to Carolina's method. 
#her method was to add the minute by minute energy expenditure over each "complete_days". The issue is that the way she 
#determined complete days has issues, so summing doesn't work well (missing data artificially reduces summed EE)

#How can I get around the issue of inaccurate complete days (i.e. complete days with missing hours or extra hours)?
#What if I take the sum of EE per hour and then calculate the mean hourly energy expenditure for each hour. 
#The issue with averages is that there is a lot of variation between hours, so a lot of detail is lost
#This would be fine for energy balance. However, it is challening for FI because it is difficult to determine how to append the manual food
#...the number of meals to include depends on the time of day when EE was recorded

#For each hour of recording, calculate the average of the minute by minute measurements for EE. The units will be kcal/hr. 
#For all hours of measurement, add together these values. This is a good estimate of the TEE over the day


# Use what is below ####
#---
#Prepare data frame for step 1 (step 1 is section below this one)
sable_TEE_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
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
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))
  
  
#---
#Step 1a: for each hour, for each ID, find the mean of the EE measurements (in kcal/hr) taken each minute
  mod_TEE_avg_hr <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, GROUP) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_TEE_hr = mean(value))
  
#Step 1b: graph hourly TEE from step 1 (units are kcal/hr) --> using "complete days" rather than "LM_complete_day"
  ggplot(mod_TEE_avg_hr, aes(x = hr, y = avg_TEE_hr, fill = GROUP)) +
    geom_bar(stat = "summary", 
             fun = "mean", 
             position = position_dodge(width = 0.9)) +
    geom_errorbar(stat = "summary", 
                  fun.data = mean_se, 
                  position = position_dodge(width = 0.9), 
                  width = 0.3) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 12)) +
    labs(
      title="Average energy expenditure each hour (kcal/hr)",
      x = "Time (Hour)",
      y = "TEE (kcal/hour)",
      fill = "Treatment group") +
    facet_wrap(~SABLE)
  
#---
#Step 2a: for the 24hrs of measurement, sum the average EE from each hour for each mouse
mod_TEE_day <- mod_TEE_avg_hr %>%
    ungroup() %>%
    group_by(SABLE, ID, GROUP) %>%
    summarise(sum_avg_TEE_hr = sum(avg_TEE_hr)) #gives a daily average TEE per mouse
  
#Step 2b: Graph daily TEE from step 2 (raw/unadjusted values)
#use df = mod_TEE_day
ggplot(mod_TEE_day, aes(x = SABLE, y = sum_avg_TEE_hr, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Daily TEE in kcal/day (sum of hourly average)",
    x = "Time point",
    y = "TEE (kcal/day)",
    fill = "Treatment group")

#MLR for TEE
#Build multiple linear regression model for TEE (use df = mod_TEE_day)
model_balance_kcal <- lmer(balance_kcal ~ SABLE * GROUP + (1 | ID), data = E_balance_kcal)
summary(model_balance_kcal)

# Calculate estimated marginal means 
emm_balance_kcal <- emmeans(model_balance_kcal, ~ SABLE * GROUP, cov.reduce = mean)
emm_balance_kcal_df <- as.data.frame(emm_balance_kcal)

# Pairwise contrasts within each GROUP
contrasts_by_group_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "GROUP")
contrasts_by_group_balance_kcal_df <- as.data.frame(contrasts_by_group_balance_kcal)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "SABLE")
contrasts_SABLE_balance_kcal_df <- as.data.frame(contrasts_by_SABLE_balance_kcal)
#---
#Step 3a: Break down by light and dark
mod_light_TEE_avg_hr <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, hr, lights) %>%
  arrange(DateTime) %>%
  summarise(avg_TEE_hr = mean(value))

#sum for light and dark period
mod_light_TEE_day <- mod_light_TEE_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(sum_avg_TEE_light = sum(avg_TEE_hr))

#Step 3b: graph TEE during lights on and off
#Light cycle on x-axis
ggplot(mod_light_TEE_day, aes(x = lights, y = sum_avg_TEE_light, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~SABLE) +
  labs(
    title="TEE during light and dark period (sum of hourly average)",
    x = "Light cycle",
    y = "TEE (kcal/period)",
    fill = "Treatment group")

#Step 3c: Break down by light and dark, but use kcal/hr as the metric for comparing light vs dark
#this is necessary because the duration of light and dark are not equal
#This is from above --> get 
mod_light_period_hrly_avg <- sable_TEE_modCS %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  arrange(DateTime) %>%
  summarise(avg_hrly_TEE = mean(value))

#Graph 3c
ggplot(mod_light_period_hrly_avg, aes(x = lights, y = avg_hrly_TEE, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~SABLE) +
  labs(
    title="Average hourly EE during light and dark period",
    x = "Light cycle",
    y = "Hourly EE (kcal/hr)",
    fill = "Treatment group")

#FI ####
sable_FI_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("FoodA_*", parameter)) %>% 
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
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Find average intake during hr 0, 1, 2, 3, 4, etc. for each ID at each SABLE time point (next step is to 
#Add these values to get 1 day of intake for each mouse at each time point
#Remove Sable FI measurements for weight cycled mice during BW loss and BW maintenance
mod_sable_FI <- sable_FI_modCS %>%
  group_by(ID, SABLE, GROUP, complete_days) %>% #should I group by complete days here?
  arrange(DateTime) %>%       # make sure data is ordered
  mutate(
    intake = lag(value) - value,       # change in food from one minute to the next
    intake = if_else(intake < 0, 0, intake)) %>%  #If mouse doesn't eat between min. x and min. x+1, intake=0
  mutate (intake = if_else(intake >2, 0, intake)) %>% #removes values that are illogically high for 1 minute
  drop_na() %>% 
  select(-KCAL_PER_GR) %>% # this is wrong, so remove and replace with KCAL_G
  ungroup () %>%
  group_by(ID, SABLE, hr, GROUP, complete_days) %>%
  summarise(hr_eaten_gr = sum(intake)) %>%
  ungroup () %>% 
group_by(ID, SABLE, hr, GROUP) %>% #get average intake for each hour across complete days 1 and 2
  summarise(avg_hr_eaten_gr = mean(hr_eaten_gr)) %>%
  ungroup () %>% 
group_by(ID, SABLE, GROUP) %>% #calculate sum of average intake during each of the 24hrs
  summarise(Day_FI_gr = sum(avg_hr_eaten_gr)) %>%
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
  mutate(Day_FI_kcal = Day_FI_gr*KCAL_G) %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice

#Manual FI (first create a df for all mice, then process this to filter for only weight cycled mice at BW loss/maintenance)

#Create df "FI_manual" -> Process FI.csv file (i.e. manual measurements of FI)
FI_manual <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
  #filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
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
      #Baseline (sable days 1-7)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
      sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance", "BW regain"))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  ungroup() 
  
FI_manual_cycled <- FI_manual %>%
  mutate(ID = factor(ID)) %>%
  filter(GROUP == "Weight cycled") %>%
  filter(SABLE %in% c("BW loss", "BW maintenance")) %>% 
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  #summarise(avg_corrected_intake_kcal = mean(corrected_intake_kcal)) %>% #In this case, this method is fine, but to be consistent I will use the other method
  mutate(INTAKE_kcal = (INTAKE_GR*3.82)) %>%
  summarise(Day_FI_kcal = mean(INTAKE_kcal)) 
  
#Join Sable measured FI and manually measured FI
Join_mod_sable_FI <- mod_sable_FI %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  select(SABLE, GROUP, ID, Day_FI_kcal)

FI_for_balance <- bind_rows(FI_manual_cycled, Join_mod_sable_FI)  %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP)

#---Energy balance ####
#EE is calculated in this script. For EE, use df= mod_TEE_day where daily EE is sum_avg_TEE_hr
#FI is calculated in this script. For FI, use df = FI_for_balance where daily FI is Day_FI_kcal

#Join the FI and the EE data frames; create new variable called E_balance_kcal
E_balance_kcal <- mod_TEE_day  %>%
  left_join(
    FI_for_balance %>% 
      select(ID, SABLE, GROUP, Day_FI_kcal),
    by = c("ID", "SABLE", "GROUP")) %>%
  mutate(balance_kcal= Day_FI_kcal - sum_avg_TEE_hr)

#Graph Energy balance (raw values) 
ggplot(E_balance_kcal, aes(x = SABLE, y = balance_kcal, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Daily energy balance (FI - TEE)",
    x = "Time point",
    y = "Energy balance (kcal)",
    fill = "Treatment group")

#Conduct MLR on energy balance
  #Multiple linear regression to determine if daily energy balance is significantly different
  #between time points and control/weight cycled mice

#Build multiple linear regression model for energy balance (use df = E_balance_kcal )
model_balance_kcal <- lmer(balance_kcal ~ SABLE * GROUP + (1 | ID), data = E_balance_kcal)
summary(model_balance_kcal)

# Calculate estimated marginal means 
emm_balance_kcal <- emmeans(model_balance_kcal, ~ SABLE * GROUP, cov.reduce = mean)
emm_balance_kcal_df <- as.data.frame(emm_balance_kcal)

# Pairwise contrasts within each GROUP
contrasts_by_group_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "GROUP")
contrasts_by_group_balance_kcal_df <- as.data.frame(contrasts_by_group_balance_kcal)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_balance_kcal <- contrast(emm_balance_kcal, method = "pairwise", by = "SABLE")
contrasts_SABLE_balance_kcal_df <- as.data.frame(contrasts_by_SABLE_balance_kcal)

#Graph estimated marginal means
barplot_emm_emm_balance_kcal <- emm_balance_kcal_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
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
    title = "Predicted average daily energy balance (kcal)",
    y = "Energy balance (FI-EE, kcal)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_emm_balance_kcal


#start delete ####
#FI during sable is calculated in script called " ". For FI, use 
mod_sable_FI2 <- mod_sable_FI %>%
  group_by(ID, complete_days, GROUP, SABLE) %>%
  arrange(DateTime) %>%       # make sure data is ordered
  mutate(
    intake = lag(value) - value,       # change in food from one minute to the next
    intake = if_else(intake < 0, 0, intake)) %>%  #If mouse doesn't eat between min. x and min. x+1, intake=0
  #filter(intake > 1) #delete this line after using it to see if there are any times when there is a huge jump between minute to minute measurements
  mutate (intake = if_else(intake >8, 0, intake)) %>% #removes values that are illogically high for 1 minute
  drop_na() %>% 
  select(-KCAL_PER_GR) %>% # this is wrong, so remove and replace with KCAL_G
  ungroup () %>% 
  group_by(ID, SABLE, hr, GROUP, complete_days) %>%
  summarise(total_eaten_gr = sum(intake), .groups = "drop") %>%
  ungroup () %>% 
  group_by(ID, SABLE, hr, GROUP) %>%
  summarise(avg_total_eaten_gr = mean(total_eaten_gr), .groups = "drop") %>%
  ungroup () %>% 
  group_by(ID, SABLE, GROUP) %>%
  summarise(Day_FI_gr = sum(avg_total_eaten_gr), .groups = "drop")

summarise(total_eaten_gr = sum(intake), .groups = "drop") %>%     # total FI per day in grams
  
  summarise(total_eaten_gr = sum(intake), .groups = "drop") %>% #average hourly intake
  summarise(total_eaten_gr = mean(total_eaten_gr), .groups = "drop") %>%
  
  mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
  mutate(total_eaten_kcal = total_eaten_gr*KCAL_G) %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>% #UPDATE: need to get average across the two complete_days...how to do this?
  summarise (avg_corrected_intake_kcal = mean(total_eaten_kcal)) %>%
  filter(!(GROUP== "Weight cycled" & SABLE %in% c("BW loss", "BW maintenance"))) #Remove measured FI for restricted mice

#End delete

#Add change in BW over the course of Sable --> quality control

  #consensus paper says to do this when looking at energy balance for quality control
  # if mice have positive energy balance during recording they should gain weight.
  #If this isn't true, there is likely an issue with the recording for a cage and
  #that cage should be removed from analysis
  #consensus method advises to use first and last BW measurement of recording for each mouse

#Prepare data frame for step 1 (step 1 is section below this one)
sable_BW_modCS <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("BodyMass_*", parameter)) %>% 
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
  filter(!ID %in% c(3715,3712)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  filter(is_complete_day == 1) %>%
  filter(complete_days %in% c(1,2)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Check if there are any BWs that are zero or less than zero

#Find first and last BW measurement for each mouse during each Sable time point


Find_sable_BW_modCS <- sable_BW_modCS %>%
  filter(value > 0) %>%
  filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  mutate(
    start_bw = first(value),
    first_DateTime = first(DateTime),
    last_DateTime = last(DateTime))




Find_sable_BW_modCS <- sable_BW_modCS %>%
  filter(value > 0) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  mutate(
    first_bw = first(value),
    last_bw  = last(value),
    bw_change = last_bw - first_bw)

Find_sable_BW_modCS <- sable_BW_modCS %>%
  filter(value > 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime, .by_group = TRUE) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop"
  )

Find_sable_BW_modCS2 <- sable_BW_modCS %>%
  filter(value > 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop"
  )

Find_sable_BW_modCS2 <- sable_BW_modCS %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  #filter(value == 0) %>%
  #filter(ID == 3706) %>%
  group_by(ID, SABLE) %>% 
  arrange(DateTime) %>% 
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop") %>%
  filter(bw_change == 0)
