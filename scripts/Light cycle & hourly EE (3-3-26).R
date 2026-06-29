#Started: 3-3-26
#Revised: 3-3-26

#Did alternative approach instead of original approach

#Data compilation for SSIB 2026 abstract and APS Summit 2026
#Author: LM
#Goal: Use summation method and per hour average method for entire duration of recording (i.e. not just 24hrs).
#Objective 1: Calculate EE using method 1 and method 2
#Objective 2: Calculate FI using method 1 and method 2
#Objective 3: Calculate EE-FI=energy balance using method 1 and method 2
#Objective 4: Determine if the result is different: 
  #4a: Inspect manually
  #4b: use MLR to look at group comparisons
  #4c: use student t test at each time point. Use ANOVA to compare time points

# Alternative approach ####
#Step 1: For each hour, find the mean of the TEE measurement each minute (units are kcal/hr)
#Step 2:For the 24hrs of measurement, sum the average from each hour --> sum is TEE 
#Step 3 (optional): break this down by light and dark, but would need to add seperately 
#and divide by either 10 or 14hrs to express in terms of each hour (since light 14 hrs / dark 10hrs)
#The consensus paper expresses EE in terms of avg hourly EE during light, dark, and total

#Conclusions:


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

#Objective 1: Calculate EE using method 1 and method 2 ####
  #For both methods, Use carolina's method for complete days

sable_EE <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
    SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                      "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                    "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
      sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  group_by(ID, DRUG,GROUP,DIET_FORMULA,SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,DIET_FORMULA,SABLE) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Method 1: 

#Method 2:



# Alternative approach ####
#Step 1: For each hour, find the mean of the TEE measurement each minute (units are kcal/hr)
#Step 2:For the 24hrs of measurement, sum the average from each hour --> sum is TEE 
#Step 3 (optional): break this down by light and dark, but would need to add seperately 
#and divide by either 10 or 14hrs to express in terms of each hour (since light 14 hrs / dark 10hrs)
#The consensus paper expresses EE in terms of avg hourly EE during light, dark, and total

#Step 1: for each hour, for each ID, find the mean of the EE measurements (in kcal/hr) taken each minute
TEE_avg_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, hr) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  summarise(avg_TEE_hr = mean(Kcal_Hr))

#Step 1b: graph hourly TEE from step 1 (units are kcal/hr)
ggplot(TEE_avg_hr, aes(x = hr, y = avg_TEE_hr, fill = GROUP)) +
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

#Step 2: for the 24hrs of measurement, sum the average EE from each hour for each mouse
TEE_day <- TEE_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP) %>%
  summarise(sum_avg_TEE_hr = sum(avg_TEE_hr)) #gives a daily average TEE per mouse

#Step 2b: Graph daily TEE from step 2
#use df = TEE_day
ggplot(TEE_day, aes(x = SABLE, y = sum_avg_TEE_hr, fill = GROUP)) +
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

#Step 3: break down by light and dark
light_TEE_avg_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, hr, lights) %>%
  arrange(DateTime) %>%
  summarise(avg_TEE_hr = mean(Kcal_Hr))

#sum for light and dark period
light_TEE_day <- light_TEE_avg_hr %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  summarise(sum_avg_TEE_light = sum(avg_TEE_hr))

#Step 3b: graph TEE during lights on and off

#Light cycle on x-axis
ggplot(light_TEE_day, aes(x = lights, y = sum_avg_TEE_light, fill = GROUP)) +
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

#SABLE on x-axis (facet_wrap with light)
ggplot(light_TEE_day, aes(x = SABLE, y = sum_avg_TEE_light, fill = GROUP)) +
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
  facet_wrap(~lights) +
  labs(
    title="TEE during light and dark period (sum of hourly average)",
    x = "Time point",
    y = "TEE (kcal/period)",
    fill = "Treatment group")

#Check if sum of light and dark (sum_avg_TEE_light) is equal to calculated daily TEE (sum_avg_TEE_hr)
#conclusion: these two are equal --> good
df_check <- light_TEE_day  %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(sum_dark_light= sum(sum_avg_TEE_light)) %>%
  left_join(
    TEE_day %>% 
      select(ID, SABLE, GROUP, sum_avg_TEE_hr),
    by = c("ID", "SABLE")) %>%
  mutate(diff= abs(sum_dark_light - sum_avg_TEE_hr))
  






