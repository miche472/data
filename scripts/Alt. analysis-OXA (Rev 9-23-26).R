#Data analysis for orexin A ICV injection dose response study in Summer of 2026

#This script is an updated version of "Orexin A analysis (Rev 8-20-26).R"
#Try new analysis approach:
  #Goal 1: Evaluate EE and locomotion relative to baseline (essentially comparing against itself to account for calibration issues)
  #Goal 2: Isolate and evaluate the first two hours post injection even though each mouse will have 
    #a different number of observations due to staggered injection times

#Start: 9-22-26
#Revised: 

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

#Libraries
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

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#Format plot (LM version 3)
format.plot_LM3 <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"), 
  panel.border = element_blank(),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black"),
  plot.title = element_text(size=17, hjust = 0.5, face="bold", vjust=2),
  legend.title=element_text(size=15, face="bold"),
  legend.text=element_text(size=13),
  axis.title.x = element_text(face="bold", size= 15),
  axis.text.x = element_text(size= 13, angle=25, vjust=0.5, hjust=0.7),
  axis.title.y = element_text(face="bold", size= 15),
  axis.text.y = element_text(size = 13))
# Define custom colors
custom_colors_OXA <- c("Baseline" = "darkgray",
  "aCSF"    = "#D9D9D9", 
                        "125pmol" = "#C7E9C0", 
                        "250pmol" = "#A1D99B", 
                        "250pmol" = "#74C476", 
                        "500pmol" = "#41AB5D", 
                        "1000pmol"= "#006D2C")

custom_colors2_OXA <- c("Baseline" = "black","aCSF" = "black", 
                        "125pmol" = "black","250pmol" = "black", 
                        "250pmol" = "black","500pmol" = "black", 
                        "1000pmol"= "black")
#custom_colors_RTIOXA <- c()


#Read in Sable data ####
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
##Filter Sable data to include only mice that are cannulated
sable_dwn_19 <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()
  
##Read in meta data for injection time
read_injection_time <- read_csv("../data/META_INJECTIONS_LG.csv")

##BW for COHORT 19
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3731,3732,3733,3735,3737,3738,3739,3740,3741))%>%
  mutate(PERIOD = case_when(ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-15" ~ "0",
                            ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-17" ~ "1",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-20" ~ "2",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-24" ~ "3",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-27" ~ "4",
                             ID %in% c(3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-29" ~ "5",
                             ID==3731 & DATE =="2026-08-03" ~ "5",
                             ID==3738 & DATE =="2026-07-29" ~ "0",
                             ID==3738 & DATE =="2026-08-03" ~ "1",
                             ID==3738 & DATE =="2026-08-11" ~ "2",
                             ID==3738 & DATE =="2026-08-14" ~ "3",
                             ID==3738 & DATE =="2026-08-17" ~ "4",
                             ID==3738 & DATE =="2026-08-20" ~ "5")) %>%
  drop_na(PERIOD) %>%
  mutate(ID=as.factor(ID),
         PERIOD=as.integer(PERIOD))

#------------------Orexin A (OXA) EE-----------------------------####

##Process injection time data 
injection_time <- read_injection_time %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = row_number()) %>%
  ungroup()

##Identify the exact injection date/time which corresponds to each chunk of recording (i.e. Period of recording) 
#This data frame has one row for each DateTime +ID...it does not yet account for the fact that multiple parameters (10)
#were measured at once for each mouse
 injection_assignment <- sable_dwn_19 %>%
  distinct(ID, DateTime) %>%
  left_join(
    injection_time,
    by = join_by(
      ID,
      DateTime >= INJECTION_DateTime
    ),
    relationship = "many-to-many"
  ) %>%
  group_by(ID, DateTime) %>%
  slice_max(
    INJECTION_DateTime,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
 
##Join injection time and sable data 
 OXA_Sable_joined <- sable_dwn_19 %>%
  left_join(
    injection_assignment,
    by = c("ID", "DateTime"))
 
##Remove INVALID injections
#Remove the first attempt at 3731, 500pmol
Sable_OXA <- OXA_Sable_joined %>%
   filter(!(VALID == "INVALID")) %>% #removes the invalid injection (3731)
   group_by(ID, PERIOD) %>%
   mutate(PERIOD = if_else(ID=="3731" & PERIOD==6, 5, PERIOD)) #assigns the valid injection for 3731 as injection 5 

##Calculate time after injection (in minutes and in hours) 
#Create a column which looks at 20min post injection
Sable_OXA_2 <- Sable_OXA %>%
  ungroup() %>%
  group_by(ID, PERIOD) %>%
  arrange(DateTime) %>%
  mutate(minutes_post_injection = as.numeric(difftime(DateTime,INJECTION_DateTime,units = "mins")),
  first_recording = min(DateTime, na.rm = TRUE), #After first actual recorded observation
  minutes_post_recording = as.numeric( difftime(DateTime, first_recording, units = "mins")))
  #hrs_post_recording = floor(minutes_post_recording / 60) + 1) %>%
  #mutate(Post_injection_plus_20min = INJECTION_DateTime + lubridate::minutes(20))

#Creat seperate dfs with locomotion and with energy expenditure
OXA_EE_1 <-Sable_OXA_2 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

OXA_loc_1 <-Sable_OXA_2 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

OXA_ped_1 <-Sable_OXA_2 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("PedMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

#Join EE and locomotion data
OXA_loc_EE <- OXA_loc_1 %>%
  left_join(
    OXA_EE_1 %>% 
      select(Kcal_Hr, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  left_join(
    OXA_ped_1 %>% 
      select(Ped_meters, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
mutate(DOSE = factor(DOSE,levels = c("aCSF", "125pmol", "250pmol", "500pmol", "1000pmol"))) 


#Add 60 minute bins
#OXA_loc_EE_60min <- OXA_loc_EE %>%
  #ungroup() %>%
  #group_by(ID, DateTime, PERIOD, DOSE) %>% 
  #arrange(DateTime) %>%
  #mutate(TEE_per_min = Kcal_Hr/60) %>%
  #mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  #filter(PERIOD %in% c(1,2,3,4,5))

##Sum of EE (orexin A) ####
#Sum of EE during 120 min post injection 
#OXA_sum_EE_hr_bins <- OXA_loc_EE_60min %>%
  #group_by(ID, recording_bin, PERIOD, DOSE) %>%
  #arrange(DateTime) %>%
  #summarise(sum_EE_60min = sum(TEE_per_min)) %>%
  #ungroup() %>%
    #left_join(
    #BW_COHORT19 %>% 
      #select(ID, BW, PERIOD), 
    #by = c("ID", "PERIOD"))

#-------------------#.
#-------------------#.
#   Baseline EE   ####
#-------------------#.
#-------------------#.
  
Baseline_EE_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

Baseline_loc_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

Baseline_ped_1 <-sable_dwn_19 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("PedMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

#Join EE and locomotion data for baseline
Baseline_loc_EE <- Baseline_loc_1 %>%
  left_join(
    Baseline_EE_1 %>% 
      select(Kcal_Hr, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  left_join(
    Baseline_ped_1 %>% 
      select(Ped_meters, ID, DateTime), 
    by = c("ID", "DateTime")) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(PERIOD = case_when(ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & 
                              date == "2026-07-16"~ "0",
                            ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) &
                              date=="2026-07-17" ~ "0",
                            ID==3738 & date=="2026-08-04"~"0",
                            ID==3738 & date=="2026-08-05"~"0")) %>%
  drop_na(PERIOD) %>%
  mutate(PERIOD= as.integer(PERIOD),
         ID=as.factor(ID)) 

#-----------------------------#.
#-----------------------------#.
#   Combine OXA and baseline  ####
#-----------------------------#.
#-----------------------------#.

Baseline_OXA_loc_EE <- bind_rows(
  Baseline_loc_EE,
  OXA_loc_EE)
#confirmed: obs in Baseline_OXA_loc_EE = obs. Baseline_loc_EE + obs. OXA_loc_EE


#Prepare BW data to be attached ####
BW_COHORT19_oxa_base <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3731,3732,3733,3735,3737,3738,3739,3740,3741))%>%
  mutate(PERIOD = case_when(ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-15" ~ "0",
                            ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-17" ~ "1",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-20" ~ "2",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-24" ~ "3",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-27" ~ "4",
                             ID %in% c(3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-29" ~ "5",
                             ID==3731 & DATE =="2026-08-03" ~ "5",
                             ID==3738 & DATE =="2026-07-29" ~ "0",
                             ID==3738 & DATE =="2026-08-03" ~ "1",
                             ID==3738 & DATE =="2026-08-11" ~ "2",
                             ID==3738 & DATE =="2026-08-14" ~ "3",
                             ID==3738 & DATE =="2026-08-17" ~ "4",
                             ID==3738 & DATE =="2026-08-20" ~ "5")) %>%

  mutate(ID=as.factor(ID),
         PERIOD=as.integer(PERIOD))  %>%
  drop_na(PERIOD)


Combined_loc_EE <- Baseline_OXA_loc_EE %>%
ungroup() %>%
    left_join(
    BW_COHORT19_oxa_base %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) 
#confirmed: obs in Combined_loc_EE = obs. in Baseline_OXA_loc_EE, but this one has an additional variable (BW)

#Calculate sum of EE during the first 2hrs post injection for each mouse after each dose
#When doing stat analysis, include BW as an additional variable
Combined_loc_EE_2 <- Combined_loc_EE %>%
  ungroup() %>%
  group_by(ID, DOSE, DateTime) %>%
  mutate(Kcal = Kcal_Hr/60) %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(flag= if_else(Ped_meters<lag(Ped_meters), "Flag", "No"))
#Need to insert some way to identify when the system was started and stopped
#to allow for ped_meters to restart


Combined_loc_EE_3 <- Combined_loc_EE_2 %>%
  group_by(ID, PERIOD) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    reset = Ped_meters < lag(Ped_meters) & Ped_meters < 1,
    offset = cumsum(replace_na(if_else(reset, lag(Ped_meters), 0), 0)),
    Ped_meters_cumulative = Ped_meters + offset
  ) %>%
  ungroup()
#LEFT OFF on Wed 9/23 ####
#I think this works. Create a check though to make sure that the maximum value for 
#cumulative ped meters occurs at the final observation during each period for each mouse


#---------------------------------------------------------#.
#---------------------------------------------------------#.
#   New analyses    ####
#---------------------------------------------------------#.
#---------------------------------------------------------#.

#In df = OXA_loc_EE there is a variable called minutes_post_injection. Filter for values <120
Two_hrs_OXA_loc_EE <- OXA_loc_EE %>%
  filter(minutes_post_injection <121) %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(total_min = n()) %>%
  summarise(Avg_EE_within_2hr_inject = mean(Kcal_Hr)) %>%
  mutate(GA_group = if_else(ID %in% c(3731,3732,3733,3735), "Low", "High"),
         GA_group = as.factor(GA_group)) %>%
  group_by(ID, DOSE)


#Graph avg EE during first 120 minutes post injection (won't have 120 observations for any of the mice)
plot_EE_2hr_kcal_hr <- ggplot(Two_hrs_OXA_loc_EE, 
       aes(x = factor(DOSE), y = Avg_EE_within_2hr_inject, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  geom_text(data = Two_hrs_OXA_loc_EE %>% group_by(ID) %>%   #label lines with ID
              slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Avg TEE (kcal/hr)", title = "Within 2 hrs post injection") +
  facet_wrap(~GA_group)
plot_EE_2hr_kcal_hr


