
#Data analysis for orexin A ICV injection dose response study in Summer of 2026

#This script is a simplified OXA dose response (Rev 8-18-26)
  #(this script is an improved, more parsimonious version)


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
  minutes_post_recording = as.numeric( difftime(DateTime, first_recording, units = "mins")),
  hrs_post_recording = floor(minutes_post_recording / 60) + 1) %>%
  
  mutate(Post_injection_plus_20min = INJECTION_DateTime + lubridate::minutes(20))

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
OXA_loc_EE_60min <- OXA_loc_EE %>%
  ungroup() %>%
  group_by(ID, DateTime, PERIOD, DOSE) %>% 
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  filter(PERIOD %in% c(1,2,3,4,5))

##Sum of EE (orexin A) ####
#Sum of EE during each 60min bin (starting when recording started) 
OXA_sum_EE_hr_bins <- OXA_loc_EE_60min %>%
  group_by(ID, recording_bin, PERIOD, DOSE) %>%
  arrange(DateTime) %>%
  summarise(sum_EE_60min = sum(TEE_per_min)) %>%
  ungroup() %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))


##Avg EE (orexin A) ####
###Avg by 60min bin
OXA_avg_EE_hr_bins <- OXA_loc_EE_60min %>%
  group_by(ID, PERIOD, DOSE, recording_bin) %>%
  arrange(DateTime) %>%
  summarise(avg_EE_60min = mean(Kcal_Hr)) %>%
  ungroup() %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) 

#---
#---
#---

#Baseline EE ####
  
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

#Use 5pm as the theoretical "start of recording" for baseline, to try to  
#roughly match circadian pattern of injection days
Baseline_loc_EE_60min <- Baseline_loc_EE %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(start_recording = case_when(
      ID == 3738 ~ lubridate::ymd_hms("2026-08-04 17:00:00"),
      TRUE ~ lubridate::ymd_hms("2026-07-16 17:00:00")),
    minutes_post_recording = as.numeric(difftime(DateTime, start_recording, units = "mins")),
    minutes_post_recording = if_else(minutes_post_recording < 0,NA_real_,minutes_post_recording)) %>% #makes minutes before the faux injection be NA rather than a negative number
  ungroup() %>%
  group_by(ID, DateTime) %>% 
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  drop_na(minutes_post_recording)

##Prepare BW data to be attached 
BW_COHORT19_Baseline <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3731,3732,3733,3735,3737,3738,3739,3740,3741))%>%
  mutate(PERIOD = case_when(ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-15" ~ "0",
                             ID==3738 & DATE =="2026-07-29" ~ "0"),
         PERIOD=as.integer(PERIOD),
         ID=as.factor(ID)) %>%
  drop_na(PERIOD)

##Sum of EE (baseline) ####
###Calculate sum of EE during each 60min bin (starts at 17:00:00 on 2026-07-16) 
Baseline_sum_EE_hr_bins <- Baseline_loc_EE_60min %>%
  group_by(ID, recording_bin, PERIOD) %>%
  arrange(DateTime) %>%
  summarise(sum_EE_60min = sum(TEE_per_min)) %>%
  ungroup() %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) %>%
  mutate(DOSE= case_when(PERIOD =="0"~"Baseline"),
         DOSE = as.factor(DOSE))


## Avg EE (baseline) ####
###Avg by 60min bin
Baseline_avg_EE_hr_bins <- Baseline_loc_EE_60min %>%
  group_by(ID, recording_bin, PERIOD) %>%
  arrange(DateTime) %>%
  summarise(avg_EE_60min = mean(Kcal_Hr)) %>%
  ungroup() %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))  %>%
  mutate(DOSE= case_when(PERIOD =="0"~"Baseline"),
         DOSE = as.factor(DOSE))

#---
#---
#---
#---------------------Merge orexin A & Baseline EE---------------------####

##---Sum of EE ####
Combined_sum_EE_hr_bins <- bind_rows(
  Baseline_sum_EE_hr_bins,
  OXA_sum_EE_hr_bins)

Summary_Combined_sum_EE_hr_bins <- Combined_sum_EE_hr_bins %>%
  group_by(ID, PERIOD, DOSE) %>%
  summarise(
    sum_EE_0_2hr = sum(sum_EE_60min[recording_bin >= 0 & recording_bin < 2], na.rm = TRUE),
    sum_EE_2_4hr = sum(sum_EE_60min[recording_bin >= 2 & recording_bin < 4], na.rm = TRUE),
    sum_EE_0_4hr = sum(sum_EE_60min[recording_bin >= 0 & recording_bin < 4], na.rm = TRUE),
    sum_EE_0_24hr = sum(sum_EE_60min[recording_bin >= 0 & recording_bin < 24], na.rm = TRUE),
    .groups = "drop") %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) %>%
  ungroup()

## -- Avg EE (kcal/hr)  ####
Combined_avg_EE_hr_bins <- bind_rows(
  Baseline_avg_EE_hr_bins,
  OXA_avg_EE_hr_bins)

Summary_Combined_avg_EE_hr_bins <- Combined_avg_EE_hr_bins %>%
  group_by(ID, PERIOD, DOSE) %>%
  summarise(
    Avg_EE_0_2hr = mean(avg_EE_60min[recording_bin >= 0 & recording_bin < 2], na.rm = TRUE),
    Avg_EE_2_4hr = mean(avg_EE_60min[recording_bin >= 2 & recording_bin < 4], na.rm = TRUE),
    Avg_EE_0_4hr = mean(avg_EE_60min[recording_bin >= 0 & recording_bin < 4], na.rm = TRUE),
    Avg_EE_0_24hr = mean(avg_EE_60min[recording_bin >= 0 & recording_bin < 24], na.rm = TRUE),
    .groups = "drop") %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))

# GRAPHS (average hourly EE) ####

## Hours 0-2 post start of recording (kcal/hr) ####
###Avg EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_2hr_avg <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_2hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "0 - 2 hrs post")
plot_EE_0_2hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_2hr_avg,
       filename="EE_0_2hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_2hr_median <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_2hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "0 - 2 hrs post")
plot_EE_0_2hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_2hr_median,
       filename="EE_0_2hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


## Hours 2-4 post start of recording (average) ####
###Avg EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_2_4hr_avg <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_2_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)",title = "2 - 4 hrs post")
plot_EE_2_4hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_2_4hr_avg,
       filename="EE_2_4hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_2_4hr_median <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_2_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "2 - 4 hrs post")
plot_EE_2_4hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_2_4hr_median,
       filename="EE_2_4hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


## Hours 0-4 post start of recording (average) ####
###Avg EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_4hr_avg <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)",title = "0 - 4 hrs post")
plot_EE_0_4hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_4hr_avg,
       filename="EE_0_4hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_4hr_median <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "0 - 4 hrs post")
plot_EE_0_4hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_4hr_median,
       filename="EE_0_4hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")



## Hours 0-24 post start of recording (kcal/hr) ####
###Avg EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_24hr_avg <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_24hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "0 - 24 hrs post")
plot_EE_0_24hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_24hr_avg,
       filename="EE_0_24hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_24hr_median <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(DOSE), y = Avg_EE_0_2hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "0 - 24 hrs post")
plot_EE_0_24hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_EE_0_24hr_median,
       filename="EE_0_24hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


#------------------------------------------------------# -
#------------------------------------------------------# -
#---------------------Locomotion-----------------------#
#------------------------------------------------------# -
#------------------------------------------------------# -
#Orexin A (OXA) Locomotion ####
#In the EE section I created a dataframe with AllMeters, Ped_Meters, & Kcal_Hr
#Use this df to analyze locomotion (df = OXA_loc_EE) 

##All_Meters ####
#Add 60 minute bins and identify where the counter resets 
        #(i.e. when recording is stopped and restarted)
OXA_loc_60min <- OXA_loc_EE %>%
  filter(PERIOD %in% c(1, 2, 3, 4, 5)) %>%
  group_by(ID, PERIOD, DOSE) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60),
    Move_m = All_meters - lag(All_meters),
    # A negative difference means the counter reset
    Move_m = if_else(Move_m < 0, NA_real_, Move_m)) %>%
  drop_na(Move_m) %>%
  ungroup()

#Calculate total meters traveled (All_Meters) within each 60min bin
OXA_sum_loc_hr_bins <- OXA_loc_60min %>%
  group_by(ID, recording_bin, PERIOD, DOSE) %>%
  summarise(Sum_move_m = sum(Move_m, na.rm = TRUE),.groups = "drop") %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))

##Ped_Meters ####
#Add 60 minute bins and identify where the counter resets 
        #(i.e. when recording is stopped and restarted)
OXA_ped_60min <- OXA_loc_EE %>%
  filter(PERIOD %in% c(1, 2, 3, 4, 5)) %>%
  group_by(ID, PERIOD, DOSE) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60),
    Ped_m = Ped_meters - lag(Ped_meters),
    Ped_m = if_else(Ped_m < 0, NA_real_, Ped_m)) %>%  # A negative difference means the counter reset
  drop_na(Ped_m) %>%
  ungroup()

#Calculate total meters traveled (All_Meters) within each 60min bin
OXA_sum_ped_hr_bins <- OXA_ped_60min %>%
  group_by(ID, recording_bin, PERIOD, DOSE) %>%
  summarise(Sum_ped_m = sum(Ped_m, na.rm = TRUE),.groups = "drop") %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))


#Baseline Locomotion ####
#In the Baseline EE section I created a dataframe with All_meters, Ped_meters, & Kcal_Hr
#Use this df to analyze locomotion (df = Baseline_loc_EE_60min) 
#Note that this df already has the 60min window created

##All_meters ####
#Add 60 minute bins and identify where the counter resets 
        #(i.e. when recording is stopped and restarted)
Baseline_loc_60min <- Baseline_loc_EE %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(start_recording = case_when(
      ID == 3738 ~ lubridate::ymd_hms("2026-08-04 17:00:00"),
      TRUE ~ lubridate::ymd_hms("2026-07-16 17:00:00")),
    minutes_post_recording = as.numeric(difftime(DateTime, start_recording, units = "mins")),
    minutes_post_recording = if_else(minutes_post_recording < 0,NA_real_,minutes_post_recording)) %>% #makes minutes before the faux injection be NA rather than a negative number
  drop_na(minutes_post_recording) %>%
  ungroup() %>%
  group_by(ID, PERIOD) %>% 
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60),
    Move_m = All_meters - lag(All_meters),
    Move_m = if_else(Move_m < 0, NA_real_, Move_m)) %>% # A negative difference means the counter reset
  drop_na(Move_m) %>%
  ungroup()  %>%
  mutate(DOSE = case_when(PERIOD =="0"~"Baseline")) %>%
  mutate(DOSE = as.factor (DOSE))

#Calculate total meters traveled (All_Meters) within each 60min bin
Baseline_sum_loc_hr_bins <- Baseline_loc_60min %>%
  group_by(ID, recording_bin, PERIOD, DOSE) %>%
  summarise(Sum_move_m = sum(Move_m, na.rm = TRUE),.groups = "drop") %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))

##Ped_meters ####
#Add 60 minute bins and identify where the counter resets 
        #(i.e. when recording is stopped and restarted)
Baseline_ped_60min <- Baseline_loc_EE %>%
  group_by(ID) %>%
  arrange(DateTime) %>%
  mutate(start_recording = case_when(
      ID == 3738 ~ lubridate::ymd_hms("2026-08-04 17:00:00"),
      TRUE ~ lubridate::ymd_hms("2026-07-16 17:00:00")),
    minutes_post_recording = as.numeric(difftime(DateTime, start_recording, units = "mins")),
    minutes_post_recording = if_else(minutes_post_recording < 0,NA_real_,minutes_post_recording)) %>% #makes minutes before the faux injection be NA rather than a negative number
  drop_na(minutes_post_recording) %>%
  ungroup() %>%
  group_by(ID, PERIOD) %>% 
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60),
    Ped_m = Ped_meters - lag(Ped_meters),
    Ped_m = if_else(Ped_m < 0, NA_real_, Ped_m)) %>% # A negative difference means the counter reset
  drop_na(Ped_m) %>%
  ungroup()  %>%
  mutate(DOSE = case_when(PERIOD =="0"~"Baseline")) %>%
  mutate(DOSE = as.factor (DOSE))

#Calculate total meters traveled (Ped_meters) within each 60min bin
Baseline_sum_ped_hr_bins <- Baseline_ped_60min %>%
  group_by(ID, recording_bin, PERIOD, DOSE) %>%
  summarise(Sum_ped_m = sum(Ped_m, na.rm = TRUE),.groups = "drop") %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD"))


#---------------------Merge orexin A & Baseline locomotion---------------------####

##---Sum of all movement (All_meters)) ####
Combined_sum_loc_hr_bins <- bind_rows(
  Baseline_sum_loc_hr_bins,
  OXA_sum_loc_hr_bins)

Summary_Combined_sum_loc_hr_bins <- Combined_sum_loc_hr_bins %>%
  group_by(ID, PERIOD, DOSE) %>%
  summarise(
    sum_loc_0_2hr = sum(Sum_move_m[recording_bin >= 0 & recording_bin < 2], na.rm = TRUE),
    sum_loc_2_4hr = sum(Sum_move_m[recording_bin >= 2 & recording_bin < 4], na.rm = TRUE),
    sum_loc_0_4hr = sum(Sum_move_m[recording_bin >= 0 & recording_bin < 4], na.rm = TRUE),
    sum_loc_0_24hr = sum(Sum_move_m[recording_bin >= 0 & recording_bin < 24], na.rm = TRUE),
    .groups = "drop") %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) %>%
  ungroup() %>%
    #3740 is outlier at 250pmol
  ungroup() %>%
  group_by(ID) %>%
  filter(!(ID==3740))

##---Sum of walking (Ped_meters)) ####
Combined_sum_ped_hr_bins <- bind_rows(
  Baseline_sum_ped_hr_bins,
  OXA_sum_ped_hr_bins)

Summary_Combined_sum_ped_hr_bins <- Combined_sum_ped_hr_bins %>%
  group_by(ID, PERIOD, DOSE) %>%
  summarise(
    sum_ped_0_2hr = sum(Sum_ped_m[recording_bin >= 0 & recording_bin < 2], na.rm = TRUE),
    sum_ped_2_4hr = sum(Sum_ped_m[recording_bin >= 2 & recording_bin < 4], na.rm = TRUE),
    sum_ped_0_4hr = sum(Sum_ped_m[recording_bin >= 0 & recording_bin < 4], na.rm = TRUE),
    sum_ped_0_24hr = sum(Sum_ped_m[recording_bin >= 0 & recording_bin < 24], na.rm = TRUE),
    .groups = "drop") %>%
    left_join(
    BW_COHORT19 %>% 
      select(ID, BW, PERIOD), 
    by = c("ID", "PERIOD")) %>%
  #3740 is outlier at 250pmol
  ungroup() %>%
  group_by(ID) %>%
  filter(!(ID==3740))


# GRAPHS (total movement or walking, meters) ####

## Hours 0-2 post start of recording (total meters) ####
### Total Distance (m) vs Dose (including baseline) ####
plot_loc_0_2hr_avg <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_2hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 2 hrs post")
plot_loc_0_2hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_2hr_avg,
       filename="loc_0_2hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median Total distance (m) vs Dose (including baseline) ####
plot_loc_0_2hr_median <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_2hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 2 hrs post")
plot_loc_0_2hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_2hr_median,
       filename="loc_0_2hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


### Walking Distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_0_2hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "0 - 2 hrs post")

### Median walking distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_0_2hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "0 - 2 hrs post")


## Hours 2-4 post start of recording (total meters) ####
### Total Distance (m) vs Dose (including baseline) ####
plot_loc_2_4hr_avg <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_2_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "2 - 4 hrs post")
plot_loc_2_4hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_2_4hr_avg,
       filename="loc_2_4hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median Total distance (m) vs Dose (including baseline) ####
plot_loc_2_4hr_median <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_2_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "2 - 4 hrs post")
plot_loc_2_4hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_2_4hr_median,
       filename="loc_2_4hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


### Walking Distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_2_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "2 - 4 hrs post")

### Median walking distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_2_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "2 - 4 hrs post")


## Hours 0-4 post start of recording (total meters) ####
### Total Distance (m) vs Dose (including baseline) ####
plot_loc_0_4hr_avg <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 4 hrs post")
plot_loc_0_4hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_4hr_avg,
       filename="loc_0_4hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median Total distance (m) vs Dose (including baseline) ####
plot_loc_0_4hr_median <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 4 hrs post")
plot_loc_0_4hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_4hr_median,
       filename="loc_0_4hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Walking Distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_0_4hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "0 - 4 hrs post")

### Median walking distance (m) vs Dose (including baseline) ####
ggplot(Summary_Combined_sum_ped_hr_bins, 
       aes(x = factor(DOSE), y = sum_ped_0_4hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_ped_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "Ambulation (m)", title = "0 - 4 hrs post")

## Hours 0-24 post start of recording (total meters) ####
### Total Distance (m) vs Dose (including baseline) ####
plot_loc_0_24hr_avg <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_24hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 24 hrs post")
plot_loc_0_24hr_avg

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_24hr_avg,
       filename="loc_0_24hr_avg.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")

### Median Total distance (m) vs Dose (including baseline) ####
plot_loc_0_24hr_median <- ggplot(Summary_Combined_sum_loc_hr_bins, 
       aes(x = factor(DOSE), y = sum_loc_0_24hr, fill = factor(DOSE))) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, whisker.linewidth = 0.7, staplewidth = 0.5) +
  geom_line(aes(group = ID), color = "gray50", linewidth = 0.7, alpha = 0.6) + # Lines connecting the same mouse across doses
  geom_jitter(aes(color = factor(DOSE)), width = 0.12, size = 2, alpha = 0.7) +   # Individual observations
  geom_text(data = Summary_Combined_sum_loc_hr_bins %>%group_by(ID) %>%slice_max(DOSE, n = 1), # ID labels
     aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "Dose", y = "All movement (m)", title = "0 - 24 hrs post")
plot_loc_0_24hr_median

#Export plot to folder called "orexin_A" 
ggsave(plot_loc_0_24hr_median,
       filename="loc_0_24hr_median.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")


#For physical activity I could look at the percent of their time spent moving as
#well as cumulative movement over the course of 24hrs after receiving each dose...
#this could be a way of looking at movement with the context of what the mouse was previously doing

#graph of EE vs BW...if they are strongly correlated then that could explain
#why I see two different groups for EE. Could also graph EE vs physical activity

# Additional analyses ####

##EE (kcal/hr) vs BW (g) ####

##Hrs 0-2 ####
#### GRAPH: EE (kcal/hr) vs BW (Hrs 0-2) ####
plot_EEvsBW_0_2hr <- ggplot(Summary_Combined_avg_EE_hr_bins, 
                                 aes(x = BW, y = Avg_EE_0_2hr, color = DOSE)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_text(aes(label = ID), hjust = -0.5, vjust = 0.5, size = 3) +
  #geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = custom_colors_OXA) +
  theme_classic() +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  theme(legend.position = "right") +
  labs(title="0 - 2 hrs post : EE vs BW",
    x = "Body Weight (g)",
    y = "Energy expenditure (kcal/hr)",
    color = "Dose") 
  #facet_wrap(~DOSE) 
plot_EEvsBW_0_2hr

#Export plot to folder called "orexin_A" 
ggsave(plot_EEvsBW_0_2hr,
       filename="EEvsBW_0_2hr.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/orexin_A")
  
####  Correlation: EE vs BW (Hrs 0-2) ####
#Pearson correlation
correlation_by_dose <- Summary_Combined_avg_EE_hr_bins %>%
  group_by(DOSE) %>%
  summarise(n = sum(complete.cases(BW, Avg_EE_0_2hr)),
    correlation_r = cor(BW, Avg_EE_0_2hr, use = "complete.obs", method = "pearson"),
    p_value = cor.test(BW, Avg_EE_0_2hr, method = "pearson")$p.value, .groups = "drop") 
correlation_by_dose

#Spearman correlation
correlation_by_dose_spearman <- Summary_Combined_avg_EE_hr_bins %>%
  group_by(DOSE) %>%
  summarise(n = sum(complete.cases(BW, Avg_EE_0_2hr)),
    correlation_rho = cor(BW, Avg_EE_0_2hr, use = "complete.obs", method = "spearman"),
    p_value = cor.test(BW, Avg_EE_0_2hr, method = "spearman")$p.value, .groups = "drop")
correlation_by_dose_spearman


##Hrs 0-24 ####
#### GRAPH: EE (kcal/hr) vs BW (Hrs 0-2) ####
ggplot(Summary_Combined_avg_EE_hr_bins, aes(x = BW, y = Avg_EE_0_24hr, color = DOSE)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_text(aes(label = ID), hjust = -0.5, vjust = 0.5, size = 3) +
  #geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = custom_colors_OXA) +
  theme_classic() +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  theme(legend.position = "right") +
  labs(title="0 - 24 hrs post : EE vs BW",
    x = "Body Weight (g)",
    y = "Energy expenditure (kcal/hr)",
    color = "Dose") 
  #facet_wrap(~DOSE)

####  Correlation: EE vs BW (Hrs 0-2) ####
#Pearson correlation
correlation_by_dose <- Summary_Combined_avg_EE_hr_bins %>%
  group_by(DOSE) %>%
  summarise(n = sum(complete.cases(BW, Avg_EE_0_2hr)),
    correlation_r = cor(BW, Avg_EE_0_2hr, use = "complete.obs", method = "pearson"),
    p_value = cor.test(BW, Avg_EE_0_2hr, method = "pearson")$p.value, .groups = "drop") 
correlation_by_dose

## EE vs physical activity ####
  # I think a correlation between these two variables would suggest that the increase 
  #in EE is due to the increase in physical activity, but I'm not certain

#Combine dataframes: Summary_Combined_avg_EE_hr_bins and Summary_Combined_sum_loc_hr_bins
Summary_Combined_sum_loc_avg_EE_hr_bins <- Summary_Combined_avg_EE_hr_bins %>%
left_join(
    Summary_Combined_sum_loc_hr_bins %>% 
      select(ID, PERIOD, DOSE, sum_loc_0_2hr, sum_loc_2_4hr, sum_loc_0_4hr, sum_loc_0_24hr), 
    by = c("ID", "PERIOD", "DOSE")) %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  filter(!(ID==3740)) #3740 is an outlier at 250pmol for locomotion, so remove it from EE df too

##Hrs 0-2 ####
#### GRAPH: EE vs All movement ####
ggplot(Summary_Combined_sum_loc_avg_EE_hr_bins, aes(x = sum_loc_0_2hr, y = Avg_EE_0_2hr, color = DOSE)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_text(aes(label = ID), hjust = -0.5, vjust = 0.5, size = 3) +
  #geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = custom_colors_OXA) +
  theme_classic() +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  theme(legend.position = "right") +
  labs(title="0 - 2 hrs post : EE vs Activity",
    x = "All movement (m)",
    y = "Energy expenditure (kcal/hr)",
    color = "Dose") +
facet_wrap(~DOSE) 

#### Correlation: EE vs All movement ####
#Spearman correlation
correlation_by_dose_spearman <- Summary_Combined_sum_loc_avg_EE_hr_bins %>%
  group_by(DOSE) %>%
  summarise(n = sum(complete.cases(sum_loc_0_2hr, Avg_EE_0_2hr)),
    correlation_rho = cor(sum_loc_0_2hr, Avg_EE_0_2hr, use = "complete.obs", method = "spearman"),
    p_value = cor.test(sum_loc_0_2hr, Avg_EE_0_2hr, method = "spearman")$p.value, .groups = "drop")
correlation_by_dose_spearman


