#Orexin A icv dose response (July/August 2026)

#Started: 7-30-26
#Revised: 8-4-26
#Need to work on standardizing post injection time so that I am looking at the same number of recorded minutes per mouse. 
#Could look at first 120 minutes of recording for each mouse rather than whether an observation was within 120 minutes of injection time
# Could also look at average minute by minute EE...I think there would be huge error using this approach though

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

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


pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand
)

#Format plot (LM version 2) ####
format.plot_LM2 <- theme(
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
custom_colors1_GLP_Exp2 <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")
custom_colors2_GLP_Exp2 <- c("Tirzepatide" = "#1E90FF", "Vehicle" = "#8B8989")
custom_colors3_GLP_Exp2 <- c("Tirzepatide" = "#104E8B", "Vehicle" = "#403d3c")

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

# 0.5-2.5hrs after injection...cumulative
# 0.5-4.5hrs after injection...cumulative 
# 0.5-24.5hrs after injection...cumulative

#TEE, RMR and SPA
#Calculate hourly basis --> use photoperiod code from Exp 2a Photo EE & RQ

#Locomotion: time spent moving and distance traveled
#Look at resulting formatting and decide how to append dose information. Could hardcode, but there may be a better approach


 
 
 
 

 
 


#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
 #Filter Sable data
sable_dwn_19 <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()

#Should I remove this differently? Maybe assign valid or invalid based on time and then remove any "invalids" ####
#Need to get rid of sable data for 3731 after the invalid injection. 
sable_dwn_OXA <- sable_dwn_19 %>%
  group_by(ID, DateTime) %>%
  arrange(DateTime) %>%
  #filter(ID== "3731" & DateTime > "2026-07-30 16:07:00" & DateTime < "2026-08-03 11:34:00") #Use this to see how many lines are recorded during this period
  filter(!(ID== "3731" & DateTime > "2026-07-30 16:07:00" & DateTime < "2026-08-03 11:34:00")) %>%
  ungroup()
  
#Read in meta data for injection time
read_injection_time <- read_csv("../data/META_INJECTIONS_LG.csv")

#Process injection time data
injection_time <- read_injection_time %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  filter(VALID == "VALID") %>% #removes INVALID injection for 3731
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = row_number()) %>%
  ungroup()

#Identify the exact injection date/time which corresponds to each chunk of recording (i.e. Period of recording)
#This data frame has one row for each DateTime +ID...it does not yet account for the fact that multiple parameters (10)
#were measured at once for each mouse
 injection_assignment <- sable_dwn_OXA %>%
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
 
#Join injection sign posts and sable data
OXA_Sable_joined <- sable_dwn_OXA %>%
  left_join(
    injection_assignment,
    by = c("ID", "DateTime")) %>%
  filter(PERIOD %in% c("1","2","3","4","5")) #Removes observations from before first injection

#---
#---
#Calculate time after injection (in minutes and in hours)
 #Make a column for whether an observation is within 24hrs of an injection --> this is a complete day within my experimental paradigm
OXA_Sable_joined_2 <- OXA_Sable_joined %>%
  ungroup() %>%
  group_by(ID, PERIOD) %>%
  arrange(DateTime) %>%
  mutate(minutes_post_injection = as.numeric(difftime(DateTime,INJECTION_DateTime,units = "mins")),
         
  first_recording = min(DateTime, na.rm = TRUE), #After first actual recorded observation
  minutes_post_recording = as.numeric( difftime(DateTime, first_recording, units = "mins")),
  hrs_post_recording = floor(minutes_post_recording / 60) + 1)


 

#Creat seperate dfs with locomotion and with energy expenditure

OXA_EE_1 <-OXA_Sable_joined_2 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

OXA_loc_1 <-OXA_Sable_joined_2 %>%
    mutate(
      Time = as_hms(format(DateTime, "%H:%M:%S")),
      lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

OXA_ped_1 <-OXA_Sable_joined_2 %>%
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

#Energy expenditure for observations 24hrs after injection ####
OXA_Minute_sum_EE_hr <- OXA_loc_EE %>%
  filter(minutes_post_recording <1440) %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(ID, PERIOD) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(DOSE, ID) %>% # group by lights to retain photoperiod in df for later use
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
  mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal)) 

#Calculate 24hr sum (add all hours together)
OXA_24hr_EE <- OXA_Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  summarise(TEE_kcal_24hr = sum(TEE_kcal),
            NEAT_kcal_24hr = sum(NEAT_kcal),
            RMR_kcal_24hr = sum(RMR_kcal),
            diff= abs(TEE_kcal_24hr - sum(NEAT_kcal_24hr + RMR_kcal_24hr)))

#Calculate mean 24hr by dose
DOSE_OXA_24hr_EE <- OXA_24hr_EE %>%
  ungroup() %>%
  group_by(DOSE) %>%
  summarise(Avg_TEE_kcal_24hr = mean(TEE_kcal_24hr),
            Avg_NEAT_kcal_24hr = mean(NEAT_kcal_24hr),
            Avg_RMR_kcal_24hr = mean(RMR_kcal_24hr),
            Avg_diff= mean(diff))

#4hrs post injection ####
OXA_4hr_Minute_sum_EE_hr <- OXA_loc_EE %>%
  filter(minutes_post_recording < 241) %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(ID, PERIOD) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(DOSE, ID, hr) %>% # group by lights to retain photoperiod in df for later use
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
  mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal)) 

#Calculate 4hr sum (add all 4 hours together)
OXA_4hr_EE <- OXA_4hr_Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  summarise(TEE_kcal_4hr = sum(TEE_kcal),
            NEAT_kcal_4hr = sum(NEAT_kcal),
            RMR_kcal_4hr = sum(RMR_kcal),
            diff= abs(TEE_kcal_4hr - sum(NEAT_kcal_4hr + RMR_kcal_4hr)))

#Calculate mean 4hr by dose
DOSE_OXA_4hr_EE <- OXA_4hr_EE %>%
  ungroup() %>%
  group_by(DOSE) %>%
  summarise(Avg_TEE_kcal_4hr = mean(TEE_kcal_4hr),
            Avg_NEAT_kcal_4hr = mean(NEAT_kcal_4hr),
            Avg_RMR_kcal_4hr = mean(RMR_kcal_4hr),
            Avg_diff= mean(diff))

# Try different approach for calculating TEE ####
OXA_2hr_Minute_sum_EE_hr <- OXA_loc_EE %>%
  filter(minutes_post_recording <121) %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(ID, PERIOD) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(DOSE, ID, hrs_post_recording) %>% # group by lights to retain photoperiod in df for later use
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
  mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal))

#Calculate 2hr sum (add all 2 hours together)
OXA_2hr_EE <- OXA_2hr_Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  summarise(TEE_kcal_2hr = sum(TEE_kcal),
            NEAT_kcal_2hr = sum(NEAT_kcal),
            RMR_kcal_2hr = sum(RMR_kcal),
            diff= abs(TEE_kcal_2hr - sum(NEAT_kcal_2hr + RMR_kcal_2hr)))

#Calculate mean 2hr EE by dose
DOSE_OXA_2hr_EE <- OXA_2hr_EE %>%
  ungroup() %>%
  group_by(DOSE) %>%
  summarise(Avg_TEE_kcal_2hr = mean(TEE_kcal_2hr),
            Avg_NEAT_kcal_2hr = mean(NEAT_kcal_2hr),
            Avg_RMR_kcal_2hr = mean(RMR_kcal_2hr),
            Avg_diff= mean(diff))

#Graph TEE (x axis is minute_post_recording and y axis is EE) ####

OXA_loc_EE_4hr <- OXA_loc_EE %>%
  filter(minutes_post_recording <241) %>%
arrange(ID, DateTime) %>%
  group_by(ID, DOSE) %>%
  mutate(
    EE_rolling_15min = rollmean(
      Kcal_Hr,
      k = 15,
      fill = NA,
      align = "right"
    )
  ) %>%
  ungroup()
#LEFT off: 
##Try a rolling average every 15 minutes and graph that in units of kcal/hr ####
##Try not removing the time immediately before recording so that i can see if there ####
#is a spike when injections are given

ggplot(OXA_loc_EE_4hr, aes(x=minutes_post_recording, y=EE_rolling_15min, group=DOSE, fill=DOSE, color=DOSE)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  #geom_errorbar(stat = "summary", 
                #fun.data = mean_se, 
               # position = position_dodge(width = 0.8), 
                #width = 0.25, linewidth = 0.65, color="black") + #454441
  #geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
            # alpha = 0.6, size = 2) +
  facet_wrap(~DOSE) 

## Try graphing the summation of EE every hr for 24 hours. May want to compare light and dark period during baseline vs after injection ####


#---
#---
#All_Meters - All physical activity ####
## 2hrs post injection ####
OXA_2hr_Minute_sum_All_Meters <- OXA_loc_EE %>%
  filter(minutes_post_recording <121) %>%
  ungroup() %>%
  group_by(ID, PERIOD, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Movement_m= All_meters-lag(All_meters)) %>%
  drop_na(Movement_m) %>%
  ungroup() %>%
  group_by(DOSE,ID) %>%
  summarise(Total_move_m = sum((Movement_m))) %>%
  ungroup() %>%
  group_by(DOSE) %>%
  summarise(avg_total_m = mean(Total_move_m))

## 4hrs post injection ####
OXA_4hr_Minute_sum_All_Meters <- OXA_loc_EE %>%
  filter(minutes_post_recording <241) %>%
  ungroup() %>%
  group_by(ID, PERIOD, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Movement_m= All_meters-lag(All_meters)) %>%
  drop_na(Movement_m) %>%
  ungroup() %>%
  group_by(DOSE,ID) %>%
  summarise(Total_move_m = sum((Movement_m))) %>%
  ungroup() %>%
  group_by(DOSE) %>%
  summarise(avg_total_m = mean(Total_move_m))

#Try calculating Ped_Meters for 0-60 minutes post start of recording
OXA_1hr_Minute_sum_All_Meters <- OXA_loc_EE %>%
  ungroup()%>%
  filter(minutes_post_recording <121) %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Movement_m= All_meters-lag(All_meters)) %>%
  drop_na(Movement_m) %>%
  ungroup() %>%
  group_by(DOSE,ID) %>%
  summarise(Total_move_m = sum((Movement_m)))

#Graph: x axis is DOSE, y axis is Ped_meters, facet_wrap(~ID)

ggplot(OXA_1hr_Minute_sum_All_Meters, aes(x=DOSE, y=Total_move_m)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  #geom_errorbar(stat = "summary", 
                #fun.data = mean_se, 
               # position = position_dodge(width = 0.8), 
                #width = 0.25, linewidth = 0.65, color="black") + #454441
  #geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
            # alpha = 0.6, size = 2) +
  facet_wrap(~ID) 


#COME back here after running macros ####
## Alternative way to remove invalid injection for 3731 ####
#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
 #Filter Sable data
sable_dwn_19 <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()

##Should I remove this differently? Maybe assign valid or invalid based on time and then remove any "invalids" ####
#Need to get rid of sable data for 3731 after the invalid injection. 
#sable_dwn_OXA <- sable_dwn_19 %>%
  #group_by(ID, DateTime) %>%
  #arrange(DateTime) %>%
  #filter(ID== "3731" & DateTime > "2026-07-30 16:07:00" & DateTime < "2026-08-03 11:34:00") #Use this to see how many lines are recorded during this period
  #filter(!(ID== "3731" & DateTime > "2026-07-30 16:07:00" & DateTime < "2026-08-03 11:34:00")) %>%
 # ungroup()
  
#Read in meta data for injection time
read_injection_time <- read_csv("../data/META_INJECTIONS_LG.csv")

#Process injection time data
injection_time <- read_injection_time %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  #filter(VALID == "VALID") %>% #removes INVALID injection for 3731
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = row_number()) %>%
  ungroup()

#Identify the exact injection date/time which corresponds to each chunk of recording (i.e. Period of recording)
#This data frame has one row for each DateTime +ID...it does not yet account for the fact that multiple parameters (10)
#were measured at once for each mouse
 injection_assignment <- sable_dwn_OXA %>%
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
 
##I think this is where I should filter out invalids in order to remove the first attempt at 3731, 500pmol ####
 #Next step assigns periods...need to have processed files though
 
#Join injection sign posts and sable data
OXA_Sable_joined <- sable_dwn_OXA %>%
  left_join(
    injection_assignment,
    by = c("ID", "DateTime")) %>%
  filter(PERIOD %in% c("1","2","3","4","5")) #Removes observations from before first injection

#Calculate time after injection (in minutes and in hours)
#Create a column which looks at 20min post injection
OXA_Sable_joined_2 <- OXA_Sable_joined %>%
  ungroup() %>%
  group_by(ID, PERIOD) %>%
  arrange(DateTime) %>%
  mutate(minutes_post_injection = as.numeric(difftime(DateTime,INJECTION_DateTime,units = "mins")),
         
  first_recording = min(DateTime, na.rm = TRUE), #After first actual recorded observation
  minutes_post_recording = as.numeric( difftime(DateTime, first_recording, units = "mins")),
  hrs_post_recording = floor(minutes_post_recording / 60) + 1) %>%
  
  mutate(Post_injection_plus_20min = INJECTION_DateTime + lubridate::minutes(20))

