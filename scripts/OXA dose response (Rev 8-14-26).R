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

#Format plot (LM version 3) ####
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
custom_colors_OXA <- c("aCSF"    = "#D9D9D9", 
                        "125pmol" = "#C7E9C0", 
                        "250pmol" = "#A1D99B", 
                        "250pmol" = "#74C476", 
                        "500pmol" = "#41AB5D", 
                        "1000pmol"= "#006D2C")
#custom_colors_RTIOXA <- c()

#***START DELTE*** ####
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

#***END DELETE*** ####

# USE (8-13-26): Alternative way to remove invalid injection for 3731 ####
##Read in Sable data ####
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
##Filter Sable data to include only mice that are cannulated ####
sable_dwn_19 <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()
  
##Read in meta data for injection time ####
read_injection_time <- read_csv("../data/META_INJECTIONS_LG.csv")

##Process injection time data ####
injection_time <- read_injection_time %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  #filter(VALID == "VALID") %>% #removes INVALID injection for 3731
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = row_number()) %>%
  ungroup()

##Identify the exact injection date/time which corresponds to each chunk of recording (i.e. Period of recording) ####
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
 
##Join injection time and sable data ####
 OXA_Sable_joined <- sable_dwn_19 %>%
  left_join(
    injection_assignment,
    by = c("ID", "DateTime"))
 
##Remove INVALID injections ####
#remove the first attempt at 3731, 500pmol
Sable_OXA <- OXA_Sable_joined %>%
   filter(!(VALID == "INVALID")) %>% #removes the invalid injection (3731)
   group_by(ID, PERIOD) %>%
   mutate(INJ_PERIOD = if_else(ID=="3731" & PERIOD==6, 5, PERIOD)) #assigns the valid injection for 3731 as injection 5 

##Calculate time after injection (in minutes and in hours) ####
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

#Add BW data ####
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3731,3732,3733,3735,3737,3738,3739,3740,3741))%>%
  mutate(PERIOD = case_when(ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-17" ~ "1",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-20" ~ "2",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-24" ~ "3",
                             ID %in% c(3731,3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-27" ~ "4",
                             ID %in% c(3732,3733,3735,3737,3739,3740,3741) & DATE =="2026-07-29" ~ "5",
                             ID==3731 & DATE =="2026-08-03" ~ "5",
                             ID==3738 & DATE =="2026-08-03" ~ "1",
                             ID==3738 & DATE =="2026-08-11" ~ "2")) %>%
  drop_na(PERIOD) %>%
  mutate(ID=as.factor(ID),
         PERIOD=as.integer(PERIOD))
  
#Calculate sum of TEE in 60 min bins post start of recording #### 
OXA_loc_EE_60min <- OXA_loc_EE %>%
  group_by(ID, PERIOD, DateTime) %>% #group by DOSE or PERIOD?
  arrange(DateTime) %>%
  mutate(TEE_per_min = Kcal_Hr/60) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  group_by(ID, PERIOD, recording_bin) %>%
  summarise(
    TEE_kcal = sum(TEE_per_min, na.rm = TRUE),
    .groups = "drop") %>%
  left_join(
    injection_time %>% 
      select(DOSE, ID, PERIOD), 
    by = c("PERIOD", "ID")) %>%
mutate(DOSE = factor(DOSE,levels = c("aCSF", "125pmol", "250pmol", "500pmol", "1000pmol"))) 

##Graph the sum of EE in 1st hour for each mouse during each period ####
plot_hr_1 <- OXA_loc_EE_60min %>%
  filter(recording_bin==0)

ggplot(plot_hr_1, aes(x=DOSE, y=TEE_kcal)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Graph the sum of EE in 2nd hour for each mouse during each period ####
plot_hr_2 <- OXA_loc_EE_60min %>%
  filter(recording_bin==1)

ggplot(plot_hr_2, aes(x=DOSE, y=TEE_kcal)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Graph the sum of EE in 3rd hour for each mouse during each period ####
plot_hr_3 <- OXA_loc_EE_60min %>%
  filter(recording_bin==2)

ggplot(plot_hr_3, aes(x=DOSE, y=TEE_kcal)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Graph the sum of EE in the first 2hrs for each mouse during each period ####
plot_2hrs_sum <- OXA_loc_EE_60min %>%
  filter(recording_bin %in% c(0,1)) %>%
  group_by(ID, DOSE) %>%
  summarise(EE_2hrs= sum(TEE_kcal)) %>%
  ungroup() %>%
  left_join(
    injection_time %>% 
      select(DOSE, ID, PERIOD), 
    by = c("DOSE", "ID")) %>%
  mutate(DOSE = factor(DOSE,levels = c("aCSF", "125pmol", "250pmol", "500pmol", "1000pmol"))) 

ggplot(plot_2hrs_sum, aes(x=DOSE, y=EE_2hrs)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

## Stat analysis for sum of EE ####
#Linear mixed model which has EE as outcome and includes DOSE, PERIOD, and ID
#Build multiple linear regression model for RMR not adjusted for BW or lean #
Stats_2hrs_sum <- plot_2hrs_sum %>%
left_join(BW_COHORT19 %>% 
      select(ID, PERIOD, BW), 
    by = c("ID", "PERIOD")) %>%
  mutate(DOSE = factor(DOSE,levels = c("aCSF", "125pmol", "250pmol", "500pmol", "1000pmol"))) 

model_EE_2hrs_sum <- lmer(EE_2hrs ~ DOSE + PERIOD + BW + (1 | ID), data = Stats_2hrs_sum)
summary(model_EE_2hrs_sum)

emmeans(model_EE_2hrs_sum, pairwise ~ DOSE)

# Calculate estimated marginal means #
emm_change_BW_g <- emmeans(model_EE_2hrs_sum, ~ Transition * GROUP, cov.reduce = mean)
emm_change_BW_g_df <- as.data.frame(emm_change_BW_g)

# Pairwise contrasts within each GROUP
contrasts_by_group_change_BW_g <- contrast(emm_change_BW_g, method = "pairwise", by = "GROUP")
contrasts_by_group_change_BW_g_df <- as.data.frame(contrasts_by_group_change_BW_g)

# Pairwise contrasts within each stage of weight cycling
contrasts_by_SABLE_change_BW_g <- contrast(emm_change_BW_g, method = "pairwise", by = "Transition")
contrasts_SABLE_change_BW_g_df <- as.data.frame(contrasts_by_SABLE_change_BW_g)


#Calculate avg. TEE (kcal/hr) in 60 min bins post start of recording #### 
OXA_loc_EE_60min_avg <- OXA_loc_EE %>%
  group_by(ID, DOSE, DateTime) %>%
  arrange(DateTime) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60)) %>%
  group_by(ID, DOSE, recording_bin) %>%
  summarise(avg_kcal_per_hr = mean(Kcal_Hr))

##Graph avg EE (kcal/hr) during 1st hour for each mouse during each period ####
plot_hr_1_avg <- OXA_loc_EE_60min_avg %>%
  filter(recording_bin==0)

ggplot(plot_hr_1_avg, aes(x=DOSE, y=avg_kcal_per_hr)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Graph avg EE (kcal/hr) during 2nd hour for each mouse during each period ####
plot_hr_2_avg <- OXA_loc_EE_60min_avg %>%
  filter(recording_bin==1)

ggplot(plot_hr_2_avg, aes(x=DOSE, y=avg_kcal_per_hr)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Graph avg EE (kcal/hr) during the first 2hrs for each mouse during each period ####
plot_2hrs_avg <- OXA_loc_EE_60min_avg %>%
  filter(recording_bin %in% c(0,1))

ggplot(plot_2hrs_avg, aes(x=DOSE, y=avg_kcal_per_hr)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  facet_wrap(~ID)

##Other EE graphs ####
#average hourly EE vs time 
OXA_loc_EE_60min_bins <- OXA_loc_EE %>%
  group_by(ID, DOSE, DateTime) %>%
  arrange(DateTime) %>%
  mutate(
    recording_bin = floor(minutes_post_recording / 60)) %>%
  filter(ID==3731) %>%
  filter(recording_bin <4)

plot_time_course_3731 <- ggplot(OXA_loc_EE_60min_bins, aes(x = DOSE, y = Kcal_Hr, fill = recording_bin, color=recording_bin)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) + facet_wrap(~recording_bin)
plot_time_course_3731

#Sum EE graph
#EE vs dose with one line for each mouse
ggplot(
  Stats_2hrs_sum,
  aes(x = DOSE, y = EE_2hrs, group = ID)
) +
  geom_line(alpha = 0.4) +
  geom_point() +
  theme_classic()

#---------------------------------.
#---------------------------------.
#---------------------------------.

# LOCOMOTION ####
##All_Meters ####
All_Meters_0hr_2hr <- OXA_loc_EE %>%
  ungroup()%>%
  filter(minutes_post_recording <121) %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Move_m= All_meters-lag(All_meters)) %>%
  drop_na(Move_m) %>%
  ungroup() %>%
  group_by(DOSE,ID) %>%
  summarise(Sum_move_m = sum((Move_m)))

#By ID
ggplot(All_Meters_0hr_2hr, aes(x=DOSE, y=Sum_move_m, color=DOSE)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  format.plot_LM3 +
    scale_color_manual(values = custom_colors_OXA) +
  facet_wrap(~ID)

#average for each dose
#remove outlier of 3740 at 250pmol
plot_All_Meters_0hr_2hr <- All_Meters_0hr_2hr %>%
  filter(!(ID==3740))
  
ggplot(plot_All_Meters_0hr_2hr, aes(x=DOSE, y=Sum_move_m, color=DOSE, group=DOSE, fill=DOSE)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
 geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
    scale_color_manual(values = custom_colors_OXA) +
  scale_fill_manual(values = custom_colors_OXA) +
  theme_bw(base_size = 14) +
 format.plot_LM3 

###Distance ####
OXA_loc_60min <- OXA_loc_EE %>%
  group_by(ID, DOSE, DateTime) %>%
  arrange(DateTime) %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Move_m= All_meters-lag(All_meters)) %>%
  drop_na(Move_m) %>%
  ungroup() %>%
  group_by(DOSE, ID) %>%
  arrange(DateTime) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  group_by(ID, DOSE, recording_bin) %>%
  summarise(Sum_move_m = sum((Move_m)))

#Example graph
plot_OXA_loc_60min <- OXA_loc_60min %>%
  #filter(!(ID==3740)) %>%
  filter(recording_bin <2)

plot_time_course <- ggplot(plot_OXA_loc_60min, aes(x = DOSE, y = Sum_move_m, fill = recording_bin, color=recording_bin)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2)  +
  facet_wrap(~recording_bin)
plot_time_course

###Percent of time spent moving ####

##Ped_Meters ####
###Distance ####
OXA_Ped_60min <- OXA_loc_EE %>%
  group_by(ID, DOSE, DateTime) %>%
  arrange(DateTime) %>%
  group_by(ID, DOSE) %>%
  arrange(DateTime) %>%
  mutate(Ped_m= Ped_meters-lag(Ped_meters)) %>%
  drop_na(Ped_m) %>%
  ungroup() %>%
  group_by(DOSE, ID) %>%
  arrange(DateTime) %>%
  mutate(recording_bin = floor(minutes_post_recording / 60)) %>%
  group_by(ID, DOSE, recording_bin) %>%
  summarise(Sum_ped_m = sum((Ped_m))) %>%
  filter(recording_bin %in% c(0,1)) %>%
  group_by(ID, DOSE) %>%
  summarise(first_2hrs_Ped_m = sum(Sum_ped_m)) %>%
  filter(!(ID==3740))

ggplot(OXA_Ped_60min, aes(x=DOSE, y=first_2hrs_Ped_m, color=DOSE, group=DOSE, fill=DOSE)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
 geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
    scale_color_manual(values = custom_colors_OXA) +
  scale_fill_manual(values = custom_colors_OXA) +
  theme_bw(base_size = 14) +
 format.plot_LM3 

###Percent of minutes during which the mouse moved each hour ####

##Cumulative sum of locomotion (for distance and time) ####

# What is energy expenditure and BW for these mice during the baseline recording period? ####
#Baseline: use 6:30pm on 2026-07-16 through 6:30pm on 2026-07-17 
#1. Compare aCSF to baseline at the same time of day --> does handling have a really big impact?
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

#functions
zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

#EE during baseline period
filter_TEE1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date %in% c("2026-07-16", "2026-07-17")) %>% 
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                                  min(DateTime),
                                                  units="hours")),
         zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()

#Locomotion during baseline period
#AllMeters
filter_Allloc1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date %in% c("2026-07-16", "2026-07-17")) %>% 
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                                  min(DateTime),
                                                  units="hours")),
         zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()

#PedMeters
filter_Pedloc1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date %in% c("2026-07-16", "2026-07-17")) %>% 
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off")) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                                  min(DateTime),
                                                  units="hours")),
         zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()

#Process EE and loc to join
filter_Allloc2 <- filter_Allloc1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_Pedloc2 <- filter_Pedloc1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

filter_TEE2 <- filter_TEE1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#Join loc and EE
filter_loc_TEE2 <- filter_Allloc2 %>%
  left_join(
    filter_TEE2 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx?
    by = c("ID", "DateTime", "sable_idx")) %>%
  left_join(
    filter_Pedloc2 %>% 
      select(Ped_meters, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx?
    by = c("ID", "DateTime", "sable_idx"))

# TEE and RMR calculated using percentile for RMR --> 
Minute_sum_EE_hr <- filter_loc_TEE2 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(ID) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(ID, hr, lights) %>% # group by lights to retain photoperiod in df for later use
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    # Total Energy Expenditure (TEE) --> minute by minute summation
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    RMR_kcal = RMR_rate * (n_obs / 60),
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE), .groups = "drop")

Daily_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))

#actually, try extracting the time into a separate column for each ID. then create a column called 
#start recording which is equal to the recording start time for aCSF day?
