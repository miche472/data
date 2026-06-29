#Goal: identify the best way to calculate RMR and NEAT

#Script started: 3-17-26
#Latest revision: 3-17-26

#New complete day method from 3-16 & 3-17-26 was the impetus for this script

#The way that I worked through this is in script: "10th percentile method for RMR (start 3-17-26)"
#This version, "10th percentile method for RMR (start 3-17-26) updated") is annotated & cleaned for clarity

#----------------------------------#
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
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))}

#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#This script uses df = filter_locom_energy3 as the source data frame. Code is at bottom of script

#----------------------------------#
#Should I use summation method for EE or mean method? 
#Checking coverage --> this tells me how many hour have less than 80% of the possible 60 observations
#Across all mice and time points I have 1917 observation hours that are part of 
#24hr perfect complete days (i.e. my method for complete day). Of these, only 84 have less than 80%
#of 60 possible observations (i.e.hrs when more than 20% of minutes were not recorded). 
#33 of these hours are from hr 18 (i.e. when recording was stopped for feeding, injections, and BW)
#Also, of 80 total "hour 18s", only 30 have 100% coverage.Perhaps I should remove hr 18 from analyses?
#In general though, summation is a better route and the amount of missing data should 
#isn't great enough to force me to use means
ID_TEE3_hr_coverage <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr) %>%
  summarise(coverage = n() / 60) %>%
filter(!coverage >= 0.8)

# In previoius scripts I used whether or not the mouse moved between two minutes to classify
#the EE during that minute as entirely NEAT or entirely RMR. This results in an 
    #overestimate of NEATsince even when I mouse is moving some of its TEE is RMR.
#A way to get around this issue is to calculate the 10th percentile of TEE for a given 
    #observation period (hr, lights on/off, or daily [aka global])--> for now, do by hr
#For min. when mouse didn't move: RMR=TEE and NEAT=0. 
#For min. when mouse moved, NEAT = TEE - (RMR=10th percentile of RMR values across that hr for the ID)

#---------------------------------------------------------------------#
# RMR calculated using percentile Percentile for RMR --> 
    #happy with this...just check if global or light/dark RMR is better ####
Minute_sum_EE_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
# Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, GROUP) %>%
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
#Re-attach GROUP
mutate(GROUP = case_when(
  ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
  ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Calculate daily sum (add all hours together)
Daily_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day))) 

#Graph RAW hourly TEE, NEAT, and RMR ####
  #use df = Minute_sum_EE_hr
#---Graph raw hourly TEE (kcal/hr)
ggplot(Minute_sum_EE_hr, aes(x = hr, y = TEE_kcal, fill = GROUP)) +
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
  facet_wrap(~SABLE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title="TEE (kcal/hour)--Raw",
    x = "Hour",
    y = "TEE (kcal/hour)",
    fill = "Treatment group")

#---Graph raw hourly RMR (kcal/hr) ####
ggplot(Minute_sum_EE_hr, aes(x = hr, y = RMR_kcal, fill = GROUP)) +
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
  facet_wrap(~SABLE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title="RMR (kcal/hour)--Raw",
    x = "Hour",
    y = "RMR (kcal/hour)",
    fill = "Treatment group")

#---Graph raw hourly NEAT (kcal/hr) ####
ggplot(Minute_sum_EE_hr, aes(x = hr, y = NEAT_kcal, fill = GROUP)) +
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
  facet_wrap(~SABLE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title="NEAT (kcal/hour)--Raw",
    x = "Hour",
    y = "NEAT (kcal/hour)",
    fill = "Treatment group")

#---
#Graph RAW daily TEE, NEAT, and RMR ####
#use df = Daily_EE
#---Graph raw daily TEE (kcal/day) ####
ggplot(Daily_EE, aes(x = SABLE, y = TEE_kcal_day, fill = GROUP)) +
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
    title="TEE (kcal/day)--Raw",
    x = "Time point",
    y = "TEE (kcal/day)",
    fill = "Treatment group")

#---Graph raw daily NEAT (kcal/day) ####
ggplot(Daily_EE, aes(x = SABLE, y = NEAT_kcal_day, fill = GROUP)) +
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
    title="NEAT (kcal/day)--Raw",
    x = "Time point",
    y = "NEAT (kcal/day)",
    fill = "Treatment group")

##---------------------------------------------####
#Multiple linear regression and bar plots for EE Summation method

#Process/format echoMRI data ####
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #tech issue with sable cages
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(day_rel = Date - first(Date),
         STATUS = case_when(
           n_measurement == 1 ~ "Baseline",
           Date %in% as.Date(c("2025-01-27", "2025-01-29", "2025-02-07")) ~ "Peak obesity",
           ID %in% c(3711, 3727) & Date == as.Date("2025-02-20") ~ "Peak obesity",
           Date == as.Date("2025-03-28") ~ "BW loss",
           Date == as.Date("2025-05-27") ~ "BW maintenance",
           Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                               "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
           TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3711 & Date == as.Date("2025-01-27"))) %>% #repeated echo
  filter(!(ID == 3727 & Date == as.Date("2025-02-07"))) %>% #repeated echo
  mutate(STATUS = factor(STATUS,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  rename(SABLE = STATUS) 

#Check number of mice in echoMRI_data
echoMRI_data %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#Combine echoMRI_data with energy expenditure by hour (TEE, RMR, NEAT)
Echo_Minute_sum_EE_hr <- Minute_sum_EE_hr  %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))

#Combine echoMRI_data with DAILY energy expenditure (TEE, RMR, NEAT)
Echo_Daily_EE <- Daily_EE  %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))


#----EE without adjustment for lean or BW (minute summation with 10th percentile RMR method)-----#
### NEAT not adj. ####
#Build multiple linear regression model for NEAT not adjusted for BW or lean #
model_NEAT_10 <- lmer(NEAT_kcal_day ~ SABLE * GROUP + (1 | ID), data = Daily_EE)
summary(model_NEAT_10)

#Calculate estimated marginal means #
emm_NEAT_10 <- emmeans(model_NEAT_10, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_10_df <- as.data.frame(emm_NEAT_10)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_10 <- contrast(emm_NEAT_10, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_10_df <- as.data.frame(contrasts_by_group_NEAT_10)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_10 <- contrast(emm_NEAT_10, method = "pairwise", by = "SABLE")
contrasts_SABLE_NEAT_10_df <- as.data.frame(contrasts_by_SABLE_NEAT_10)

# Bar plot - Graph predicted NEAT not adjusted for BW or lean #
barplot_emm_NEAT_10 <- emm_NEAT_10_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "NEAT (kcal/day) unadjusted",
    y = "NEAT (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold", size=15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11))
barplot_emm_NEAT_10
#Results: At BW regain, control and weight cycled mice almost had significantly different NEAT
  #linear mixed model has p= 0.073

###-----RMR not adj. -----####
#unadjusted RMR: Build multiple linear regression model for RMR not adjusted for BW or lean #
#Build multiple linear regression model
model_RMR_10 <- lmer(RMR_kcal_day ~ SABLE * GROUP + (1 | ID), data = Daily_EE)
summary(model_RMR_10)

# Calculate estimated marginal means #
emm_RMR_10 <- emmeans(model_RMR_10, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_10_df <- as.data.frame(emm_RMR_10)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_10 <- contrast(emm_RMR_10, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_10_df <- as.data.frame(contrasts_by_group_RMR_10)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_10 <- contrast(emm_RMR_10, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RMR_10_df <- as.data.frame(contrasts_by_SABLE_RMR_10)

# Bar plot - Graph predicted RMR unadjusted for lean or BW #
barplot_emm_RMR_10 <- emm_RMR_10_df %>%
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
    title = "RMR (kcal/day) unadjusted",
    y = "RMR (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold", size=15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11))
barplot_emm_RMR_10
#Conclusion: At BW regain, control mice have significantly higher 
    #RMR compared to weight cycled mice p=0.0056

###-----TEE not adj.-----####
#unadjusted TEE: Build multiple linear regression model for TEE not adjusted for BW or lean #
#Build multiple linear regression model
model_TEE_10 <- lmer(TEE_kcal_day ~ SABLE * GROUP + (1 | ID), data = Daily_EE)
summary(model_TEE_10)

# Calculate estimated marginal means #
emm_TEE_10 <- emmeans(model_TEE_10, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_10_df <- as.data.frame(emm_TEE_10)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_10 <- contrast(emm_TEE_10, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_10_df <- as.data.frame(contrasts_by_group_TEE_10)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_10 <- contrast(emm_TEE_10, method = "pairwise", by = "SABLE")
contrasts_SABLE_TEE_10_df <- as.data.frame(contrasts_by_SABLE_TEE_10)

# Bar plot - Graph predicted TEE unadjusted for lean or BW #
barplot_emm_TEE_10 <- emm_TEE_10_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Total energy expenditure (kcal/day) unadjusted",
    y = "TEE (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold", size=15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11))
barplot_emm_TEE_10
#Results of MLR (unadjusted): At BW regain, control mice had significanlty higher TEE
#compared to weight cycled mice (p=0.00495)







#--------------------------------------------------------------#
#--------------------------------------------------------------#
# Code for df = filter_locom_energy_3 ####

##------- Get only the All_meters parameter --> value column has All_meters data
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
  group_by(ID) 

#Get complete day 1 and 2 by explicitly definining the zt_time and SABLE_DAY
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
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_loc3 <- filter_loc2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2))

#Check number of mice in filter_loc3
filter_loc3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

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
  group_by(ID) 

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
filter_EE2 <- filter_EE1 %>% 
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
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_EE3 <- filter_EE2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2)) #remove observations not from complete day 1 or 2

#Check number of mice in filter_EE3
filter_EE3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

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

#---
#Join filter_locom and filter_energy into a df called filter_locom_energy
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_locom_energy3 <- filter_locom3 %>%
  left_join(
    filter_energy3 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#T test for BW
mod_TEE_day_regain <- mod_TEE_day %>% 
  filter(SABLE=="BW regain") %>%
  ungroup() %>% 
  group_by(GROUP)

#Summary statistics
mod_TEE_day_regain %>%
  group_by(GROUP) %>%
  get_summary_stats(sum_avg_TEE_hr, type = "mean_sd")

#T test
res <- t.test(sum_avg_TEE_hr ~ GROUP, data = mod_TEE_day_regain)
res
#--------------------------------------------------------------#
#Idea --> scratched ####
#Calculate RMR for entire light and dark cycle (rather than each hr using 10th percentile method)
#This isn't a great idea actually because 8am and 6pm could throw off RMR for the whole light cycle.
#Also, the fact that the light and dark cycle are different durations complicates things.
ID_TEE_lightdark <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(SABLE, ID, GROUP, lights) %>%
  
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  
  ungroup() %>%
  group_by(SABLE, ID, lights) %>%
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE), #Calculate the sum of TEE during the dark and light cycle
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE), # RMR light/dark using 10th percentile
    RMR_kcal = RMR_rate * (n_obs / 60), # Total RMR energy in the light/dark cycle
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE), # NEAT:subtract RMR only during movement
    .groups = "drop") %>%
  group_by(ID, SABLE) %>%
  mutate(LD_TEE = sum(TEE_kcal))