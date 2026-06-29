#Original Goal: identify the best way to calculate RMR and NEAT
#Subsequent inclusions: graphing Raw and estimates for TEE, RMR, and NEAT
#Linear mixed models for TEE, RMR, and NEAT with 10% method
#Photo period based TEE, RMR, and NEAT (calculation, graphs, linear mixed models)-> not completely polished

#Script started: 3-17-26
#Latest revision: 4-3-26

#New complete day method from 3-16 & 3-17-26 was the impetus for this script
#The way that I worked through this is in script: "10th percentile method for RMR (start 3-17-26)"
#"10th percentile method for RMR (start 3-17-26) updated") is annotated & cleaned for clarity

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
#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
custom_colors2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")
custom_colors3 <- c("on" = "#ffd145","off" = "#2c62a8")

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

#---Graph raw daily RMR (kcal/day) ####
ggplot(Daily_EE, aes(x = SABLE, y = RMR_kcal_day, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3, color="black") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  labs(
    title="RMR (kcal/day)--Raw",
    x = "Time point",
    y = "RMR (kcal/day)",
    fill = "Treatment group",
    color= "Treatment group")

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

##---------------------------------------------##
#Multiple linear regression and bar plots for EE Summation method ####

#Process/format echoMRI data 
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




#Try making box plots to visually look for outliers 
ggplot(emm_TEE_10_df, aes(x = SABLE, y = emmean, fill = GROUP)) + 
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")


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

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
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

#--------------------------------------------------------------#

#EE, RMR, and NEAT broken down by light and dark cycle
#cumulative EE using the 10th percentile method for RMR 

#Photo periods: EE (cumulative)-> TEE, NEAT, and RMR ####
#distinguishing between light and dark photo periods
Photo_EE <- Minute_sum_EE_hr %>%
  mutate(lights= if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% #add photo-period
  ungroup() %>%
  group_by(ID, SABLE, GROUP, lights) %>%
  summarise(TEE_kcal_photo = sum(TEE_kcal), 
            NEAT_kcal_photo = sum(NEAT_kcal),
            RMR_kcal_photo = sum(RMR_kcal))

#Sum light and dark photo periods to confirm that this matches the daily totals in 
#df=Daily_EE. Matches --> good, I can proceed
check_photo_EE <- Photo_EE %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(TEE_kcal_check = sum(TEE_kcal_photo), 
            NEAT_kcal_check = sum(NEAT_kcal_photo),
            RMR_kcal_check = sum(RMR_kcal_photo))
  
#Scale dark and light period so they are both 12 hrs in duration
#1) multiply dark by 1.2 (10hrs-->12hrs) and multiply light by 0.857 (14hrs->12hrs)
#2) sum original light and dark
#3) sum scaled light and dark
#4) calculate difference between original and scaled 

Scaled_Photo_EE <- Photo_EE %>%
  ungroup() %>%
  group_by(ID, SABLE, lights) %>%
  mutate(scaled_TEE = if_else(lights=="on", TEE_kcal_photo*0.8571429, TEE_kcal_photo*1.2),
         scaled_NEAT = if_else(lights=="on", NEAT_kcal_photo*0.8571429, NEAT_kcal_photo*1.2),
         scaled_RMR = if_else(lights=="on", RMR_kcal_photo*0.8571429, RMR_kcal_photo*1.2))

#Check that scaling doesn't have a massive impact on total daily TEE, RMR, and NEAT
check_Scaled_Photo_EE <- Photo_EE %>%
  ungroup() %>%
  group_by(ID, SABLE, lights) %>%
  mutate(scaled_TEE = if_else(lights=="on", TEE_kcal_photo*0.8571429, TEE_kcal_photo*1.2),
         scaled_NEAT = if_else(lights=="on", NEAT_kcal_photo*0.8571429, NEAT_kcal_photo*1.2),
         scaled_RMR = if_else(lights=="on", RMR_kcal_photo*0.8571429, RMR_kcal_photo*1.2)) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(TEE_kcal_check = sum(TEE_kcal_photo), 
            NEAT_kcal_check = sum(NEAT_kcal_photo),
            RMR_kcal_check = sum(RMR_kcal_photo),
            scaled_TEE_check = sum(scaled_TEE),
            scaled_NEAT_check=sum(scaled_NEAT),
            scaled_RMR_check =sum(scaled_RMR),
            TEE_original_scaled = abs(TEE_kcal_check - scaled_TEE_check),
            NEAT_original_scaled = abs(NEAT_kcal_check - scaled_NEAT_check),
            RMR_original_scaled = abs(RMR_kcal_check - scaled_RMR_check))



##---Graph raw EE by photo period (RAW EE) ####
###TEE ####
#### Graph: TEE light (RAW) ####
Light_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="on")

plot_raw_TEE_light <-ggplot(Light_Scaled_Photo_EE, aes(x = SABLE, y = scaled_TEE, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  ylim(0,10)+
  labs(title="Light cycle (12hrs): TEE",
    y = "TEE light cycle (kcal)")
plot_raw_TEE_light

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_TEE_light,
       filename="Raw_TEE_light_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph: TEE Dark (RAW) ####
Dark_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="off")

plot_raw_TEE_dark <-ggplot(Dark_Scaled_Photo_EE, aes(x = SABLE, y = scaled_TEE, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
              fun.data = mean_se, 
              position = position_dodge(width = 0.8), 
              width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  ylim(0,10)+
  labs(title="Dark cycle (12hrs): TEE",
    y = "TEE dark cycle (kcal)")
plot_raw_TEE_dark

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_TEE_dark,
       filename="Raw_TEE_dark_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph TEE: light/dark at Sable time points (RAW) ####
ggplot(Scaled_Photo_EE, aes(x = lights, y = scaled_TEE, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  labs(title="TEE (12hrs dark/12hrs light)",
    x = "Time point",
    y = "TEE (kcal/photo period)",
    fill = "Treatment group",
    color= "Treatment group") +
    facet_wrap(~SABLE)

#------#
###RMR ####
####Graph: RMR Dark (RAW) ####
Dark_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="off")

plot_raw_RMR_dark <-ggplot(Dark_Scaled_Photo_EE, aes(x = SABLE, y = scaled_RMR, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  #ylim(0,10)+
  labs(title="Dark cycle (12hrs): RMR",
       y = "RMR dark cycle (kcal)")
plot_raw_RMR_dark

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_RMR_dark,
       filename="Raw_RMR_dark_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph: RMR Light (RAW) ####
Light_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="on")

plot_raw_RMR_light <-ggplot(Light_Scaled_Photo_EE, aes(x = SABLE, y = scaled_RMR, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  #ylim(0,10)+
  labs(title="RMR: Light (12hrs)",
       y = "RMR light cycle (kcal)")
plot_raw_RMR_light

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_RMR_light,
       filename="Raw_RMR_light_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph RMR light/dark at Sable time points (RAW) ####
ggplot(Scaled_Photo_EE, aes(x = lights, y = scaled_RMR, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  labs(
    title="RMR (12hrs dark/12hrs light)",
    x = "Time point",
    y = "RMR (kcal/photo period)",
    fill = "Treatment group",
    color= "Treatment group") +
  facet_wrap(~SABLE)

#---#
###NEAT ####
####Graph: NEAT Dark (RAW) ####
Dark_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="off")

plot_raw_NEAT_dark <-ggplot(Dark_Scaled_Photo_EE, aes(x = SABLE, y = scaled_NEAT, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  ylim(0,2)+
  labs(title="Dark cycle (12hrs): NEAT",
       y = "NEAT dark cycle (kcal)")
plot_raw_NEAT_dark

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_NEAT_dark,
       filename="Raw_NEAT_dark_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph: NEAT Light (RAW) ####
Light_Scaled_Photo_EE <- Scaled_Photo_EE %>%
  filter(lights=="on")

plot_raw_NEAT_light <-ggplot(Light_Scaled_Photo_EE, aes(x = SABLE, y = scaled_NEAT, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  ylim(0,2)+
  labs(title="NEAT: Light (12hrs)",
       y = "NEAT light cycle (kcal)")
plot_raw_NEAT_light

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_NEAT_light,
       filename="Raw_NEAT_light_plot.png", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

####Graph NEAT light/dark at Sable time points (RAW) ####
ggplot(Scaled_Photo_EE, aes(x = lights, y = scaled_NEAT, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank(),
        axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  format.plot+
  labs(
    title="NEAT (12hrs dark/12hrs light)",
    x = "Time point",
    y = "NEAT (kcal/photo period)",
    fill = "Treatment group",
    color= "Treatment group") +
  facet_wrap(~SABLE)

#Export plot to folder called "APS_figures" 
ggsave(plot_photo_per_NEAT,
       filename="NEAT(12hrs-12hrs).png", 
       width = 4, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#---#
##Stats: Linear mixed model ####
### MLR: TEE (not adj.) ####
#Build multiple linear regression model for TEE not adjusted for BW or lean #
model_TEE_photo <- lmer(scaled_TEE ~ SABLE * GROUP * lights + (1 | ID), data = Scaled_Photo_EE)
summary(model_TEE_photo)

#Calculate estimated marginal means #
emm_TEE_photo <- emmeans(model_TEE_photo, ~ SABLE * GROUP * lights, cov.reduce = mean)
emm_TEE_photo_df <- as.data.frame(emm_TEE_photo)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_photo <- contrast(emm_TEE_photo, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_photo_df <- as.data.frame(contrasts_by_group_TEE_photo)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_photo <- contrast(emm_TEE_photo, method = "pairwise", by = "SABLE")
contrasts_SABLE_TEE_photo_df <- as.data.frame(contrasts_by_SABLE_TEE_photo)


#-------------------#
### MLR: RMR (not adj.) ####
#Build multiple linear regression model for NEAT not adjusted for BW or lean #
model_RMR_photo <- lmer(scaled_RMR ~ SABLE * GROUP * lights + (1 | ID), data = Scaled_Photo_EE)
summary(model_RMR_photo)

#Calculate estimated marginal means #
emm_RMR_photo <- emmeans(model_RMR_photo, ~ SABLE * GROUP * lights, cov.reduce = mean)
emm_RMR_photo_df <- as.data.frame(emm_RMR_photo)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_photo <- contrast(emm_RMR_photo, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_photo_df <- as.data.frame(contrasts_by_group_RMR_photo)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_photo <- contrast(emm_RMR_photo, method = "pairwise", by = "SABLE")
contrasts_SABLE_RMR_photo_df <- as.data.frame(contrasts_by_SABLE_RMR_photo)


#-------------------#
### MLR: NEAT (not adj.) ####
#Build multiple linear regression model for NEAT not adjusted for BW or lean #
model_NEAT_photo <- lmer(scaled_NEAT ~ SABLE * GROUP * lights + (1 | ID), data = Scaled_Photo_EE)
summary(model_NEAT_photo)

#Calculate estimated marginal means #
emm_NEAT_photo <- emmeans(model_NEAT_photo, ~ SABLE * GROUP * lights, cov.reduce = mean)
emm_NEAT_photo_df <- as.data.frame(emm_NEAT_photo)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_photo <- contrast(emm_NEAT_photo, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_photo_df <- as.data.frame(contrasts_by_group_NEAT_photo)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_photo <- contrast(emm_NEAT_photo, method = "pairwise", by = "SABLE")
contrasts_SABLE_NEAT_photo_df <- as.data.frame(contrasts_by_SABLE_NEAT_photo)

##---Graph estimated EE by photo period (emmeans energy expenditure from mixed model) ####
###TEE estimates ####

#### Graph: TEE Dark (estimates) ####
#filter emm_TEE_photo_df to only include rows when lights = "off"
Dark_emm_TEE_photo_df <- emm_TEE_photo_df %>%
  filter(lights=="off")

plot_Dark_emm_TEE_photo <-ggplot(Dark_emm_TEE_photo_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8)+
  labs(
    title = "TEE (12hrs dark) ",
    y = "TEE (kcal/dark cycle)")
plot_Dark_emm_TEE_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Dark_emm_TEE_photo,
       filename="emm_TEE_dark_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: TEE Light (estimates) ####
#filter emm_TEE_photo_df to only include rows when lights = "on"
Light_emm_TEE_photo_df <- emm_TEE_photo_df %>%
  filter(lights=="on")

plot_Light_emm_TEE_photo <-ggplot(Light_emm_TEE_photo_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "TEE (12hrs light)",
    y = "TEE (kcal/light cycle)")
plot_Light_emm_TEE_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Light_emm_TEE_photo,
       filename="emm_TEE_light_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: TEE Control (estimates) ####
#filter emm_TEE_photo_df to only include control mice
Control_emm_TEE_photo_df <- emm_TEE_photo_df %>%
  filter(GROUP=="Control")

plot_Control_emm_TEE_photo <-ggplot(Control_emm_TEE_photo_df, aes(x = SABLE, y = emmean, fill = lights)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "Control: TEE 12hrs light/dark",
    y = "TEE (kcal/photo period)")
plot_Control_emm_TEE_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Control_emm_TEE_photo,
       filename="emm_TEE_photo_Control_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: TEE Weight cycled (estimates) ####
#filter emm_TEE_photo_df to only include weight cycled mice
Cycled_emm_TEE_photo_df <- emm_TEE_photo_df %>%
  filter(GROUP=="Weight cycled")

plot_Cycled_emm_TEE_photo <-ggplot(Cycled_emm_TEE_photo_df, aes(x = SABLE, y = emmean, fill = lights)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "Weight cycled: TEE 12hrs light/dark",
    y = "TEE (kcal/photo period)")
plot_Cycled_emm_TEE_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Cycled_emm_TEE_photo,
       filename="emm_TEE_photo_Cycled_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

###RMR estimates ####

#### Graph: RMR Dark (estimates) ####
#filter emm_RMR_photo_df to only include rows when lights = "off"
Dark_emm_RMR_photo_df <- emm_RMR_photo_df %>%
  filter(lights=="off")

plot_Dark_emm_RMR_photo <-ggplot(Dark_emm_RMR_photo_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8)+
  labs(
    title = "RMR (12hrs dark) ",
    y = "RMR (kcal/dark cycle)")
plot_Dark_emm_RMR_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Dark_emm_RMR_photo,
       filename="emm_RMR_dark_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: RMR Light (estimates) ####
#filter emm_RMR_photo_df to only include rows when lights = "on"
Light_emm_RMR_photo_df <- emm_RMR_photo_df %>%
  filter(lights=="on")

plot_Light_emm_RMR_photo <-ggplot(Light_emm_RMR_photo_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "RMR (12hrs light)",
    y = "RMR (kcal/light cycle)")
plot_Light_emm_RMR_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Light_emm_RMR_photo,
       filename="emm_RMR_light_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: RMR Control (estimates) ####
#filter emm_RMR_photo_df to only include control mice
Control_emm_RMR_photo_df <- emm_RMR_photo_df %>%
  filter(GROUP=="Control")

plot_Control_emm_RMR_photo <-ggplot(Control_emm_RMR_photo_df, aes(x = SABLE, y = emmean, fill = lights)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "Control: RMR 12hrs light/dark",
    y = "RMR (kcal/photo period)")
plot_Control_emm_RMR_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Control_emm_RMR_photo,
       filename="emm_RMR_photo_Control_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#### Graph: RMR Weight cycled (estimates) ####
#filter emm_RMR_photo_df to only include weight cycled mice
Cycled_emm_RMR_photo_df <- emm_RMR_photo_df %>%
  filter(GROUP=="Weight cycled")

plot_Cycled_emm_RMR_photo <-ggplot(Cycled_emm_RMR_photo_df, aes(x = SABLE, y = emmean, fill = lights)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  ylim(0,8) +
  labs(
    title = "Weight cycled: RMR 12hrs light/dark",
    y = "RMR (kcal/photo period)")
plot_Cycled_emm_RMR_photo

#Export plot to folder called "APS_figures" 
ggsave(plot_Cycled_emm_RMR_photo,
       filename="emm_RMR_photo_Cycled_plot.png", 
       width = 6, 
       height = 5, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#-------------------#
#-------------Start ignore-------------#

#Actually, I think that I should use the model that has lights included
#MLR using df with just dark period
# df = Dark_Scaled_Photo_EE

#TEE
#Build multiple linear regression model for NEAT not adjusted for BW or lean #
model_TEE_dark <- lmer(scaled_TEE ~ SABLE * GROUP + (1 | ID), data = Dark_Scaled_Photo_EE)
summary(model_TEE_dark)

#Calculate estimated marginal means #
emm_TEE_dark <- emmeans(model_TEE_dark, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_dark_df <- as.data.frame(emm_TEE_dark)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_dark <- contrast(emm_TEE_dark, method = "pairwise", by = "GROUP")
contrasts_by_group_TEE_dark_df <- as.data.frame(contrasts_by_group_TEE_dark)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_dark <- contrast(emm_TEE_dark, method = "pairwise", by = "SABLE")
contrasts_SABLE_TEE_dark_df <- as.data.frame(contrasts_by_SABLE_TEE_dark)

#-------------End ignore-------------#



