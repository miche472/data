
#Started: 3-3-26
#Revised: 3-3-26

#Data compilation for SSIB 2026 abstract and APS Summit 2026
#Author: LM
#Goal: Determine if weight cycling corresponds to a decrease in the proportion of TEE contributed by RMR

#Conclusion based on Objectives 2 and 3 (MLR for RMR/TEE, NEAT/TEE, and RMR/TEE adjusted for lean mass):
  #Objective 2:Proportion of TEE comprised of RMR and NEAT isn't compelling. 
  #Objective 3:Proportion of TEE comprised of RMR and NEAT isn't compelling even 
    #after adjusting for lean mass

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

#Get complete day 1 and 2 by explicitly definining the zt_time and SABLE_DAY
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

#---
#1) Calculate Energy expenditure: TEE, NEAT, RMR ####
#df with minute by minute EE and AllMeters is "filter_locom_energy3

#Total energy expenditure (total kcal/day --> minute summation method)
ID_TEE3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, LM_complete_day) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(TEE_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(TEE_teske = sum(TEE_per_min), .groups="drop")

#NEAT (kcal/day --> minute summation method)
ID_NEAT3 <- filter_locom_energy3 %>%
  ungroup() %>%
  group_by(SABLE, ID, GROUP, LM_complete_day) %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>% #only keep data from minutes when the mouse moved
  mutate(NEAT_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(NEAT_teske = sum(NEAT_per_min), .groups="drop")

#Resting metabolic rate (kcal/day --> minute summation method)
ID_RMR3 <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==0) %>% #only keep data from minutes when the mouse didn't move
  mutate(RMR_per_min= (Kcal_Hr/60)) %>% #interpolate EE per minute from EE per hr
  summarise(RMR_teske = sum(RMR_per_min), .groups="drop")

#---
#2a) Calculate proportion of TEE comprised of RMR and NEAT for each mouse at each time point ####
TEE_composition <- ID_RMR3 %>%
  left_join(
    ID_NEAT3 %>% 
      select(NEAT_teske, ID, GROUP, SABLE), by = c("ID", "GROUP", "SABLE"))  %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  mutate(NEAT_plus_RMR = NEAT_teske + RMR_teske) %>% #check --> good
  left_join(
    ID_TEE3 %>% 
      select(TEE_teske, ID, GROUP, SABLE), by = c("ID", "GROUP", "SABLE")) %>%
  mutate(diff_TEE = TEE_teske - NEAT_plus_RMR, #check --> good
         RMR_TEE = RMR_teske/TEE_teske, #proportion of TEE that is RMR
         NEAT_TEE = NEAT_teske/TEE_teske, #proportion of TEE that is NEAT
         Whole = RMR_TEE + NEAT_TEE) #check --> good

#---
#2b) Conduct MLR for RMR/TEE --> does this change during weight cycling? ####
  #Use df called TEE_composition

#Build multiple linear regression model for RMR/TEE (raw values)
model_RMR_TEE <- lmer(RMR_TEE ~ SABLE * GROUP + (1 | ID), data = TEE_composition)
summary(model_RMR_TEE)

# Calculate estimated marginal means 
emm_RMR_TEE <- emmeans(model_RMR_TEE, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_TEE_df <- as.data.frame(emm_RMR_TEE)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_TEE <- contrast(emm_RMR_TEE, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_TEE_df <- as.data.frame(contrasts_by_group_RMR_TEE)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_TEE <- contrast(emm_RMR_TEE, method = "pairwise", by = "SABLE")
contrasts_SABLE_RMR_TEE_df <- as.data.frame(contrasts_by_SABLE_RMR_TEE)

#Conclusions:RMR/TEE ####
  #Within GROUP contrasts: at BW maintenance, weight cycled mice have a 22.9% higher proportion of TEE 
    #comprised of RMR compared to control mice (p=0.0028)
  #Within time point contrasts: 8 comparisons that are sig. different (3 cycled and 5 control), 
    #but not notable comparisons

#Bar plot - Graph predicted RMR/TEE 
barplot_emm_RMR_TEE <- emm_RMR_TEE_df %>%
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
    title = "Proportion of TEE comprised of RMR",
    y = "RMR/TEE",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_TEE

#---
#2c) Condult MLR for RMR/NEAT --> does this change during weight cycling? ####
  #Use df called TEE_composition

#Build multiple linear regression model for NEAT/TEE (raw values)
model_NEAT_TEE <- lmer(NEAT_TEE ~ SABLE * GROUP + (1 | ID), data = TEE_composition)
summary(model_NEAT_TEE)

# Calculate estimated marginal means 
emm_NEAT_TEE <- emmeans(model_NEAT_TEE, ~ SABLE * GROUP, cov.reduce = mean)
emm_NEAT_TEE_df <- as.data.frame(emm_NEAT_TEE)

# Pairwise contrasts within each GROUP
contrasts_by_group_NEAT_TEE <- contrast(emm_NEAT_TEE, method = "pairwise", by = "GROUP")
contrasts_by_group_NEAT_TEE_df <- as.data.frame(contrasts_by_group_NEAT_TEE)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_NEAT_TEE <- contrast(emm_NEAT_TEE, method = "pairwise", by = "SABLE")
contrasts_SABLE_NEAT_TEE_df <- as.data.frame(contrasts_by_SABLE_NEAT_TEE)

#Conclusions:NEAT/TEE ####
  #Within GROUP contrasts: at BW maintenance, weight cycled mice have a 22.9% lower proportion of TEE 
    #comprised of NEAT compared to control mice (p=0.0028)
  #Within time point contrasts: 8 comparisons that are sig. different (3 cycled and 5 control), 
    #but not notable comparisons

#3) Independent of lean mass, does composition of TEE change with weight cycling? ####
    # Include lean mass in the MLR models. 

#Prepare echo data for merging with TEE data
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
  filter(!(ID == 3711 & Date == as.Date("2025-01-27"))) %>% #repeated
  filter(!(ID == 3727 & Date == as.Date("2025-02-07"))) %>% #repeated
  mutate(STATUS = factor(STATUS,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>% 
  rename(SABLE = STATUS) 

#Check number of mice in echoMRI_data
echoMRI_data %>% group_by(SABLE) %>% summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#Combine echoMRI_data with TEE_composition
Echo_TEE_composition <- TEE_composition   %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))

# Multiple linear regression for RMR/TEE, including Lean as a covariate
#Build multiple linear regression model for RMR/TEE (raw values)
model_RMR_TEE_lean <- lmer(RMR_TEE ~ SABLE * GROUP + Lean + (1 | ID), data = Echo_TEE_composition)
summary(model_RMR_TEE_lean)

# Calculate estimated marginal means 
emm_RMR_TEE_lean <- emmeans(model_RMR_TEE_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_TEE_lean_df <- as.data.frame(emm_RMR_TEE_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_RMR_TEE_lean <- contrast(emm_RMR_TEE_lean, method = "pairwise", by = "GROUP")
contrasts_by_group_RMR_TEE_lean_df <- as.data.frame(contrasts_by_group_RMR_TEE_lean)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RMR_TEE_lean <- contrast(emm_RMR_TEE_lean, method = "pairwise", by = "SABLE")
contrasts_SABLE_RMR_TEE_lean_df <- as.data.frame(contrasts_by_SABLE_RMR_TEE_lean)

#Conclusions:RMR/TEE, adjusted for lean mass ####
#Within GROUP contrasts: at BW maintenance, weight cycled mice have a 23.5% lower proportion of TEE 
  #comprised of RMR compared to control mice (p=0.006)
#Within time point contrasts: 5 comparisons that are sig. different (1 cycled and 4 control).
  #RMR contributed significantly less to TEE at BW regain compared to BW loss. Difference is 21.9% (p= 0.038).
  #None of the comparisons were conceptually notable.

#Bar plot - Graph predicted RMR/TEE, adjusted for lean mass
barplot_emm_RMR_TEE_lean <- emm_RMR_TEE_lean_df %>%
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
    title = "Proportion of TEE comprised of RMR, adjusted for lean mass",
    y = "RMR/TEE adjusted for lean mass",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_RMR_TEE_lean

#-----------------#
#For graphing
#----------------Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

