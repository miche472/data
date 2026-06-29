
#Took components of 10th percentile method for RMR (Rev. 4-3-26) and modified
#the aesthetics of the graphs (RAW and emmeans) to be suitable for APS poster

#Started: 4-3-26
#Revised: 4-3-26

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

#---
#------------Create dfs necessary to make graphs------------ ####


# Code for df = filter_locom_energy_3 

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

## Calculate sum of EE measured each minute ####
# RMR calculated using percentile Percentile for RMR  
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

##Calculate daily sum (add all hours together) ####
Daily_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day))) 

#-------------------------------------------------#
# Graphs ####
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
# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
custom_colors2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")

#------------Graphs: Raw daily EE------------####
## Graph: Raw daily TEE ####
#update: angle text labels or increase width between bars or decrease font size
#set y axis limit; remove gridlines

plot_raw_TEE_daily <- ggplot(Daily_EE, aes(x = SABLE, y = TEE_kcal_day, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(#legend.position = "right",
    legend.position = "none",
    #plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    #legend.title = element_text(size = 12, face="bold"),
    #legend.text = element_text(size = 12),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20, face="bold"),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  labs(
    #title="Total energy expenditure (kcal/day)",
    y = "Total energy expenditure (kcal/day)")
plot_raw_TEE_daily

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_TEE_daily,
       filename="Raw_TEE_daily_plot.png", 
       width = 6.8, 
       height = 6.8, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#---
## Graph: Raw daily RMR ####
plot_raw_RMR_daily <- ggplot(Daily_EE, aes(x = SABLE, y = RMR_kcal_day, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(#legend.position = "right",
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    #legend.title = element_text(size = 12, face="bold"),
    #legend.text = element_text(size = 12),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20, face="bold"),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  labs(
    #title="Resting metabolic rate (kcal/day)",
    y = "Resting metabolic rate (kcal/day)")
plot_raw_RMR_daily

#Export plot to folder called "APS_figures" 
ggsave(plot_raw_RMR_daily,
       filename="Raw_RMR_daily_plot.png", 
       width = 6.8, 
       height = 6.8, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#-------------------------------------------------#
# Graphs: emmeans daily EE ####
#adjusted only for multiple measurements from same mouse and SABLE*GROUP interaction

## Emmeans daily TEE ####
### Model (TEE)####
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
# Convert to a data frame
contrasts_SABLE_TEE_10_df <- as.data.frame(contrasts_by_SABLE_TEE_10)

###Graph (TEE) ####
# Bar plot - Graph predicted TEE unadjusted for lean or BW #
plot_emm_TEE_daily <-ggplot(emm_TEE_10_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "none",
    #plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 22, face="bold"),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  scale_y_continuous(breaks = seq(0, 14, by = 3),
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    #title = "Total energy expenditure (kcal/day)",
    y = "Total energy expenditure (kcal/day)")
plot_emm_TEE_daily

#Export plot to folder called "APS_figures" 
ggsave(plot_emm_TEE_daily,
       filename="emm_TEE_daily_plot.png", 
       width = 7, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#----------------------#
## Emmeans daily RMR ####
### Model (RMR)####
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

###Graph (RMR) ####
# Bar plot - Graph predicted RMR unadjusted for lean or BW #
plot_emm_RMR_daily <-ggplot(emm_RMR_10_df, aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  # error bars using SE
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
    axis.title.y = element_text(size= 22, face="bold"),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  scale_y_continuous(breaks = seq(0, 14, by = 3),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    #title = "Resting metabolic rate (kcal/day)",
    y = "Resting metabolic rate (kcal/day)")
plot_emm_RMR_daily

#Export plot to folder called "APS_figures" 
ggsave(plot_emm_RMR_daily,
       filename="emm_RMR_daily_plot.png", 
       width = 7, 
       height = 7, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")