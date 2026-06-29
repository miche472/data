#RQ for light and dark period
#Used NZO_Figure7b-RMR_correctedbyLean (LM, 10-30).R" as a starting point for the script
#Used "Hourly_TEE_LM, 12-2.R" for the code that separated light and dark periods

#Pending: how can I do contrasts that will assess difference between light and dark period at each time point?
  #I would want to do this for weight cycled mice and then control mice

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

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#### RQ: filter sable_dwn to have only values for the parameter RQ. Re-attach other variables ####

RQ_data <- sable_dwn %>% 
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
  filter(grepl("RQ_*", parameter)) %>% 
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
  
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  group_by(ID, SABLE) %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

RQ_data <- RQ_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )
#filter (ID == "3706") %>%
#filter (SABLE == "BW regain")

#### Attach echoMRI_data to sable_TEE_data --> name new df as sable_TEE_adj ####
#echo info is from NZO_Figure7 - TEE_correctedbyLean (rev. LM).R on (LM accessed on 10-16-25)

#Process echoMRI info for NZO mice
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  dplyr::select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join Lean, Fat, and Weight info into TEE dataset
RQ_data_adj <- RQ_data %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat),
    by = c("ID", "SABLE")
  )

#RQ during the light/dark period for each mouse at each time point

#### RQ during the Light and Dark period for each mouse during each Sable time point ####
#For each mouse at each Sable time point, average the minute by minute TEE (kcal_hr)
#which are during the dark period. If this doesn't work I can average the hourly values

RQ_data_lightcycle <- RQ_data_adj %>%
  #mutate(hour = hour(DateTime)) %>%
  group_by(ID, SABLE, lights) %>%
  summarise(
    RQ_lightcycle = mean(value, na.rm = TRUE),
    GROUP       = first(GROUP),
    DRUG        = first(DRUG),
    #complete_days = first(complete_days),
    .groups = "drop"
  ) %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))) %>%
  filter(!ID == 3714)  #After graphing it looks like this was an erroneous measurement during at peak obesity

#####----- Barplot - NZO RQ during light/dark ####

ggplot(RQ_data_lightcycle, aes(x = SABLE, y = RQ_lightcycle, fill = GROUP)) + 
  stat_summary( # mean bars
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  stat_summary( # error bars (mean ± SE)
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(  # individual data points
    aes(color = GROUP), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  labs(title = "Light/dark cycle RQ", 
       y = "Respiratory quotient (RQ)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  format.plot +
  #ylim(0,1)+
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color="black", size=14),
    plot.title = element_text(hjust = 0.5, face = "bold")) +
  #geom_text_repel(aes(label = ID), size = 3, alpha = 0.7) +
  facet_wrap(~lights)

#Barplot -- only weight cycled mice or only control mice
plot_RQ_WC <- RQ_data_lightcycle %>%
  filter(GROUP=="Weight cycled")

ggplot(plot_RQ_WC, aes(x = SABLE, y = RQ_lightcycle, fill = lights)) + 
  stat_summary( # mean bars
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  stat_summary( # error bars (mean ± SE)
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(  # individual data points
    aes(color = lights), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  labs(title = "RQ for weight cycled mice", 
       y = "Respiratory quotient (RQ)", 
       x= "Time point", 
       color = "Lights",
       fill = "Lights") +
  format.plot +
  #ylim(0,1)+
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color="black", size=14),
    plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_text_repel(aes(label = ID), size = 3, alpha = 0.7)

####----- Light period: Linear mixed model, emmeans, contrasts -----####
#Lights ON
# Build linear mixed model for RQ (lights on)
RQ_data_light_on <- RQ_data_lightcycle %>%
  filter(lights=="on")
  
model_RQ_light <- lmer(RQ_lightcycle ~ SABLE * GROUP + (1 | ID), data = RQ_data_light_on)
summary(model_RQ_light)

#Confirm the number of mice in the data frame is correct...15 mice (3714 removed)
n_distinct(RQ_data_light_on$ID) 

# Calculate estimated marginal means (emmeans) 
emm_RQ_light <- emmeans(model_RQ_light, ~ SABLE * GROUP, cov.reduce = mean)
emm_RQ_light_df <- as.data.frame(emm_RQ_light)

# Pairwise contrasts within each GROUP
contrasts_by_group_RQ_light <- contrast(emm_RQ_light, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_RQ_light_df <- as.data.frame(contrasts_by_group_RQ_light)
contrasts_RQ_light_df

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RQ_light <- contrast(emm_RQ_light, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RQ_light_df <- as.data.frame(contrasts_by_SABLE_RQ_light)
contrasts_SABLE_RQ_light_df

####----- Dark period: Linear mixed model, emmeans, contrasts -----####
#Lights OFF
# Build linear mixed model for RQ (lights off)
RQ_data_light_off <- RQ_data_lightcycle %>%
  filter(lights=="off")

model_RQ_dark <- lmer(RQ_lightcycle ~ SABLE * GROUP + (1 | ID), data = RQ_data_light_off)
summary(model_RQ_dark)

#Confirm the number of mice in the data frame is correct...15 mice (3714 removed)
n_distinct(RQ_data_light_off$ID) 

# Calculate estimated marginal means (emmeans) 
emm_RQ_dark <- emmeans(model_RQ_dark, ~ SABLE * GROUP, cov.reduce = mean)
emm_RQ_dark_df <- as.data.frame(emm_RQ_dark)

# Pairwise contrasts within each GROUP
contrasts_by_group_RQ_dark <- contrast(emm_RQ_dark, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_RQ_dark_df <- as.data.frame(contrasts_by_group_RQ_dark)
contrasts_RQ_dark_df

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RQ_dark <- contrast(emm_RQ_dark, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RQ_dark_df <- as.data.frame(contrasts_by_SABLE_RQ_dark)
contrasts_SABLE_RQ_dark_df

####----- Dark vs light period--Weight cycled mice: Linear mixed model, emmeans, contrasts -----####
#Weight cycled mice
# Build linear mixed model for RQ (lights off)
RQ_data_WC <- RQ_data_lightcycle %>%
  filter(GROUP=="Weight cycled")

model_RQ_WC <- lmer(RQ_lightcycle ~ SABLE * lights + (1 | ID), data = RQ_data_WC)
summary(model_RQ_WC)

#Confirm the number of mice in the data frame is correct...8 mice (weight cycled NZO - 3714)
n_distinct(RQ_data_WC$ID) 

# Calculate estimated marginal means (emmeans) 
emm_RQ_WC <- emmeans(model_RQ_WC, ~ SABLE * lights, cov.reduce = mean)
emm_RQ_WC_df <- as.data.frame(emm_RQ_dark)

# Pairwise contrasts within each lights
contrasts_by_group_RQ_WC <- contrast(emm_RQ_WC, method = "pairwise", by = "lights")
# Convert to a data frame
contrasts_RQ_WC_df <- as.data.frame(contrasts_by_group_RQ_WC)
contrasts_RQ_WC_df

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_RQ_WC <- contrast(emm_RQ_WC, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_RQ_WC_df <- as.data.frame(contrasts_by_SABLE_RQ_WC)
contrasts_SABLE_RQ_WC_df
