
#Major topics: % body weght change, ΔBW over time, Energy efficiency (ΔBW/kcal) ####
    #Date started: 2-17-26
    #Revised:3-20-26 
#Revision on 2-23-26 changed cumulative FI to use INTAKE_GR rather than corrected_intake_kcal
#Revision on 3-10-26 --> changed labeling of x axis to "obesity development, Weight loss, 
#weight maintenance, and weight regain". The format of "Baseline->peak obesity etc. is in revision from 2-23-26

#This script is a version of "Energy efficiency & %BW (Rev. 3-20-26).R" that has
#statistical analyses for change in BW over stages and average daily FI during stages;
#there are also graphs associated with the linear mixed model statistical analyses

#Objectives:
    #Objective 1: Percent change in BW between each BW measurement, relative to first BW measurement
    #Objective 2: All measurements of BW and FI over the course of the study 
    #Objective 3: Energy efficiency during each transition period
    #Objective 4: Consider ΔBW (g) during each transition period (independent of FI)
    #Objective 5: Calculate percent change in BW
    #Objective 6: calculate duration of each stage and FI during each stage (cumulative FI/days in stage)

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(patchwork)

#-------------------------------------------------------#
#Objective 1: Percent change in BW over time ####

#Plot percent change over time for each mouse ID (include all measurements, not just sign posts) ####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
rename(COMMENTS_BW = COMMENTS) %>%
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"),
    day_rel = DATE - first(DATE)) %>%
  mutate(
    STATE = case_when(
      ##Baseline: First day of LFD/first day of obesity development
      #Date is the first date after start of LFD for which there was a BW and FI measurement
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
      ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
      ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
      
      ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
      #First day of calorie restriction for Weight Cycled mice. 
      #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
      ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
      
      ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
      #Sable day 15
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
      
      ##BW maintenance: End of Sable recording for BW maintenance
      #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
      ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
      ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
      
      ##BW regain: Final day of experiment --> day of sac (End of regain)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!(ID == 3726 & DATE == as.Date("2025-04-28"))) %>%  #repeated
  mutate(STATE = factor(STATE, # Make STATE an ordered factor
      levels = c("baseline", "peak obesity", "BW loss", "BW maintenance", "BW regain")))

#Format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        #   strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

#Plot percent change in BW (relative to starting BW...BW on Sable day 7) at each subsequent BW measurement
plot_BW_percent <- BW_data %>%
  ggplot(aes(x = day_rel, y = bw_rel, group = ID, fill=GROUP)) +
  # individual trajectories
  geom_line(alpha = 0.3) +   
  geom_point(size = 2, alpha = 0.3) +  
  # mean ± SD ribbon (need to set 'fill' separately from 'color')
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1), 
    geom = "ribbon", aes(group = GROUP, fill = GROUP), 
    alpha = 0.2, color = NA) +
  # mean solid line
  stat_summary(
    fun = mean, geom = "line", aes(group = GROUP, color = GROUP), 
    linewidth = 1.2) +
  # mean dashed line (optional, if you want to keep it too)
  # stat_summary(fun = mean, geom = "line", aes(group = DRUG, color = DRUG), 
  #              size = 1.2, linetype = "dashed") +
  theme_minimal() +
  labs(y = "Body weight relative to baseline (%)", color = "Treatment group", fill = "Treatment group") +
  #facet_wrap(~ID) + #To see each mouse separately
  #facet_wrap(~GROUP) + #To see Control and Weight cycled mice separately
  format.plot+
  labs(title="Percent change in body weight",
       color="Treatment group",
       fill="Treatment group",
       x="Days",
       y="Body weight relative to baseline (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_BW_percent

#-----------------------------------------------------------------------------#
#Objective 2: All measurements of BW and FI over the course of the study ####

#Prepare a df to look at every measurement of FI and BW (df_all_BW_FI)
#The data frame called df_all_BW_FI in this script is based on a df called df1 in "FI & BW (started 2-16-26).R"
  #To create df_all_BW_FI, first run the chunk of code below which combines FI.csv and BW.csv to create FI_BW_joined
  #df_all_BW_FI retains measurements from all days whereas df1 only retains measurements from sign post days
df_all_BW_FI <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW)) %>% 
  mutate(day_rel = DATE - first(DATE),
         FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
         STATE = case_when(
        ##Baseline: First day of LFD/first day of obesity development
           #Date is the first date after start of LFD for which there was a BW and FI measurement
           ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
           ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
           ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
           ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
           
        ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
           #First day of calorie restriction for Weight Cycled mice. 
           #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
           ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                                             3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                                             3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
           
        ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
           #Sable day 15
           ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
           ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           
        ##BW maintenance: End of Sable recording for BW maintenance
           #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           
        ##BW regain: Final day of experiment --> day of sac (End of regain)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  #in this script we want to retain data from all dates, not just phase "sign post" dates
  mutate(STATE = factor(STATE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 

#-------------------------------------------------------#
#Objective 3: Energy efficiency during each transition period ####
#Calculate cumulative FI and BW change during each transition period -->

#Create df1 directly in this script (originally created in "FI & BW (Started 2-16-26).R")
    #Prepare BW.csv
    BW_to_join <- read_csv("../data/BW.csv") %>% 
      filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
      filter(!ID %in% c(3712, 3715)) %>% # died during study
      group_by(ID) %>% 
      arrange(DATE) %>%
      rename(COMMENTS_BW = COMMENTS)

    #Prepare FI.csv
    FI_to_join <- read_csv("../data/FI.csv") %>% 
      filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
      filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% #remove time when fed chow
      filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
      filter(!ID %in% c(3712, 3715)) %>% # died during study
      group_by(ID) %>% 
      arrange(DATE) %>% 
      rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
      select(-DIET_FORMULA.y) %>% 
      filter(!is.na(corrected_intake_gr)) %>% 
      mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>%
      rename(COMMENTS_FI = COMMENTS)
    
#Create FI_BW_joined
    #Join prepared BW and FI data frames
    FI_BW_joined <- FI_to_join %>%
      left_join(
        BW_to_join %>% 
          select(ID, DATE, BW, COMMENTS_BW),
        by = c("ID", "DATE"))

# Create df1 
  #Adds variables: GROUP, DRUG, STATE, day_rel, FI_rel, FI_cum to joined BW & FI data
df1 <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(day_rel = DATE - first(DATE),
         FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
         STATE = case_when(
        ##Baseline: First day of LFD/first day of obesity development
           #Date is the first date after start of LFD for which there was a BW and FI measurement
           ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
           ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
           ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
           ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
           
        ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
           #First day of calorie restriction for Weight Cycled mice. 
           #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
           ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                     3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                     3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
           
        ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
           #Sable day 15
           ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
           ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           
        ##BW maintenance: End of Sable recording for BW maintenance
           #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           
        ##BW regain: Final day of experiment --> day of sac (End of regain)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(STATE)) %>%   #Remove measurements that aren't from sign post dates (i.e. STATE)
  mutate(STATE = factor(STATE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 

#FI by transition period: Use df1 to calculate cumulative FI for each transition period ####
# Summarize cumulative FI per ID and STATE
    FI_stage_summary <- df1 %>%
      group_by(ID, GROUP, DRUG, STATE) %>%
      summarise(FI_cum_end = max(FI_cum_INTAKE_kcal, na.rm = TRUE), .groups = "drop") %>%
      # Reshape into wide format: one row per ID, columns = each STATUS
      pivot_wider(
        names_from = STATE,
        values_from = FI_cum_end) %>%
      # Calculate kcal consumed between stages
      mutate(
        kcal_baseline_to_peak = `Peak obesity` - Baseline,
        kcal_peak_to_loss = `BW loss` - `Peak obesity`,
        kcal_loss_to_maint = `BW maintenance` - `BW loss`,
        kcal_maint_to_regain = `BW regain` - `BW maintenance`)
    
    # Convert to long format
    FI_stage_long <- FI_stage_summary %>%
      select(ID, DRUG, GROUP,
             kcal_baseline_to_peak,
             kcal_peak_to_loss,
             kcal_loss_to_maint,
             kcal_maint_to_regain) %>%
      pivot_longer(
        cols = starts_with("kcal_"),
        names_to = "Transition",
        values_to = "kcal") %>%
      mutate(
        Transition = factor(
          Transition,
          levels = c("kcal_baseline_to_peak",
                     "kcal_peak_to_loss",
                     "kcal_loss_to_maint",
                     "kcal_maint_to_regain"),
          labels = c("Obesity development",         #formerly Baseline → Peak obesity
                     "Weight loss",                 #formerly Peak obesity → BW loss
                     "Weight maintenance",          #formerly BW loss → BW maintenance
                     "Weight regain")))             #formerly BW maintenance → BW regain

#Duration of each transition: use df1 to calculate duration (days) for each transition period ####
# Summarize duration (days) per ID and STATE
    Days_stage_summary <- df1 %>%
      group_by(ID, GROUP, DRUG, STATE) %>%
      summarise(max_day_rel = max(day_rel, na.rm = TRUE), .groups = "drop") %>%
      # Reshape into wide format: one row per ID, columns = each STATUS
      pivot_wider(
        names_from = STATE,
        values_from = max_day_rel) %>%
      # Calculate days between stages
      mutate(
        Days_baseline_to_peak = `Peak obesity` - Baseline,
        Days_peak_to_loss = `BW loss` - `Peak obesity`,
        Days_loss_to_maint = `BW maintenance` - `BW loss`,
        Days_maint_to_regain = `BW regain` - `BW maintenance`)
    
    # Convert to long format
    Days_stage_long <- Days_stage_summary %>%
      select(ID, DRUG, GROUP,
             Days_baseline_to_peak,
             Days_peak_to_loss,
             Days_loss_to_maint,
             Days_maint_to_regain) %>%
      pivot_longer(
        cols = starts_with("Days_"),
        names_to = "Transition",
        values_to = "Duration_day") %>%
      mutate(
        Transition = factor(
          Transition,
          levels = c("Days_baseline_to_peak",
                     "Days_peak_to_loss",
                     "Days_loss_to_maint",
                     "Days_maint_to_regain"),
          labels = c("Obesity development",    #formerly Baseline → Peak obesity
                     "Weight loss",            #formerly Peak obesity → BW loss
                     "Weight maintenance",     #formerly BW loss → BW maintenance
                     "Weight regain")))        #formerly BW maintenance → BW regain
    
# df2: Join duration (days) and cummulative FI during each phase
  FI_duration_joined <- FI_stage_long %>%
      left_join(
        Days_stage_long %>% 
          select(Transition, ID, Duration_day),
        by = c("ID", "Transition"))
  
  df2 <- FI_duration_joined %>%
      group_by(GROUP, Transition, ID) %>%
      mutate(Daily_kcal = kcal/as.numeric(Duration_day))
    
    #Verify that there are 22 mice in each of the 4 transition periods 
    df2 %>% 
      group_by(Transition) %>%
      summarise(n_ID = n_distinct(ID)) #this we have 22 NZO in each transition period

#ΔBW during each transition: use df1 to calculate delta_BW (g) for each transition period ####
# Summarize Δ BW (g) per ID and STATE
  BW_stage_summary <- df1 %>% #df1 has 5 values for each mouse (one per time point)
  group_by(ID, GROUP, STATE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,         #time period (5 total)
    values_from = BW_end) %>%   #BW at end of each time period
  # Calculate change in BW between time points for each ID
  mutate(
    BWg_baseline_to_peak = `Peak obesity` - Baseline,
    BWg_peak_to_loss = `BW loss` - `Peak obesity`,
    BWg_loss_to_maint = `BW maintenance` - `BW loss`,
    BWg_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
BW_stage_long <- BW_stage_summary %>%
  select(ID, GROUP,
         BWg_baseline_to_peak,
         BWg_peak_to_loss,
         BWg_loss_to_maint,
         BWg_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("BWg_"),
    names_to = "Transition",       #name of new column which lists all transition periods for all mice
    values_to = "delta_BW_g") %>%  #name of new column with delta BW values for each transition period
  mutate(
    Transition = factor(
      Transition,
      levels = c("BWg_baseline_to_peak",
                 "BWg_peak_to_loss",
                 "BWg_loss_to_maint",
                 "BWg_maint_to_regain"),
      labels = c("Obesity development",    #formerly Baseline → Peak obesity
                 "Weight loss",            #formerly Peak obesity → BW loss
                 "Weight maintenance",     #formerly BW loss → BW maintenance
                 "Weight regain")))        #formerly BW maintenance → BW regain

#Combine cumulative FI, duration, ΔBW for each phase -> Add BW_stage_long to df2 
    #df2 has cumulative FI and # of days in each transition period
FI_duration_BW_joined <- df2 %>%
  left_join(
    BW_stage_long %>% 
      select(Transition, ID, delta_BW_g),
    by = c("ID", "Transition"))

# Calculate change in BW per kcal consumed (phase_BW_FI) for each ID during each Transition
df4 <- FI_duration_BW_joined %>%
  group_by(GROUP, Transition, ID) %>%
  mutate(phase_BW_FI = delta_BW_g/kcal)

#For df4: Verify that there are 22 mice in each of the 4 transition periods 
    df4 %>% 
      group_by(Transition) %>%
      summarise(n_ID = n_distinct(ID)) #we have 22 NZO in each transition period

# Define custom colors
    custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
    custom_colors2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")
#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove minor grid lines only
  panel.grid.major = element_blank(), # remove major grid lines only
  axis.line = element_line(color = "black")) # keep axis lines
    
## (use) Graph Energy Efficiency for Control and Weight cycled mice during each Transition ####
    #mean change in BW (g) per kcal consumed for each GROUP during every Transition
ggplot(df4, aes(x = Transition, y = phase_BW_FI, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
        fun = "mean", 
        position = position_dodge(width = 0.8), width=0.73) +
      geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
                 alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
        fun.data = mean_se, 
        position = position_dodge(width = 0.9), 
        width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  format.plot+
  theme(legend.position = "right",
            plot.title = element_text(size=16, hjust = 0.5, face="bold"),
            legend.title = element_text(size = 12, face="bold"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
            axis.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 12)) +
  labs(
    title="Energy efficiency (ΔBW/cumulative intake)",
    x = "Stage of weight cycle",
    y = "ΔBW (g) per kcal consumed",
    fill = "Treatment group",
    color = "Treatment group")
    
#(use!) Graph Energy Efficiency for ONLY Obesity development & Weight regain ####
df4_efficient <- df4 %>%
  filter(Transition %in% c("Obesity development", "Weight regain"))
    
ggplot(df4_efficient, aes(x = Transition, y = phase_BW_FI, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), 
           width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
            alpha = 0.6, 
            size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8),
                width = 0.3, 
                color="black") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  #theme_minimal() + #removes the lines that make the graph a box
  format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  labs(title="Figure 2C: Energy efficiency",
    #x = "Stage of weight cycle",
    y = "ΔBW (g) per kcal consumed",
    fill = "Treatment group",
    color = "Treatment group")

# START ignore #----------
#POSSIBLE ALTERNATIVE/reverse method for energy efficiency: kcals of intake required to cause a 1g change in BW
alt_efficiency <- FI_duration_BW_joined %>%
      group_by(GROUP, Transition, ID) %>%
      mutate(phase_FI_BW = kcal/delta_BW_g)

#Graph for alternative method for 
ggplot(alt_efficiency, aes(x = Transition, y = phase_FI_BW, fill = GROUP)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary",fun.data = mean_se,position = position_dodge(width = 0.9),width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  format.plot+
  theme(legend.position = "right",
            plot.title = element_text(size=16, hjust = 0.5, face="bold"),
            legend.title = element_text(size = 12, face="bold"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
            axis.title = element_text(face = "bold"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 12)) +
      labs(title="Alt. energy efficiency (Cumulative intake/ΔBW)",
        x = "Stage of weight cycle",
        y = "Cumulative intake (kcal)/ΔBW (g)",
        fill = "Treatment group")   
    
# END ignore #----------
    
    
# Objective 4: Consider ΔBW (g) during each transition period (independent of FI) ####
#df4 has delta_BW_grams which is the change in BW for each ID during each transition period
## (use!) Graph ΔBW (g) for Control and Weight cycled mice--> may help explain how weight cycled mice at more but lost weight ####
ggplot(df4, aes(x = Transition, y = delta_BW_g, fill = GROUP, color=GROUP)) +
      geom_bar(stat = "summary", 
               fun = "mean", 
               position = position_dodge(width = 0.8), 
               width=0.73) +
      geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 2) +
      geom_errorbar(stat = "summary", 
                    fun.data = mean_se, 
                    position = position_dodge(width = 0.8), 
                    width = 0.3, color="black") +
      theme_bw(base_size = 14) +
      scale_fill_manual(values = custom_colors) +
      scale_color_manual(values = custom_colors2) +
     # theme_minimal() + #removes the lines that make the graph a box
      format.plot+
      theme(legend.position = "right",
            plot.title = element_text(size=16, hjust = 0.5, face="bold"),
            legend.title = element_text(size = 12, face="bold"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
            strip.text = element_text(face = "bold", size = 12),
            axis.title.x = element_blank()) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
      labs(
        title="Figure 2A: Change in body weight (g)",
        #x = "Stage of weight cycle",
        y = "Change in body weight (g)",
        fill = "Treatment group",
        color= "Treatment group")
#This shows change in BW in grams, but percent change in BW may be more useful since there was a wide range in BW
    
#Objective 5: Calculate percent change in BW ####
  #Need to group by ID and state and then do--> bw_rel = 100 * (BW - first(BW)) / first(BW)
  #Then filter out measurements that don't have a STATE?...would this give me % change between start and end of a transition period?
#Two methods for doing this are shown below
    
#Percent change relative to the "start" of each STATE rather than to Baseline BW
  #This method works, but using the lag(BW) approach is risky
    BW_try <- FI_BW_joined %>%
      ungroup() %>%
      group_by(ID) %>% 
      arrange(DATE) %>% 
      mutate(
        GROUP = case_when(
          ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
          ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
        DRUG = case_when(
          ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
          ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
      mutate(day_rel = DATE - first(DATE),
            STATE = case_when(
            ##Baseline: First day of LFD/first day of obesity development
               #Date is the first date after start of LFD for which there was a BW and FI measurement
               ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
               ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
               ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
               ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
               
            ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
               #First day of calorie restriction for Weight Cycled mice. 
               #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
               ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                         3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
               ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                         3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
               ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
               
            ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
               #Sable day 15
               ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
               ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
               ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
               ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
               ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
               
            ##BW maintenance: End of Sable recording for BW maintenance
               #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
               ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
               ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
               ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
               ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
               ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
               ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
               ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
               ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
               
            ##BW regain: Final day of experiment --> day of sac (End of regain)
               ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
               ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
               ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
               ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
               ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
               ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
               ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
               TRUE ~ NA_character_)) %>%
      filter(!is.na(STATE)) %>%   #Retain all measurements
      mutate(STATE = factor(STATE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>%
      ungroup() %>%
      group_by(ID) %>%
      arrange(DATE) %>%
    mutate(bw_rel = 100 * (BW - lag(BW)) / lag(BW)) 
    
#Method 2: for percent change in BW
    #Instead of using lag(BW) try modifying the code that creates BW_stage_summary.
    #Modification: calculate % change rather than change in grams between phases
    #ΔBW during each transition: use df1 to calculate delta_BW (g) for each transition period
    # Summarize Δ BW (g) per ID and STATE
    BW_stage_summary_percent <- df1 %>% #df1 has 5 values for each mouse (one per time point)
      group_by(ID, GROUP, STATE) %>%
      summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
      # Reshape into wide format: one row per ID, columns = each STATUS
      pivot_wider(
        names_from = STATE,         #time period (5 total)
        values_from = BW_end) %>%   #BW at end of each time period
      # Calculate change in BW between time points for each ID
      group_by(ID) %>%
      mutate(
        BWp_baseline_to_peak = 100*((`Peak obesity` - Baseline)/Baseline), 
        BWp_peak_to_loss = 100*((`BW loss` - `Peak obesity`)/`Peak obesity`),
        BWp_loss_to_maint = 100*((`BW maintenance` - `BW loss`)/`BW loss`),
        BWp_maint_to_regain = 100*((`BW regain` - `BW maintenance`)/`BW maintenance`))
    
# Convert to long format
    BW_stage_long_percent <- BW_stage_summary_percent %>%
      select(ID, GROUP,
             BWp_baseline_to_peak,
             BWp_peak_to_loss,
             BWp_loss_to_maint,
             BWp_maint_to_regain) %>%
      pivot_longer(
        cols = starts_with("BWp_"),
        names_to = "Transition",       #name of new column which lists all transition periods for all mice
        values_to = "delta_BW_perce") %>%  #name of new column with delta BW values for each transition period
      mutate(
        Transition = factor(
          Transition,
          levels = c("BWp_baseline_to_peak",
                     "BWp_peak_to_loss",
                     "BWp_loss_to_maint",
                     "BWp_maint_to_regain"),
          labels = c("Obesity development",    #formerly Baseline → Peak obesity
                     "Weight loss",            #formerly Peak obesity → BW loss
                     "Weight maintenance",     #formerly BW loss → BW maintenance
                     "Weight regain")))        #formerly BW maintenance → BW regain
    
# df5: Add BW_stage_long to df4. df4 has cumulative FI, # of days, and ΔBW (g) in each transition period
    
df5 <- df4 %>%
      left_join(
        BW_stage_long_percent %>% 
          select(Transition, ID, delta_BW_perce),
        by = c("ID", "Transition"))

## (use!) Graph percent change in BW for Control and Weight cycled mice ####
ggplot(df5, aes(x = Transition, y = delta_BW_perce, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
        #axis.title = element_text(face = "bold"),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  labs(
    title="Figure 2A: Change in body weight (%)",
    #x = "Stage of weight cycle", 
    y = "Change in body weight (%)",
    fill = "Treatment group", color="Treatment group")

##Stat analysis: %BW change by stage --> located in script called 
              #"Change in BW during weight cycle (start 3-19-26).R"

#------------------------------------------#
#Objective 6: FI during stages of weight cycle ####
#Is there more fluctuation in daily FI for ad libitum mice compared to during BW loss? ###
Explore_flux <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  group_by(ID, GROUP) %>%
  arrange(DATE) %>%
  filter(DATE>= as.Date("2025-02-24") & DATE<=as.Date("2025-06-27"))
  
#Graph daily FI over the course of BW loss (look at daily fluctuation in each group)
  ggplot(Explore_flux, aes(x = DATE, y = INTAKE_GR, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "right",
          plot.title = element_text(size=16, hjust = 0.5, face="bold"),
          legend.title = element_text(size = 12, face="bold"),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size= 13, angle = 45, hjust = 1, color="black"),
          axis.title = element_text(face = "bold")) +
  labs(title="Daily food intake (kcal)",
    x = "Date",
    y = "Food intake (kcal/day)",
    fill = "Treatment group") +
  facet_wrap(~GROUP) +
  #facet_wrap(~ID)
#--------------------------------- #
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
  custom_colors <- c(
    "Control" = "#FAAC41",              
    "Weight cycled" = "#3498DB")
  
## (use!) Graph FI (raw values) -> during the stages ####
ggplot(df4, aes(x = Transition, y = Daily_kcal, fill = GROUP, color=GROUP)) +
    geom_bar(stat = "summary", 
             fun = "mean", 
             position = position_dodge(width = 0.8),
             width=0.73) +
    geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
               alpha = 0.6, 
               size = 2) +
    geom_errorbar(stat = "summary", 
                  fun.data = mean_se, 
                  position = position_dodge(width = 0.8), 
                  width = 0.3,
                  color="black") +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors2) +
    #theme_minimal() + #removes the lines that make the graph a box
    format.plot+
    labs(title="Figure 2B: Energy intake (kcal/day)",
      y = "Energy intake (kcal/day)",
      fill = "Treatment group",
      color = "Treatment group") +
    theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank())
  
## (use) Stat analysis: FI during stages (linear mixed model for daily FI) ####
  #Build multiple linear regression model
  model_FI_daily_stage <- lmer(Daily_kcal ~ Transition * GROUP + (1 | ID), data = df4)
  summary(model_FI_daily_stage)
  
  # Calculate estimated marginal means 
  emm_FI_daily_stage <- emmeans(model_FI_daily_stage, ~ Transition * GROUP, cov.reduce = mean)
  emm_FI_daily_stage_df <- as.data.frame(emm_FI_daily_stage)
  
  # Pairwise contrasts within each GROUP
  contrasts_by_group_FI_daily_stage <- contrast(emm_FI_daily_stage, method = "pairwise", by = "GROUP")
  contrasts_by_group_FI_daily_stage_df <- as.data.frame(contrasts_by_group_FI_daily_stage)
  
  # Pairwise contrasts within each stage of weight cycling
  contrasts_by_SABLE_FI_daily_stage <- contrast(emm_FI_daily_stage, method = "pairwise", by = "Transition")
  contrasts_SABLE_FI_daily_stage_df <- as.data.frame(contrasts_by_SABLE_FI_daily_stage)
  
  ##Graph predicted daily FI by stage (emmeans for MLR) ####
  barplot_emm_FI_daily_stage <- emm_FI_daily_stage_df %>%
    ggplot(aes(x = Transition, y = emmean, fill = GROUP)) +
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
      title = "MLR: daily FI (kcal/day) during weight cycle",
      y = "Food intak (kcal/day)",
      x = "Time point",
      fill = "Treatment group") +
    format.plot +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face="bold", size=15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 11))
  barplot_emm_FI_daily_stage
#Conclusion: During the weight regain stage, FI was not significantly different 
  #for weight cycled and control mice (emmeans, p=0.91).
  #Weight cycled mice had comparable daily FI during obesity development and weight regain (emmeans, p=0.12)
  
  
#--------------------------------- #
  #Work in progress ####
#Cumulative FI during BW regain
#All mice were in regain phase for the same number of days (28), so I don't need to correct 
  #for the number of days
#Create new data frame starting with df called FI_BW_joined
#Create a new type of day_rel called "DAY" --> first day of injections Day =0 and sac day = 28
  For_regain <- FI_BW_joined %>%
    ungroup() %>%
    group_by(ID) %>% 
    arrange(DATE) %>% 
    mutate(
      GROUP = case_when(
        ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
        ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
      DRUG = case_when(
        ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
        ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
    mutate(
      #day_rel = DATE - first(DATE),
           FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
           STATE = case_when(
             ##Baseline: First day of LFD/first day of obesity development
             #Date is the first date after start of LFD for which there was a BW and FI measurement
             ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
             ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
             ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
             ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
             
             ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
             #First day of calorie restriction for Weight Cycled mice. 
             #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
             ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                       3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
             ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                       3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
             ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
             
             ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
             #Sable day 15
             ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
             ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
             ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
             ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
             ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
             
             ##BW maintenance: End of Sable recording for BW maintenance
             #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
             ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
             ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
             ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
             ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
             ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
             ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
             ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
             ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
    
             ##BW regain: Final day of experiment --> day of sac (End of regain)
             ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
             ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
             ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
             ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
             ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
             ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
             ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
             TRUE ~ NA_character_)) %>%
mutate(STATE = factor(STATE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>%
    
For_regain2 <- For_regain %>%
filter((DATE >= "2025-06-12" & ID %in% c(3706, 3707, 3708))) & (DATE >= "2025-06-13" & ID %in% c(3709, 3710, 3711)))) %>%
mutate(try = max(DATE)-DATE,
       try = as.numeric(try),
       try_30 = max(try)-30) %>%
    
filter(try)
    
filter
  
mutate(DAY = case_when(
  ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ 0,
  ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ 0,
  ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ 0,
  ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ 0,
  ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ 0,
  ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ 0,
  ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
  ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
  
  ##BW regain: Final day of experiment --> day of sac (End of regain)
  ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
  ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
  ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
  ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
  ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
  ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
  ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
  TRUE ~ NA_character_))

mutate(DAY = ifelse((SABLE="BW maintenance"), 0, DATE-start(DATE))

For_regain2 <- For_regain %>%
    ungroup() %>%
    group_by(ID) %>%
    arrange(DATE)  %>%
    mutate(DAY = DATE - first(DATE),
           FI_cumsum_kcal =cumsum((INTAKE_GR*3.82)))
   
  