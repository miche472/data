
#Major topics: Use INTAKE_GR rather than corrected_intake_kcal to calculate FI ####

#Date started: 2-19-26
#Revised

#Objective was to explore:Does cumulative FI (calculated as sum of INTAKE_GR)/total number of days, 
#yield a more accurate quantification of FI compared to corrected_intake_kcal?
  #If you sum corrected_intake_kcal I don't think you actually get cumulative intake. 
  #I think you need to add INTAKE_GR to get cumulative intake. 

#libraries----
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

#--------------------------#
# start delete ####
#First Try calculating cumulative FI using INTAKE_GR
try_cum_FI <- FI_BW_joined %>%
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
  #filter(DATE>= as.Date("2025-02-24") & DATE<=as.Date("2025-06-27")) %>%
  #summarise(sum_INTAKE = sum(INTAKE_GR, na.rm = TRUE), .groups = "drop")
mutate(FI_cum_try =cumsum(INTAKE_GR)) %>%
mutate(INTAKE__GR_rel = INTAKE_GR - first(INTAKE_GR),
       day_rel = DATE - first(DATE),
       FI_cum =cumsum(corrected_intake_kcal),
       mutate(bw_rel = 100 * (BW - first(BW)) / first(BW))
       
       #try graphing this method
       ggplot(try_cum_FI, aes(x = GROUP, y = sum_INTAKE)) +
         geom_bar(stat = "summary", 
                  fun = "mean", 
                  position = position_dodge(width = 0.9)) +
         geom_errorbar(stat = "summary", 
                       fun.data = mean_se, 
                       position = position_dodge(width = 0.9), 
                       width = 0.3) +

# end delete ####
#--------------------------#
       
df_hi <- FI_BW_joined %>%
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
        day_rel = DATE - first(DATE),
        cum_INTAKE_GR = cumsum(INTAKE_GR),
        INTAKE_GR_rel = cum_INTAKE_GR - first(cum_INTAKE_GR),
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
       FI_stage_summary2 <- df_hi %>%
         group_by(ID, GROUP, DRUG, STATE) %>%
         summarise(INTAKE_cum_end_GR = max(cum_INTAKE_GR, na.rm = TRUE), .groups = "drop") %>%
         # Reshape into wide format: one row per ID, columns = each STATUS
         pivot_wider(
           names_from = STATE,
           values_from = INTAKE_cum_end_GR) %>%
         # Calculate kcal consumed between stages
         mutate(
           GR_baseline_to_peak = `Peak obesity` - Baseline,
           GR_peak_to_loss = `BW loss` - `Peak obesity`,
           GR_loss_to_maint = `BW maintenance` - `BW loss`,
           GR_maint_to_regain = `BW regain` - `BW maintenance`)
       
       # Convert to long format
       FI_stage_long2 <- FI_stage_summary2 %>%
         select(ID, DRUG, GROUP,
                GR_baseline_to_peak,
                GR_peak_to_loss,
                GR_loss_to_maint,
                GR_maint_to_regain) %>%
         pivot_longer(
           cols = starts_with("GR_"),
           names_to = "Transition",
           values_to = "GR") %>%
         mutate(
           Transition = factor(
             Transition,
             levels = c("GR_baseline_to_peak",
                        "GR_peak_to_loss",
                        "GR_loss_to_maint",
                        "GR_maint_to_regain"),
             labels = c("Baseline → Peak obesity",
                        "Peak obesity → BW loss",
                        "BW loss → BW maintenance",
                        "BW maintenance → BW regain"))) %>%
         mutate(INTAKE_kcal = GR*3.82)
       
#Duration of each transition: use df1 to calculate duration (days) for each transition period ####
       # Summarize duration (days) per ID and STATE
       Days_stage_summary2 <- df_hi %>%
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
       Days_stage_long2 <- Days_stage_summary2 %>%
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
             labels = c("Baseline → Peak obesity",
                        "Peak obesity → BW loss",
                        "BW loss → BW maintenance",
                        "BW maintenance → BW regain")))
       
       # df2: Join duration (days) and cummulative FI during each phase
       FI_duration_joined2 <- FI_stage_long2 %>%
         left_join(
           Days_stage_long2 %>% 
             select(Transition, ID, Duration_day),
           by = c("ID", "Transition"))
       
       df22 <- FI_duration_joined2 %>%
         group_by(GROUP, Transition, ID) %>%
         mutate(Daily_INTAKE_kcal = INTAKE_kcal/as.numeric(Duration_day))
       
       #Verify that there are 22 mice in each of the 4 transition periods 
       df22 %>% 
         group_by(Transition) %>%
         summarise(n_ID = n_distinct(ID)) #this we have 22 NZO
       
#ΔBW during each transition: use df1 to calculate delta_BW (g) for each transition period ####
       # Summarize Δ BW (g) per ID and STATE
       BW_stage_summary2 <- df_hi %>% #df1 has 5 values for each mouse (one per time point)
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
       BW_stage_long2 <- BW_stage_summary2 %>%
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
             labels = c("Baseline → Peak obesity",
                        "Peak obesity → BW loss",
                        "BW loss → BW maintenance",
                        "BW maintenance → BW regain")))
       
       #Combine FI, BW change, and duration of each transition phase
       # Add BW_stage_long to df2. df2 has cumulative FI and # of days in each transition period
       
       FI_duration_BW_joined2 <- df22 %>%
         left_join(
           BW_stage_long2 %>% 
             select(Transition, ID, delta_BW_g),
           by = c("ID", "Transition"))
       
       # Calculate change in BW per kcal consumed (phase_BW_FI) for each ID during each Transition
       df42 <- FI_duration_BW_joined2 %>%
         group_by(GROUP, Transition, ID) %>%
         mutate(phase_BW_FI = delta_BW_g/INTAKE_kcal)
       
       #For df4: Verify that there are 22 mice in each of the 4 transition periods 
       df42 %>% 
         group_by(Transition) %>%
         summarise(n_ID = n_distinct(ID)) #we have 22 NZO in each transition period

#Graph Energy Efficiency for Control and Weight cycled mice during each Transition ####
       #mean change in BW (g) per kcal consumed for each GROUP during every Transition
       ggplot(df42, aes(x = Transition, y = phase_BW_FI, fill = GROUP)) +
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
         labs(
           title="Energy efficiency during weight cycling (ΔBW /kcal)",
           x = "Phase of weight cycle",
           y = "ΔBW (g) per kcal consumed",
           fill = "Treatment group")
       
#Graph FOOD INTAKE for Control and Weight cycled mice during each Transition ####
       #Total FI during each phase, by Treatment group
       ggplot(df42, aes(x = Transition, y = INTAKE_kcal, fill = GROUP)) +
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
         labs(
           title="INTAKE during weight cycling (INTAKE_GR method)",
           x = "Phase of weight cycle",
           y = "Intake (kcal)",
           fill = "Treatment group")
       
#Graph Daily FOOD INTAKE for Control and Weight cycled mice during each Transition ####
  #Daily food intake calculated from INTAKE_GR rather than corrected_intake_kcal
       ggplot(df42, aes(x = Transition, y = Daily_INTAKE_kcal, fill = GROUP)) +
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
         labs(
           title="Daily FI during weight cycling (INTAKE_GR method)",
           x = "Phase of weight cycle",
           y = "Average daily intake (kcal)",
           fill = "Treatment group")