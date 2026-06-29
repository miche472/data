# Major topics: Combine BW and FI csv files; developing a way to accurately indicate 
    # transitions between time periods (i.e. summarize phases)

#Started: 2-16-26
#Revised: 2-19-26

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

#join BW and FI data
FI_BW_joined <- FI_to_join %>%
  left_join(
    BW_to_join %>% 
      select(ID, DATE, BW, COMMENTS_BW),
    by = c("ID", "DATE"))

#Add variables: GROUP, DRUG, STATE and day_rel
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
  mutate(day_rel = DATE - first(DATE)) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATE = case_when(
##Baseline: First day of LFD (first day of obesity development/start of ad lib LFD)
    #Start of ad libitum LFD
    ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
    ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
    ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
    ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
  
##Peak obesity: End peak obesity period (last day of Peak obesity sable)
    #Sable day 10 or 11
      #ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "Peak obesity", 
      #ID %in% c(3719, 3726) & DATE == as.Date("2025-02-10") ~ "Peak obesity",
      #ID %in% c(3711, 3718) & DATE == as.Date("2025-03-03") ~ "Peak obesity",
    #First day of restricted diet (for weight cycled mice. Date matching for Control mice)
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
    #Sable day 19
    #STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "BW maintenance", 
    #STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719) & DATE == as.Date("2025-06-08") ~ "BW maintenance",
    #STRAIN == "NZO/HlLtJ" & ID %in% c(3726) & DATE == as.Date("2025-06-12") ~ "BW maintenance",

    #First day of injections for all mice and restoration of ad libitum LFD for weight cycled mice
    ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
    ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
    ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
    ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
    ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
    ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
    ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
    ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
  
##BW regain: Final Sac (End of regain)
    ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
    ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
    ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
    ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
    ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
    ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
    ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
TRUE ~ NA_character_)) %>%
  filter(!is.na(STATE)) %>%
  mutate(STATE = factor(STATE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 

#Apply logic from FI course to make phases. The difference is that the dates that I have created are more accurate
  #than what is in the FI course script

#Verify that there are 22 mice in each of the 5 states
df1 %>% 
  group_by(STATE) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO per STATE


# Summarize cumulative FI per ID and STATE ----
FI_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `Peak obesity` - Baseline,
    kcal_peak_to_loss = `BW loss` - `Peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format for plotting ----
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
      labels = c("Baseline → Peak obesity",
                 "Peak obesity → BW loss",
                 "BW loss → BW maintenance",
                 "BW maintenance → BW regain")))

#DURATION: calculate the number of days in each transition period ####
Days_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(max_day_rel = max(day_rel, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = max_day_rel) %>%
  # Calculate kcal consumed between stages
  mutate(
    Days_baseline_to_peak = `Peak obesity` - Baseline,
    Days_peak_to_loss = `BW loss` - `Peak obesity`,
    Days_loss_to_maint = `BW maintenance` - `BW loss`,
    Days_maint_to_regain = `BW regain` - `BW maintenance`)

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
      labels = c("Baseline → Peak obesity",
                 "Peak obesity → BW loss",
                 "BW loss → BW maintenance",
                 "BW maintenance → BW regain")))

# Join duration (days) and FI during each phase
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


#How is it possible that the food restricted mice ate more during BW loss? 
#We saw that the control mice lost a lot of weight despite no food restriction
#so maybe looking at food efficiency is necessary to make sense of this.
#2-18-26: I think that this could be due to the use of corrected_FI during this time period
#2-20-26: I am now certain that my thought on 2-18-26 was correct. Refer to:
#FI using INTAKE_GR (Rev...).R for accurate/corrected approach

# Summarize change in BW per ID and STATE ----
#Used "paperNZOC57courseFI.R" to guide connecting of BW change and FI
#This section is essentially taken directly from the script cited and only slightly modified
#by me. The way that cummulative FI is calculated and that energy efficiency is calculated
#is incorrect, so the outcomes below are illogical. Don't use this section.

FI_stage_LM <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
summarise(
  FI_cum_end = max(FI_cum, na.rm = TRUE),
  day_end = max(day_rel, na.rm = TRUE),
  .groups = "drop")

FI_stage_LM_LM <- FI_stage_LM %>%
  left_join(
    df1 %>% 
      select(ID, STATE, BW),
    by = c("ID", "STATE"))

# Create transition columns
df_transitions_LM <- FI_stage_LM_LM %>%
  arrange(ID, day_end) %>%
  group_by(ID) %>%
  mutate(
    STATE_TRANSITION = paste(STATE, "->", lead(STATE)),
    BW_change = lead(BW) - BW,
    delta_days = lead(day_end) - day_end,
    delta_FI_cum = lead(FI_cum_end) - FI_cum_end,) %>%
  # Remove last row (no transition after final STATE)
  filter(!is.na(BW_change)) %>%
  ungroup()
  
df_transitions_LM %>% #for purposes of viewing data
  select(ID, STATE_TRANSITION, BW_change, delta_days)

df_transitions_LM <- df_transitions_LM %>% #for purposes of selecting data to proceed with
  select(ID, GROUP, DRUG, STATE_TRANSITION, BW_change, delta_days, delta_FI_cum)

df_transitions_LM_LM <- df_transitions_LM %>%
  mutate(
    delta_days_num = as.numeric(gsub(" days", "", delta_days)),
    FI_kcal_per_day = delta_FI_cum / delta_days_num,
    FI_kcal_per_BWchange = ifelse(abs(BW_change) > 0.001, delta_FI_cum / BW_change, NA),
    energy_efficiency = ((BW_change /delta_FI_cum)/delta_days_num),
    DRUG = factor(DRUG)) 

# Define order for the transitions
STATE_order <- c(
  "Baseline -> Peak obesity",
  "Peak obesity -> BW loss",
  "BW loss -> BW maintenance",
  "BW maintenance -> BW regain")

# Prepare data
df_transitions_LM_LM_LM <- df_transitions_LM_LM %>%
  mutate(
    delta_days_num = as.numeric(gsub(" days", "",delta_days)), #
    FI_kcal_per_day = delta_FI_cum / delta_days_num,
    STATE_TRANSITION = factor(STATE_TRANSITION, levels = STATE_order)) 

# Bar plot of daily food intake during each of the 4 phases
ggplot(df_transitions_LM_LM_LM, aes(x = STATE_TRANSITION, y = FI_kcal_per_day, fill = GROUP)) +
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
    x = "Phase transition",
    y = "Food intake (kcal/day)",
    fill = "Treatment group")

# Bar plot with FI_kcal_per_BWchange and free y-axis
ggplot(df_transitions_LM_LM_LM, aes(x = STATE_TRANSITION, y = FI_kcal_per_BWchange, fill = GROUP)) +
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
    x = "Status transition",
    y = "Food intake (kcal / g BW change)",
    fill = "Treatment group")
