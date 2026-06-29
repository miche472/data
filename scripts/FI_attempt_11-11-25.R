#Average daily food intake for each diet group during each phase
#1. Calculate cumulative food intake for the phase. 
#2. Calculate total number of days in the phase
#3. Divide cumulative food intake for the phase/number of days in phase
#4. Graph each diet condition at each time point

#
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


####---- Calculate cumulative food intake for each phase ----####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
           day_rel == 0 ~ "Baseline", 
           STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "Peak obesity",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) & DATE == as.Date("2025-04-04") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725) & DATE == as.Date("2025-04-18") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-04-25") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
           STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3715, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-21") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-25") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID == 3726 & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719, 3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_
         )) %>% 
  filter(!is.na(STATUS)) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()

FI_data %>% 
  group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO

# Summarize cumulative FI per ID and STATUS ----
FI_stage_summary <- FI_data %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATUS,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages (CS way)
  mutate(
    kcal_baseline_to_peak = `Peak obesity` - Baseline,
    kcal_peak_to_loss = `BW loss` - `Peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`)

FI_stage_summary #now we can plot this

# Convert to long format for plotting ----
FI_stage_long <- FI_stage_summary %>%
  select(ID, STRAIN, DRUG, GROUP,
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

# Graph cumulative FI for each phase ---- 
ggplot(FI_stage_long, aes(x = Transition, y = kcal, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3) +
  labs(
    x = "Transition",
    y = "Calories consumed (kcal)",
    fill = "Treatment group") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)) 

#### Calculate number of days in each phase ####

#in FI_data for each mouse ID order the observations by date
#then calculate: "days = day_rel -lag (day_rel)" 

####---- Calculate cumulative food intake for each phase ----####
#Create FI_data and add cohort # for each mouse based on ID
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
           day_rel == 0 ~ "Baseline", 
           STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "Peak obesity",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) & DATE == as.Date("2025-04-04") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725) & DATE == as.Date("2025-04-18") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-04-25") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
           STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3715, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-21") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-25") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID == 3726 & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719, 3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>% 
  filter(!is.na(STATUS)) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()

FI_data %>% 
  group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO

####start delete ####
# Summarize cumulative FI per ID and STATUS ----
FI_stage_summary <- FI_data %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>% 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATUS,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages (CS way)
  mutate(
    kcal_baseline_to_peak = `Peak obesity` - Baseline,
    kcal_peak_to_loss = `BW loss` - `Peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`)

FI_stage_summary #now we can plot this

# Convert to long format for plotting ----
FI_stage_long <- FI_stage_summary %>%
  select(ID, STRAIN, DRUG, GROUP,
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
      labels = c("Baseline → Peak Obesity",
                 "Peak Obesity → BW Loss",
                 "BW Loss → BW Maintenance",
                 "BW Maintenance → BW Regain"))) %>%
  mutate(
    COHORT = case_when(
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) ~ 3,
      ID %in% c(3716, 3717, 3718, 3719, 3720,3721, 3722, 3723, 3724, 3725) ~ 4,
      ID %in% c(3726, 3727, 3728, 3729) ~ 5))
####end delete ####

#### start delete ####
####-try---#####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
           day_rel == 0 ~ "Baseline", 
           STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "Peak obesity",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) & DATE == as.Date("2025-04-04") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725) & DATE == as.Date("2025-04-18") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-04-25") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
           STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3715, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-21") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-25") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID == 3726 & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719, 3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_
         )) %>% 
  filter(!is.na(STATUS)) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()
####end delete ####

#### try (actual) ####
# Summarize duration per ID and STATUS ----

Days_transition_long <- FI_data %>%
  group_by(ID) %>%
  mutate(day= day_rel -lag(day_rel)) %>%
  select(ID, DRUG, GROUP, STATUS, day) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(Transition = 
  case_when(
    STATUS=="Peak obesity" ~ "Baseline → Peak obesity",
    STATUS=="BW loss" ~ "Peak obesity → BW loss",
    STATUS=="BW maintenance" ~ "BW loss → BW maintenance",
    STATUS=="BW regain" ~ "BW maintenance → BW regain")) %>%
  select(ID, DRUG, GROUP, day, Transition) %>%
  drop_na(Transition) #this gets rid of the Baseline values since they are actually reflected by the Baseline -> Peak obesity
    
#### Combine days and food intake into one data frame ####
  
FI_Days <- FI_stage_long %>%
  left_join(
    Days_transition_long %>% select(ID, DRUG, GROUP, day, Transition),
    by = c("ID", "Transition"))
  
Time_FI_calc <-FI_Days %>%
  rename(GROUP = GROUP.x) %>% #There is no differences between columns x and y. 
  select(-GROUP.y) %>% 
  rename(DRUG = DRUG.x) %>% #There is no differences between columns x and y. 
  select(-DRUG.y) %>% 
  mutate (day = as.numeric(day)) %>% 
  ungroup() %>%
  group_by(ID, Transition) %>%
  mutate(mouse_FI_day = kcal/day) %>%
  ungroup() 

#Average kcal/day during each transition for each GROUP
Avg_FI <- Time_FI_calc %>%
  group_by(Transition, GROUP) %>%
  summarise(
    mean_kcal = mean(mouse_FI_day, na.rm = TRUE),
    sem_kcal  = sd(mouse_FI_day, na.rm = TRUE) / sqrt(n()),
    n = n(), .groups = "drop") #drop is the same as ungroup
  

#T-test to compare GROUPs at each Transition
  
  
#Repeated measures ANOVA to compare within GROUP between Transition periods
  
  
#Bar graph of average calories consumed per day at each Transition for each group
  
  
  
  
  ####start delete -------####
  ungroup() %>%
  pivot_wider(
    names_from = STATUS,
    values_from = day) %>%
 
#this part creates columns with em
mutate(
  day_baseline_to_peak = `Peak obesity` - Baseline,
  day_peak_to_loss = `BW loss` - `Peak obesity`,
  day_loss_to_maint = `BW maintenance` - `BW loss`,
  day_maint_to_regain = `BW regain` - `BW maintenance`)

mutate(
  Transition = factor(Transition,
                      levels = c("Baseline → Peak obesity", "Peak obesity → BW Loss","BW loss → BW maintenance","BW maintenance → BW regain"),
                      labels = c("Baseline → Peak obesity","Peak obesity → BW loss","BW loss → BW maintenance","BW maintenance → BW regain"))) %>% 

####end delete ####
  