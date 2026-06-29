#NZO food intake
  #Cumulative between time periods (1-26-26) AND during Sable measurements (1-27-26)

#Started: 1-26-26
#Latest revision: 1-27-26
#Starting source: CS's script called paperNZOC57courseFI.R

#1-26-26: Removed C57 and updated dates for time points using my Excel file called BWstages
  #added time points called 30 days loss and 20% loss (cumulative intake peak obesity)
#1-27-26: did not work on the part from 1-26-26 on 1-27-26. Calculated FI at peak obesity for 3706. Used two approaches
#and cross referenced with manual measurements. The two approaches that I used both worked.
  #Left off on line 495 "Extend..."
#1-__-26: Next step -> extend code that I wrote for 3706 at peak obesity to include all IDs at peak obesity
#Note: keep working on the approach that I started on 1-27-26. Work from 1-26-26 is still incomplete

#### Started with Lines 21 through 99 in paperNZOC57courseFI.R and modified to remove C57 and ammend dates for NZO peak obesity
#BW loss and BW maintenance
#FI CSV data import RTIOXA 47 ----
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>%
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
           day_rel == 0 ~ "baseline", 
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714, 
                                             3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724,
                                             3725, 3726, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "peak obesity",
           STRAIN == "NZO/HlLtJ" & ID == 3727 & DATE == as.Date("2025-03-10") ~ "peak obesity",
           
           STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID == 3714 & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           
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
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()

####Based on lines 102 -123 ####

FI_data %>% 
  group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO per STATUS

# Summarize cumulative FI per ID and STATUS ----
FI_stage_summary <- FI_data %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATUS,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `peak obesity` - baseline,
    kcal_peak_to_loss = `BW loss` - `peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`
  )

FI_stage_summary #now we can plot this

#### based on lines 126 through 150 ####
# Convert to long format for plotting ----
FI_stage_long <- FI_stage_summary %>%
  select(ID, STRAIN, DRUG, GROUP,
         kcal_baseline_to_peak,
         kcal_peak_to_loss,
         kcal_loss_to_maint,
         kcal_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("kcal_"),
    names_to = "transition",
    values_to = "kcal"
  ) %>%
  mutate(
    transition = factor(
      transition,
      levels = c("kcal_baseline_to_peak",
                 "kcal_peak_to_loss",
                 "kcal_loss_to_maint",
                 "kcal_maint_to_regain"),
      labels = c("Baseline → Peak Obesity",
                 "Peak Obesity → BW Loss",
                 "BW Loss → BW Maintenance",
                 "BW Maintenance → BW Regain")
    )
  )

#### based on lines 152 through 180 ####
#Note: this is not adjusted for the number of days that each mouse was in each phase
# this is an issue because there was a large range in the time that each mouse was in each phase
#I should make a variable called time_in_phase which is the day_rel(STATUS 1) - day_rel(STATUS 2)
# Plot ---- 
ggplot(FI_stage_long, aes(x = transition, y = kcal, fill = GROUP)) +
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
  #facet_wrap(~ STRAIN*GROUP, scales = "free_y") +
  labs(
    x = "Transition",
    y = "Calories Consumed (kcal)",
    fill = "Feeding Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) 

#### Attempt at different approach by LM. Look at cumulative FI just during the period of BW decrease
#Try this out

#### LM: Look at cumulative FI during the first 30 days of food restriction. How does this compare
  # to what was calculated above for the arbitrary period of BW loss (arbitrary because BW maintenance is weird)
FI_data2 <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>%
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )) %>%
mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
       day_rel = DATE - first(DATE),
       FI_cum =cumsum(corrected_intake_kcal),
       STATUS = case_when(
         day_rel == 0 ~ "baseline", 
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714, 
                                           3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724,
                                           3725, 3726, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "peak obesity",
         STRAIN == "NZO/HlLtJ" & ID == 3727 & DATE == as.Date("2025-03-10") ~ "peak obesity",
         
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714, 
                                           3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724,
                                           3725, 3726, 3728, 3729) & DATE == as.Date("2025-03-19") ~ "30 days loss",
         STRAIN == "NZO/HlLtJ" & ID == 3727 & DATE == as.Date("2025-04-09") ~ "30 days loss",
         
         STRAIN == "NZO/HlLtJ" & ID == 3729 & DATE == as.Date("2025-03-17") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID == 3728 & DATE == as.Date("2025-03-21") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3725) & DATE == as.Date("2025-03-25") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID == 3723 & DATE == as.Date("2025-03-26") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3722) & DATE == as.Date("2025-03-28") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID == 3727 & DATE == as.Date("2025-04-04") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-04-07") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3721) & DATE == as.Date("2025-04-09") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID == 3720 & DATE == as.Date("2025-04-12") ~ "20 loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 
                                           3717, 3718, 3719, 3726) & DATE == as.Date("2025-04-14") ~ "20 loss",
         
         STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
         
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
         STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
         STRAIN == "NZO/HlLtJ" & ID == 3714 & DATE == as.Date("2025-06-19") ~ "BW maintenance",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
         STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
         
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
                         levels = c("baseline", "peak obesity", "30 days loss", "20 loss", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()


FI_data2 %>% 
  group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO per STATUS

# Summarize cumulative FI per ID and STATUS ----
FI_stage_summary2 <- FI_data2 %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATUS,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `peak obesity` - baseline,
    kcal_peak_to_loss = `BW loss` - `peak obesity`,
    kcal_peak_to_30_days = `30 days loss` - `peak obesity`,
    kcal_peak_to_20_percent = `20 loss` - `peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`
  )

FI_stage_summary2 #now we can plot this

# Convert to long format for plotting ----
FI_stage_long2 <- FI_stage_summary2 %>%
  select(ID, STRAIN, DRUG, GROUP,
         kcal_baseline_to_peak,
         kcal_peak_to_loss,
         kcal_peak_to_30_days,
         kcal_peak_to_20_percent,
         kcal_loss_to_maint,
         kcal_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("kcal_"),
    names_to = "transition",
    values_to = "kcal"
  ) %>%
  mutate(
    transition = factor(
      transition,
      levels = c("kcal_baseline_to_peak",
                 "kcal_peak_to_loss",
                 "kcal_peak_to_30_days",
                 "kcal_peak_to_20_percent",
                 "kcal_loss_to_maint",
                 "kcal_maint_to_regain"),
      labels = c("Baseline → Peak Obesity",
                 "Peak Obesity → BW Loss",
                 "Peak Obesity → 30 days Loss",
                 "Peak Obesity → 20% Loss",
                 "BW Loss → BW Maintenance",
                 "BW Maintenance → BW Regain")
    )
  )

# Plot ---- 
FI_stage_long3 <-FI_stage_long2 %>%
  filter(DRUG == "vehicle")
ggplot(FI_stage_long3, aes(x = transition, y = kcal, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "median",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3) +
  #facet_wrap(~ STRAIN*GROUP, scales = "free_y") +
  labs(
    x = "Transition",
    y = "Calories Consumed (kcal)",
    fill = "Feeding Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) 

#Try looking at FI for all mice between the day after Sable day 10 (i.e. peak obesity sable) and before acclimation for sable day 11
mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
       day_rel = DATE - first(DATE),
       FI_cum =cumsum(corrected_intake_kcal),
       STATUS = case_when(
         day_rel == 0 ~ "baseline", 
         STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714, 
                                           3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724,
                                           3725, 3726, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "peak obesity",
         STRAIN == "NZO/HlLtJ" & ID == 3727 & DATE == as.Date("2025-03-10") ~ "peak obesity",

####
##### 1-27-26 (calculating Sable FI) ####
#Start of 1/27attempt using section from locomotion_LM_12-2.R to do initial processing of sable_dwn to select just FoodA
sable_FI_data <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light"),
    SABLE = case_when(
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                               "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                               "SABLE_DAY_7") ~ "Baseline",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_)) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("FoodA_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  group_by(ID, DRUG,GROUP,DIET_FORMULA,SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,DIET_FORMULA,SABLE) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter( is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain")))
#This created a df from sable_dwn with only FoodA under the "value" column

####This method worked (manually measured intake =10.42g and Sable measured 10.89g) ####
#Filter for just 3706 on the first complete day of SABLE=Peak obesity
sable_FI_3706 <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  filter(complete_days == "1") %>%
  arrange(hr) %>%
  ungroup() %>%
  group_by(hr) %>%
  mutate (hr_FI = first(value) - value) #%>%
  #mutate (hour_FI = max(hr_FI))
  #could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706 <- sable_FI_3706 %>%
  group_by(ID, GROUP, DRUG, SABLE, hr) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706 <- hr_FI_3706 %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
    .groups = "drop")
  
#Filter for just 3706 on the second complete day of SABLE=Peak obesity
sable_FI_3706_2 <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  filter(complete_days == "2") %>%
  arrange(hr) %>% #should I arrange by DateTime instead since hrs 20-23 actually occurred before hrs 0-19?
  ungroup() %>%
  group_by(hr) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_2 <- sable_FI_3706_2 %>%
  group_by(ID, GROUP, DRUG, SABLE, hr) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706_2 <- hr_FI_3706_2 %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
            .groups = "drop")

####To get the sum of FI during day 1 and day 2 using one chunk of code I grouped by complete_days ####
#this approach worked (sum of days 1 and 2 = 10.89)
#Filter for 3706 and sable day 9
sable_FI_3706_bothdays <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == "1") %>% # commenting this out leaves both complete days in
  arrange(hr, date) %>% 
  ungroup() %>%
  group_by(hr, date) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_bothdays <- sable_FI_3706_bothdays %>%
  group_by(ID, GROUP, DRUG, SABLE, hr, complete_days) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),
    .groups = "drop") 

sum_3706_bothdays <- hr_FI_3706_bothdays %>%
  summarise(day_cum_FI = sum(FI_cum_hr),  
            .groups = "drop")

#Here I used the same approach as immediately above, but I arranged by DateTime. I think this will be better
  #for when I want to expand this code to be used for multiple time points

sable_FI_3706_cum <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == 1) %>% # commenting this out leaves both complete days in
  arrange(DateTime) %>%
  ungroup() %>%
  group_by(hr, date) %>%
  mutate (hr_FI = first(value) - value) #%>%
#mutate (hour_FI = max(hr_FI))
#could slice off the last entry for each hr if I use the mutate line to find the hrly FI

hr_FI_3706_cum <- sable_FI_3706_cum %>%
  group_by(ID, GROUP, DRUG, SABLE, hr, complete_days) %>%
  summarise(
    FI_cum_hr = max(hr_FI, na.rm = TRUE),.groups = "drop") 

sum_3706_cum <- hr_FI_3706_cum %>%
  summarise(day_cum_FI = sum(FI_cum_hr),.groups = "drop")

#### Different approach: Sum minute by minute FI to obtain FI for complete day 1 and 2  ####
#this works
sable_FI_3706_minutes <- sable_FI_data %>%
  group_by(ID, complete_days, DRUG, SABLE) %>%
  filter(ID == 3706) %>%
  filter(SABLE == "Peak obesity") %>%
  #filter(complete_days == "2") %>% #Use this line to calculate just one of the two complete days
  arrange(ID, DateTime, hr, complete_days) %>%       # make sure data is ordered
  group_by(ID, complete_days, SABLE, GROUP, DRUG) %>% 
  mutate(
    intake = lag(value) - value,       # change in meters
    intake = if_else(intake < 0, 0, intake)) %>%  
    #the line above says that if the calculated intake is less than zero then put 0 in the cell. 
    #Otherwise, put the intake (g) which was calculated into the cell.
  drop_na() %>% 
  summarise(
    total_eaten_gr = sum(intake),      # total FI per day in grams
    .groups = "drop") 

#### Extend the code that I used for just 3706 at peak obesity for all mouse IDs at peak obesity ####

  
