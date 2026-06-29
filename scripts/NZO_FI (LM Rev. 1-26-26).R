#NZO food intake
#Started: 1-26-26
#Revised 1-26-26
#Starting source: CS's script called paperNZOC57courseFI.R
#1-26-26: Removed C57 and updated dates for time points using my Excel file called BWstages
  #added time points called 30 days loss and 20% loss (cumulative intake peak obesity)


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
                                           3717, 3718, 3719, 3726) & DATE == as.Date("2025-04-15") ~ "20 loss",
         
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
