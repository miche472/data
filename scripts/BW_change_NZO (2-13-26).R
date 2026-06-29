#2.   Food efficiency (left join data frames below) ####
#For Baseline, Peak obesity, and BW regain look at BW change during sable recording period
#Cumulative FI during sable recording/delta BW from first BW to final BW measurement
#BW measure could be manual from BW.csv or from sable_dwn

#last revised on 2-13-26 (see note about where I left off)

BW_manual <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
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
  mutate(
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23")) %>% 
  mutate(
    SABLE = case_when(
      sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
      sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
      sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
      sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
      sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain", 
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% #This removes days that aren't during sable recording
  mutate(SABLE = factor(SABLE, 
                        levels = c("Baseline", "Peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) 
  
BW_manual %>% 
  group_by(STRAIN,SABLE) %>%
  summarise(n_ID = n_distinct(ID)) 


BW_stage_summary <- BW_manual %>%
  group_by(ID, STRAIN, GROUP, DRUG, SABLE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each SABLE
  pivot_wider(
    names_from = SABLE,
    values_from = BW_end) %>%
  # Calculate kcal consumed between stages
  mutate(
    gram_baseline_to_peak = `Peak obesity` - Baseline,
    gram_peak_to_loss = `BW loss` - `Peak obesity`,
    gram_loss_to_maint = `BW maintenance` - `BW loss`,
    gram_maint_to_regain = `BW regain` - `BW maintenance`)
  
BW_stage_long <- BW_stage_summary %>%
  select(ID, STRAIN, DRUG, GROUP,
         gram_baseline_to_peak,
         gram_peak_to_loss,
         gram_loss_to_maint,
         gram_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("gram_"),
    names_to = "Transition",
    values_to = "BW_change_gram") %>%
mutate(
  Transition = factor(
    Transition,
    levels = c("gram_baseline_to_peak",
               "gram_peak_to_loss",
               "gram_loss_to_maint",
               "gram_maint_to_regain"),
    labels = c("Baseline → Peak Obesity",
               "Peak Obesity → BW Loss",
               "BW Loss → BW Maintenance",
               "BW Maintenance → BW Regain")))

# Plot ---- 
ggplot(BW_stage_long, aes(x = Transition, y = BW_change_gram, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3) +
  labs(title= "Change in body weight between indirect calorimetry measurements (NZO)",
    x = "Weight cycle phase",
    y = "Change in body weight (grams)",
    fill = "Treatment group") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)) 

#Now try calculating the change between each BW measurement. Then sum the total of the changes in BW between measurements
#First try to get the difference from measurement to measurement without delineating between phases

BW_delta_cum <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
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
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_))
  #Once I am able to get Sable appended to each measurement, I Could add SABLE phases here. 
    #Then group by SABLE and ID. This might give me the sum of change in BW for each ID
  
 BW_working <- BW_delta_cum %>%
  mutate(day_rel = DATE - first(DATE)) %>%
          #filter (COMMENTS == "DAY_7_SABLE_DAY_1_LFD_AD_LIB") #filters to help with manually finding ranges for day_rel
          #filter(sable_idx %in% c("SABLE_DAY_10", "SABLE_DAY_11"))
   
          #filter(sable_idx %in% c("SABLE_DAY_10", "SABLE_DAY_11")) # start of Peak obesity --> BW loss (control mice)
          #filter(COMMENTS=="RESTRICTED_DAY_1") # start of Peak obesity --> BW loss (weight cycled mice)
          #filter(sable_idx %in% c("SABLE_DAY_14", "SABLE_DAY_15")) # end of Peak obesity --> BW loss (all mice)
          
          #filter(sable_idx %in% c("SABLE_DAY_14", "SABLE_DAY_15")) #start of BW loss --> BW maintenance
          #filter(sable_idx %in% c("SABLE_DAY_18", "SABLE_DAY_19")) #end of BW loss --> BW maintenance (control mice)
          #filter(COMMENTS=="DAY_1_INJECTIONS") #end of BW loss --> BW maintenance (restricted mice)
          
          #filter(COMMENTS=="DAY_1_INJECTIONS") # start of BW maintenance --> BW regain
          #filter(COMMENTS=="DAY_4_SABLE_AND_SAC") # end of BW maintenance --> BW regain
#### LEFT OFF working on this section (2-13-26) ####
   mutate(STATUS = case_when(
       day_rel == 0 ~ "Baseline",
       ID == 3706 & day_rel > 24 & day_rel <94 ~ "Baseline → Peak Obesity",
       ID == 3706 & day_rel > 24 & day_rel <94 ~ "Baseline → Peak Obesity",
       
       ID == 3706 & day_rel > & day_rel < ))
 
 
#this works, but I need to do the part above first 
mutate(BW_measure_change = BW-lag(BW)) %>%
  #mutate(cum_sum_BW = sum(BW_measure_change))
  drop_na() %>%
summarise(
  BW_sum = sum(BW_measure_change), .groups = "drop")
  
Phase = case_when(
  sable_idx = 7 



#--------------------------------------------------------------------
#### delete code below...this approach doesn't work ####
#For 3706, graph BW over time
BW_3706 <- BW_manual %>%
  mutate(ID = factor(ID)) %>%
  mutate(
    sable_day_numeric = case_when(
     sable_idx == "SABLE_DAY_1" ~ 1,  
     sable_idx == "SABLE_DAY_2" ~ 2,
     sable_idx == "SABLE_DAY_3" ~ 3,
     sable_idx == "SABLE_DAY_4" ~ 4,
     sable_idx == "SABLE_DAY_5" ~ 5,
     sable_idx == "SABLE_DAY_6" ~ 6,
     sable_idx == "SABLE_DAY_7" ~ 7, # End baseline phase
     sable_idx == "SABLE_DAY_8" ~ 8,   
     sable_idx == "SABLE_DAY_9" ~ 9,
     sable_idx == "SABLE_DAY_10" ~ 10,
     sable_idx == "SABLE_DAY_11" ~ 11, #End Peak obesity phase
     sable_idx == "SABLE_DAY_12" ~ 12, 
     sable_idx == "SABLE_DAY_13" ~ 13,
     sable_idx == "SABLE_DAY_14" ~ 14,
     sable_idx == "SABLE_DAY_15" ~ 15, #End BW loss phase
     sable_idx == "SABLE_DAY_16" ~ 16, 
     sable_idx == "SABLE_DAY_17" ~ 17,
     sable_idx == "SABLE_DAY_18" ~ 18,
     sable_idx == "SABLE_DAY_19" ~ 19, #End BW maintenance phase
     sable_idx == "SABLE_DAY_20" ~ 20, 
     sable_idx == "SABLE_DAY_21" ~ 21,
     sable_idx == "SABLE_DAY_22" ~ 22,
     sable_idx == "SABLE_DAY_23" ~ 23)) %>% #End BW regain phase
  filter(ID == "3706") %>%
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>% #order data by date for each mouse ID
  mutate(
    sable_day_filled = zoo::na.locf(sable_day_numeric, na.rm = FALSE)) %>% #extends the most recent value into NA cells
  ungroup() %>%
  mutate(phase = case_when(
      sable_day_filled <= 7  ~ "Baseline",
      sable_day_filled <= 11 ~ "Peak obesity",
      sable_day_filled <= 15  ~ "BW loss",
      sable_day_filled < 19 ~ "BW maintenance",
      sable_day_filled > 19 & sable_day_filled <= 23 ~ "BW regain",
      TRUE~ NA_character_))


  
