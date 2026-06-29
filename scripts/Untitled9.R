#To get average intake during BW loss and then BW maintenance use SABLE_DAY_"X" as sign posts
#For Baseline: average corrected_intake_kcal from DATE of SABLE_DAY_1 through SABLE_DAY_7
#For Peak obesity: average corrected_intake_kcal from DATE of day after SABLE_DAY_7 through SABLE_DAY_11
#For BW loss: average corrected_intake_kcal from DATE of day after SABLE_DAY_11 (first day of BW loss caloric restriction) through 
    #DATE of SABLE_DAY_15 (last day of SABLE= BW loss)
#For BW maintenance: average corrected_intake_kcal from DATE of day after SABLE_DAY_15 (first day of BW maintenance caloric restriction) 
    #through DATE of SABLE_DAY_19 (last day of SABLE= BW maintenace)
#For BW regain: average corrected_intake_kcal from DATE after SABLE_DAY_19 through DATE of SABLE_DAY_23

#I want to calculate average_corrected_intake_kcal for each mouse ID during the following time periods:
#Baseline (the avearge of SABLE_DAY_1 through SABLE_DAY_7)...call this SABLE= "Baseline"
#Baseline -> Peak obesity...call this SABLE= "Peak obesity"
#Peak obesity -> BW loss...call this SABLE= "BW loss"
#BW loss -> BW maintenace...call this SABLE= "BW maintenance"
#BW maintenance -> BW regain...call this SABLE= "BW regain"

FI_data_LM <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #This chow and LFD identified, as applicable 
  select(-DIET_FORMULA.y) %>% #this assumes that all mice had LFD for all measurements. this is wrong so delete the column.
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(
    SABLE_DAY = case_when(
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
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23")) 


FI_3706 <- FI_data_LM %>%
  filter(ID=="3706") %>%
mutate(
  date_sable = case_when( #create a column called "SABLE" and assign all dates that fall between
    SABLE_DAY %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline"))

if_then(DATE="2025-07-05")
    
    
))
date=

#start by doing this just for 3706


mutate(
  SABLE_DAY = case_when(
    #Baseline (sable days 1-7)
    STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
    

#####try ####
    sable_dates <- FI_data_LM %>%
      filter(!is.na(SABLE_DAY)) %>%
      select(ID, SABLE_DAY, DATE) %>%
      arrange(ID, DATE)

sable_ranges <- sable_dates %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    next_DATE = lead(DATE),
    range_start = DATE,
    range_end = next_DATE
  ) %>%
  filter(!is.na(next_DATE)) %>%   # remove last SABLE day (no end)
  select(ID, SABLE_DAY, range_start, range_end)
    
