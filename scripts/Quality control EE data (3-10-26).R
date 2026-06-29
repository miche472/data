
#Quality control of EE data --> look at energy balance vs BW change during sable
#If energy balance is positive then mice should have gained weigh during sable
#This approach to quality control for EE data was suggested in the consenus paper for iCAL

#Script started: 3-10-26
#Author: Lauren Michels

#Conclusion: Sable food sensors were not calibrated and therefore the measurements of 
#FI in sable are not completely reliable. Manual FI while mice were in Sable were not 
#consistent for all mice so they aren't a perfect stand in. The same goes for BW in Sable
#In the future, this script will be a useful guide for EE data quality control

#Extract from FI code:
#Dates corresponding to each sable day (1 to 23) for each mouse
#Sable days corresponding to each sable phase
#Ideally there will be at least 160 observations (16 mice x 5 periods x 2 BW measures)
BW_data <- read_csv("../data/BW.csv") %>% 
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
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
    sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
    TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance", "BW regain"))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  ungroup() 

#Once data is loaded, group by ID and arrange by date. mutate(delta_BW= Weight - first(Weight) 
#This will result in 5 observations/mouse unless a mouse was weighed on >2 sable days during a single time point

BW_change_data <- BW_data %>%
  ungroup() %>%
  mutate(ID = factor(ID)) %>%
  group_by(ID, SABLE) %>%
  arrange(DATE)%>%
  summarise(
    first_bw = first(BW),
    last_bw  = last(BW),
    bw_change = last_bw - first_bw,
    .groups = "drop") %>%
  mutate(direction_bw = if_else((bw_change >0), "Gained", "Lost"))

#Quality control check:
#Alternative method daily energy balance (manual food for all mice during sable)
#Confirm that energy balance direction and BW change direction are the same (gain = positive, loss=negative)
  
#df = Alt_E_balance_kcal (energy balance is balance_kcal)

Alt_E_balance_kcal2 <- Alt_E_balance_kcal %>%
  ungroup() %>%
  mutate(ID = factor(ID)) %>%
  group_by(ID, SABLE) %>%
  mutate(direction_EE = if_else((balance_kcal >0), "Positive", "Negative"))
  
#Join df for BW and for energy balance
Quality_control <- Alt_E_balance_kcal2  %>%
  left_join(
    BW_change_data %>% 
      select(ID, SABLE, first_bw, last_bw, bw_change, direction_bw),
    by = c("ID", "SABLE"))
  
  
#Graph energy balance (FI-TEE) vs change in BW for each ID at each SABLE
ggplot(Quality_control, aes(x = bw_change, y = balance_kcal, fill = GROUP)) +
  #geom_bar(stat = "summary", 
           #fun = "mean", 
           #position = position_dodge(width = 0.9)) +
  #geom_errorbar(stat = "summary", 
               # fun.data = mean_se, 
                #position = position_dodge(width = 0.9), 
                #width = 0.3) +
  # individual mouse points
  geom_point(
    color = "black",
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
    size = 2)+
  # mouse ID labels
  geom_text_repel(
    aes(label = ID, color = "black"),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9),
    size = 3,
    show.legend = FALSE) +
  theme_bw(base_size = 14) +
  facet_wrap(~SABLE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  scale_fill_manual(values = custom_colors) +
  labs(
    title="Energy balance (alt) vs. BW change",
    x = "Change in BW (g)",
    y = "Energy balance (kcal)",
    fill = "Treatment group")

  
  
#Use day rel to get the number of days in sable for each mouse for each sable measurement for
these days
#filter for only sable days and sum intake
FI_manual <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(INTAKE_GR >0 & INTAKE_GR <31) %>% #removes 1-29-25 measurements 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  #filter(!is.na(INTAKE_GR)) %>% 
  #mutate(INTAKE_GR = replace_na(INTAKE_GR, 0),) %>%
  mutate(
    day_rel = DATE - first(DATE),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(
    sable_idx = case_when(
      #Baseline (sable days 1-7)
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-12") ~ "SABLE_DAY_1", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-19") ~ "SABLE_DAY_1",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-26") ~ "SABLE_DAY_1",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-03") ~ "SABLE_DAY_1",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-13") ~ "SABLE_DAY_2", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-20") ~ "SABLE_DAY_2",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-27") ~ "SABLE_DAY_2",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-04") ~ "SABLE_DAY_2",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-14") ~ "SABLE_DAY_3", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-21") ~ "SABLE_DAY_3",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-28") ~ "SABLE_DAY_3",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-05") ~ "SABLE_DAY_3",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-15") ~ "SABLE_DAY_4", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-22") ~ "SABLE_DAY_4",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-29") ~ "SABLE_DAY_4",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-06") ~ "SABLE_DAY_4",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-16") ~ "SABLE_DAY_5", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-23") ~ "SABLE_DAY_5",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-11-30") ~ "SABLE_DAY_5",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-07") ~ "SABLE_DAY_5",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-17") ~ "SABLE_DAY_6", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-24") ~ "SABLE_DAY_6",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-01") ~ "SABLE_DAY_6",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-08") ~ "SABLE_DAY_6",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-18") ~ "SABLE_DAY_7", 
      ID %in% c(3713, 3714, 3716, 3717, 3718, 3719) & DATE == as.Date("2024-11-25") ~ "SABLE_DAY_7",
      ID %in% c(3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-02") ~ "SABLE_DAY_7",
      ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-09") ~ "SABLE_DAY_7",
      #Peak obesity (Sable days 8-11) 
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-05") ~ "SABLE_DAY_8", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_8",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-12") ~ "SABLE_DAY_8",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-02-28") ~ "SABLE_DAY_8",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-06") ~ "SABLE_DAY_9", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-09") ~ "SABLE_DAY_9",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-13") ~ "SABLE_DAY_9",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-01") ~ "SABLE_DAY_9",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-07") ~ "SABLE_DAY_10", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-10") ~ "SABLE_DAY_10",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-14") ~ "SABLE_DAY_10",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-02") ~ "SABLE_DAY_10",
      ID %in% c(3706, 3707, 3709, 3713, 3716, 3717) & DATE == as.Date("2025-02-08") ~ "SABLE_DAY_11", 
      ID %in% c(3708, 3710, 3714, 3719, 3720, 3721, 3726) & DATE == as.Date("2025-02-11") ~ "SABLE_DAY_11",
      ID %in% c(3722, 3723, 3724, 3725, 3728, 3729) & DATE == as.Date("2025-02-15") ~ "SABLE_DAY_11",
      ID %in% c(3711, 3718, 3727) & DATE == as.Date("2025-03-03") ~ "SABLE_DAY_11",
      #BW loss (Sable days 12-15)
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-08") ~ "SABLE_DAY_12", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-12") ~ "SABLE_DAY_12",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-16") ~ "SABLE_DAY_12",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-20") ~ "SABLE_DAY_12",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-24") ~ "SABLE_DAY_12",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-09") ~ "SABLE_DAY_13", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-13") ~ "SABLE_DAY_13",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-17") ~ "SABLE_DAY_13",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-21") ~ "SABLE_DAY_13",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-25") ~ "SABLE_DAY_13",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-10") ~ "SABLE_DAY_14", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-14") ~ "SABLE_DAY_14",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-18") ~ "SABLE_DAY_14",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-22") ~ "SABLE_DAY_14",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-26") ~ "SABLE_DAY_14",
      ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "SABLE_DAY_15", 
      ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "SABLE_DAY_15",
      ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "SABLE_DAY_15",
      ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "SABLE_DAY_15",
      ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "SABLE_DAY_15",
      #BW maintenance (Sable days 16-19) 
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-01") ~ "SABLE_DAY_16", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-05") ~ "SABLE_DAY_16",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-09") ~ "SABLE_DAY_16",
      ID %in% c(3729) & DATE == as.Date("2025-06-13") ~ "SABLE_DAY_16",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-02") ~ "SABLE_DAY_17", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-06") ~ "SABLE_DAY_17",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-10") ~ "SABLE_DAY_17",
      ID %in% c(3729) & DATE == as.Date("2025-06-14") ~ "SABLE_DAY_17",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-03") ~ "SABLE_DAY_18", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-07") ~ "SABLE_DAY_18",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-11") ~ "SABLE_DAY_18",
      ID %in% c(3729) & DATE == as.Date("2025-06-15") ~ "SABLE_DAY_18",
      ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713) & DATE == as.Date("2025-06-04") ~ "SABLE_DAY_19", 
      ID %in% c(3714, 3716, 3717, 3718, 3719, 3720, 3721) & DATE == as.Date("2025-06-08") ~ "SABLE_DAY_19",
      ID %in% c(3722, 3723, 3724, 3725, 3726, 3727, 3728) & DATE == as.Date("2025-06-12") ~ "SABLE_DAY_19",
      ID %in% c(3729) & DATE == as.Date("2025-06-16") ~ "SABLE_DAY_19",
      #BW regain (Sable days 20-23) 
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-05") ~ "SABLE_DAY_20", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_20",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-11") ~ "SABLE_DAY_20",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_20",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_20",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-18") ~ "SABLE_DAY_20",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_20",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-06") ~ "SABLE_DAY_21", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_21",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-12") ~ "SABLE_DAY_21",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_21",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_21",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-19") ~ "SABLE_DAY_21",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_21",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "SABLE_DAY_22", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_22",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-13") ~ "SABLE_DAY_22",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-15") ~ "SABLE_DAY_22",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_22",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-20") ~ "SABLE_DAY_22",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_22",
      ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-08") ~ "SABLE_DAY_23", 
      ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "SABLE_DAY_23",
      ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "SABLE_DAY_23",
      ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "SABLE_DAY_23",
      ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-17") ~ "SABLE_DAY_23",
      ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "SABLE_DAY_23",
      ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "SABLE_DAY_23",
      TRUE ~ NA_character_)) %>% 
  mutate(SABLE = case_when(
    sable_idx %in% c("SABLE_DAY_1", "SABLE_DAY_2", "SABLE_DAY_3", "SABLE_DAY_4", "SABLE_DAY_5", "SABLE_DAY_6", "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8", "SABLE_DAY_9", "SABLE_DAY_10", "SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12", "SABLE_DAY_13", "SABLE_DAY_14", "SABLE_DAY_15") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_16", "SABLE_DAY_17", "SABLE_DAY_18", "SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20", "SABLE_DAY_21", "SABLE_DAY_22", "SABLE_DAY_23") ~ "BW regain",
    TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>% 
  mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss","BW maintenance", "BW regain"))) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
  mutate(day_rel = as.numeric(day_rel))
  ungroup() 
  
#Now calcualte number of days during each sable period and intake during each sable period
  filtered_FI_manual <- FI_manual %>%
    group_by(ID, SABLE) %>%
    arrange(DATE) %>%
    summarise(
          start_sable= min(day_rel),
           end_sable = max(day_rel),
           total_days = end_sable - start_sable,
           sum_INTAKE_GR = sum(INTAKE_GR),
          INTAKE_per_day_GR= sum_INTAKE_GR/total_days,
          .groups = "drop") %>%
    mutate(KCAL_G = if_else(SABLE=="Baseline", 3.1, 3.82)) %>%
    mutate(INTAKE_per_day_kcal = INTAKE_per_day_GR*KCAL_G)
  
  
  summarise(
    first_bw = first(value),
    last_bw  = last(value),
    first_DateTime = first(DateTime),
    last_DateTime  = last(DateTime),
    bw_change = last_bw - first_bw,
    .groups = "drop")
  
  #what if I try removing the negative values for EE
  
  sable_TEE_modCS <- sable_dwn %>% 
    filter(COHORT %in% c(3, 4, 5)) %>%   
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light")) %>% 
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
    group_by(ID, complete_days) %>% 
    mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
    ungroup() %>% 
    filter(!ID %in% c(3715,3712)) %>% 
    filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>%
    filter(is_complete_day == 1) %>%
    filter(complete_days %in% c(1,2)) %>%
    mutate(
      GROUP = case_when(
        ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
        ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
      DRUG = case_when(
        ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
        ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
    mutate(SABLE = factor(SABLE,levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))) %>%
    filter(value >0)
  
  
  #---
  #Step 1a: for each hour, for each ID, find the mean of the EE measurements (in kcal/hr) taken each minute
  mod_TEE_avg_hr <- sable_TEE_modCS %>%
    ungroup() %>%
    group_by(SABLE, ID, hr, GROUP) %>%
    arrange(DateTime) %>%     # make sure rows are in time order
    summarise(avg_TEE_hr = mean(value))
  
  #Step 1b: graph hourly TEE from step 1 (units are kcal/hr) --> using "complete days" rather than "LM_complete_day"
  ggplot(mod_TEE_avg_hr, aes(x = hr, y = avg_TEE_hr, fill = GROUP)) +
    geom_bar(stat = "summary", 
             fun = "mean", 
             position = position_dodge(width = 0.9)) +
    geom_errorbar(stat = "summary", 
                  fun.data = mean_se, 
                  position = position_dodge(width = 0.9), 
                  width = 0.3) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 12)) +
    labs(
      title="Average energy expenditure each hour (kcal/hr)",
      x = "Time (Hour)",
      y = "TEE (kcal/hour)",
      fill = "Treatment group") +
    facet_wrap(~SABLE)