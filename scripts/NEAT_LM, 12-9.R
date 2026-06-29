# Calculate NEAT

#Revised: 12-9-25
#Source data is a direct read in of the rds file created by data_proc
      #name of this is sable_downsampled_data.rds
#Script directly calculates NEAT by adding up the energy expenditure recorded 
      #during minutes when the mouse moved (indicated by an increase in AllMeters)
#Goal: calculate NEAT so that I can add it to the calculated RMR and see how close 
      #it is to the calculated TEE. This is following up from the meeting with Jen Teske

#### All SABLE time points for all NZO mice ####
#------- Get only the kcal_hr parameter --> value column has kcal_hr data only
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

filter_EE <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#------- Get only the AllMeters parameter --> value column has AllMeters data
filter_loc <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  #For NZO_Figure7b-RRM_correctedbyLean I used both complete days 1 and 2...complete_days %in% c(1,2)
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#---
#In df filter_loc use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE use mutate to make a column called kcal_hr_ using the value column data

filter_locom <- filter_loc %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_energy <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)


#Join filter_locom and filter_energy into a df called filter_locom_energy

# Left join Lean, Fat, and Weight info into TEE data set
filter_locom_energy <- filter_locom %>%
  left_join(
    filter_energy %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Calculate daily NEAT for each mouse ID at each SABLE time point
ID_NEAT <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
  summarise(NEAT_teske = sum(EE_per_min), .groups="drop")

#Calculate average daily NEAT within each GROUP at each SABLE time point
GROUP_NEAT <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
  summarise(NEAT_teske = sum(EE_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  #summarise(GROUP_EE_teske = mean(EE_teske), .groups="drop")
  summarise(GROUP_NEAT = mean(NEAT_teske), 
          SD_GROUP_NEAT = sd(NEAT_teske),
          SE_GROUP_NEAT = sd(NEAT_teske)/sqrt(n()))

#### End of NEAT calculation for all IDs at all SABLE time points ####

####Use the same approach to calculate TEE (rather than NEAT) ####

#Calculate daily TEE for each mouse ID at each SABLE time point
ID_TEE <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE, ID, GROUP) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
  summarise(EE_teske = sum(EE_per_min), .groups="drop")

#Calculate average daily TEE within each GROUP at each SABLE time point
GROUP_TEE <- filter_locom_energy %>%
  ungroup() %>%
  arrange(DateTime) %>%     # make sure rows are in time order...not totally necessary but good for consistency with NEAT calc.
  group_by(SABLE, ID, GROUP) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
  summarise(TEE_teske = sum(EE_per_min), .groups="drop") %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  #summarise(GROUP_EE_teske = mean(EE_teske), .groups="drop")
  summarise(GROUP_TEE = mean(TEE_teske), 
            SD_GROUP_TEE = sd(TEE_teske),
            SE_GROUP_TEE = sd(TEE_teske)/sqrt(n()))

#Left join to combine TEE and NEAT for individual IDs (join ID_Neat with ID_TEE)
TEE_NEAT <- ID_TEE %>%
  left_join(
    ID_NEAT %>% 
      select(NEAT_teske, ID, GROUP, SABLE),
    by = c("ID", "GROUP", "SABLE"))

#### Body composition (process echomri.csv for attachment) ####
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    SABLE = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_)) %>% 
  filter(!is.na(SABLE)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make SABLE an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain")))

# Left join body comp to TEE_NEAT dataframe to make TEE_NEAT_adj
TEE_NEAT_adj <- TEE_NEAT %>%
  left_join(
    echoMRI_data %>% 
      select(ID, SABLE, Lean, Weight, Fat),
    by = c("ID", "SABLE"))

#Calculate average daily NEAT and TEE within each GROUP at each SABLE time point
GROUP_TEE_NEAT <- TEE_NEAT_adj %>%
  ungroup() %>%
  group_by(SABLE, GROUP) %>%
  summarise(GROUP_TEE = mean(EE_teske), 
            SD_GROUP_TEE = sd(EE_teske),
            SE_GROUP_TEE = sd(EE_teske)/sqrt(n()),
            GROUP_NEAT = mean(NEAT_teske), 
            SD_GROUP_NEAT = sd(NEAT_teske),
            SE_GROUP_NEAT = sd(NEAT_teske)/sqrt(n()))


#Confirm that data frame has correct number of mice
n_distinct(TEE_NEAT_adj$ID) #good we have 16 animals

#### Build linear mixed model for TEE (adjusted by lean) ####
model_TEE_NEAT_lean <- lmer(EE_teske ~ SABLE * GROUP + Lean + (1 | ID), data = TEE_NEAT_adj)
summary(model_TEE_NEAT_lean)

#### Calculate estimated marginal means ####
emm_TEE_NEAT_lean <- emmeans(model_TEE_NEAT_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_NEAT_lean_df <- as.data.frame(emm_TEE_NEAT_lean)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_NEAT_lean <- contrast(emm_TEE_NEAT_lean, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_TEE_NEAT_lean_df <- as.data.frame(contrasts_by_group_TEE_NEAT_lean)

#Filter for significant contrasts
Sig_contrasts_TEE_NEAT_lean <- contrasts_TEE_NEAT_lean_df %>%
  filter(p.value <=0.05)
Sig_contrasts_TEE_NEAT_lean

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_NEAT_lean <- contrast(emm_TEE_NEAT_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_TEE_NEAT_lean_df <- as.data.frame(contrasts_by_SABLE_TEE_NEAT_lean)
contrasts_SABLE_TEE_NEAT_lean_df

#### Bar plot - Graph predicted TEE adjusted for lean mass ####
#Format plot (optional)
#scale_color_manual <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))
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

barplot_emm_TEE_NEAT_lean <- emm_TEE_NEAT_lean_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
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
    title = "Total energy expenditure (AllMeters method) adjusted for FFM",
    y = "Adjusted TEE (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_TEE_NEAT_lean


#Try adjusting by BW from manual BW measurement that occurred on the relavent sable day (get from Sable_BW)

Manual_BW <- Sable_BW %>%
  mutate(ID = as.factor(ID)) %>%
filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
rename(Manual_BW = BW_sable)

#left join onto TEE_NEAT_adj
TEE_NEAT_adj_manualBW <- TEE_NEAT_adj %>%
  left_join(
    Manual_BW %>% 
      select(ID, SABLE, Manual_BW),
    by = c("ID", "SABLE"))

#### Build linear mixed model for TEE (adjusted by manually measured BW) ####
model_TEE_NEAT_BW <- lmer(EE_teske ~ SABLE * GROUP + Manual_BW + (1 | ID), data = TEE_NEAT_adj_manualBW)
summary(model_TEE_NEAT_BW)

#### Calculate estimated marginal means ####
emm_TEE_NEAT_BW <- emmeans(model_TEE_NEAT_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_NEAT_BW_df <- as.data.frame(emm_TEE_NEAT_BW)

# Pairwise contrasts within each GROUP
contrasts_by_group_TEE_NEAT_BW <- contrast(emm_TEE_NEAT_BW, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_TEE_NEAT_BW_df <- as.data.frame(contrasts_by_group_TEE_NEAT_BW)

#Filter for significant contrasts
Sig_contrasts_TEE_NEAT_BW <- contrasts_TEE_NEAT_BW_df %>%
  filter(p.value <=0.05)
Sig_contrasts_TEE_NEAT_BW

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_TEE_NEAT_BW <- contrast(emm_TEE_NEAT_BW, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_TEE_NEAT_BW_df <- as.data.frame(contrasts_by_SABLE_TEE_NEAT_BW)
contrasts_SABLE_TEE_NEAT_BW_df

#### left off here ####
#calculate adjusted NEAT for manual BW and lean
#determine how to do ANCOVA for adjustment 

#---------------------------------------------------------------
# I used this when I was figuring out how to calculate NEAT using Jen Teske's approach #
#### NEAT energy expenditure for 3706 at all SABLE time points ####
#------- Get only the kcal_hr parameter --> value column has kcal_hr data only
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

filter_EE <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#------- Get only the AllMeters parameter --> value column has AllMeters data
filter_loc <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  #For NZO_Figure7b-RRM_correctedbyLean I used both complete days 1 and 2...complete_days %in% c(1,2)
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#---Rename variables to allow for joining of the two data frames...locomotion and EE for 3706

filter_loc_3706 <- filter_loc %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) %>%
  filter(ID=="3706")

filter_EE_3706 <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value) %>%
  filter(ID=="3706")
  

#Join filter_loc_3706 and filter_EE_3706 into a df called filter_loc_EE_3706allsable

# Left join Lean, Fat, and Weight info into TEE data set
filter_loc_EE_3706allsable <- filter_loc_3706 %>%
  left_join(
    filter_EE_3706 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#Do filtering, grouping, summarizing, etc. --> create df NEAT_3706_allsable
  #the "summarise" step gives physical activity EE (i.e. NEAT) for all time points for 3706

NEAT_3706_allsable <- filter_loc_EE_3706allsable %>%
  ungroup() %>%
  #filter(SABLE == "BW loss") %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  group_by(SABLE) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
summarise(EE_teske = sum(EE_per_min), .groups="drop")

#### End of code for 3706 at all sable time points ####


#### Code for 3706 at BW loss SABLE time point ####

#------- Get only the kcal_hr parameter --> value column has kcal_hr data only
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

filter_EE <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#------- Get only the AllMeters parameter --> value column has AllMeters data
filter_loc <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
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
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  # remove dead mice, keep both complete days, remove mice with cage issues
  #For NZO_Figure7b-RRM_correctedbyLean I used both complete days 1 and 2...complete_days %in% c(1,2)
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days==2) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#---
#In df filter_loc use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE use mutate to make a column called kcal_hr_ using the value column data

filter_loc_3706 <- filter_loc %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) %>%
  filter(ID=="3706")

filter_EE_3706 <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value) %>%
  filter(ID=="3706")


#Join filter_loc and filter_EE into a df called filter_loc_EE

# Left join Lean, Fat, and Weight info into TEE data set
filter_loc_EE_3706 <- filter_loc_3706 %>%
  left_join(
    filter_EE_3706 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

NEAT_3706 <- filter_loc_EE_3706 %>%
  ungroup() %>%
  group_by(SABLE) %>%
  #rename(SABLE = SABLE.x) %>% #This chow and LFD identified, as applicable 
  filter(SABLE == "BW loss") %>%
  arrange(DateTime) %>%     # make sure rows are in time order
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing=0)) %>%
  filter(move==1) %>%
  ungroup() %>%
  #group_by(DateTime) %>%
  mutate(EE_per_min= (Kcal_Hr/60)) %>%
  summarise(EE_teske = sum(EE_per_min), .groups="drop")

#### end of code for just 3706 at BW loss ####

