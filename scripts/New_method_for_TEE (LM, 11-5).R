#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels
library(slider)
library(lubridate)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#### Resting metabolic rate (RMR): Identify 30min with lowest avg. TEE ####
#Use the df created in this chunk (sable_RMR_data) for the code that identifies the 30min window
#This code is basically creating sable_TEE_data, but the steps that calculate avg tee are deleted

sable_RMR_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  
  # remove dead mice, keep both complete days, remove mice with cage issues
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  ####added today --calc tee ####
#try to get the average kcal_hr across completeday=1 and 2 for mouse with ID 3706
sable_EE_data_mean <- sable_RMR_data %>% #will need to change sable_RMR_data to use sable_dwn as the source, but have a diff name for the df created
  
ungroup() %>%
  group_by(ID,SABLE)%>%
  #summarise(mean_tee = mean(value)) %>%
  mutate(mean_tee = mean(value))
  
  group_by(ID, SABLE) %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_RMR_data <- sable_RMR_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )

###############-----Left off here...start of chat attempt--------------#################
library(dplyr)
library(tidyr)
library(lubridate)

# --- Helper function ---
zt_time <- function(hr){
  if_else(hr >= 20 & hr <= 23, hr - 20, hr + 4)
}

# --- Load data ---
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds")

# --- Process and compute daily average TEE ---
sable_daily_TEE2 <- sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%
  # Map sable_idx to timepoint labels
  mutate(
    SABLE = case_when(
      sable_idx %in% paste0("SABLE_DAY_", 1:7) ~ "Baseline",
      sable_idx %in% paste0("SABLE_DAY_", 8:11) ~ "Peak obesity",
      sable_idx %in% paste0("SABLE_DAY_", 12:15) ~ "BW loss",
      sable_idx %in% paste0("SABLE_DAY_", 16:19) ~ "BW maintenance",
      sable_idx %in% paste0("SABLE_DAY_", 20:23) ~ "BW regain",
      TRUE ~ NA_character_
    ),
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")
  ) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>%
  
  # Compute ZT time and identify complete days
  group_by(ID, SABLE) %>%
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_day = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>%
  ungroup() %>%
  
  # Keep only full days (ZT0–ZT23)
  group_by(ID, SABLE, complete_day) %>%
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>%
  ungroup() %>%
  filter(is_complete_day == 1, complete_day %in% c(1,2)) %>%
  
  # Remove problematic mice
  filter(!ID %in% c(3709, 3712, 3715, 3717, 3718, 3723, 3724, 3725)) %>%
  
  # --- STEP 1: Compute mean TEE per complete day ---
  group_by(ID, SABLE, complete_day) %>%
  summarise(mean_TEE_kcal_hr = mean(value, na.rm = TRUE), .groups = "drop") %>%
  
  # --- STEP 2: Average across complete days to get one value per ID × SABLE ---
  group_by(ID, SABLE) %>%
  summarise(avg_daily_TEE_kcal_hr = mean(mean_TEE_kcal_hr, na.rm = TRUE), .groups = "drop") %>%
  
  # --- Add GROUP and DRUG info ---
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    ),
    SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))
  )

#####--####
# --- Libraries ---
library(dplyr)
library(tidyr)
library(lubridate)

# --- Helper function ---
zt_time <- function(hr){
  if_else(hr >= 20 & hr <= 23, hr - 20, hr + 4)
}

# --- Load data ---
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds")

# --- Process and compute TEE per day and per time point ---
sable_daily_TEE <- sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%
  
  # Map sable_idx to timepoint labels
  mutate(
    SABLE = case_when(
      sable_idx %in% paste0("SABLE_DAY_", 1:7) ~ "Baseline",
      sable_idx %in% paste0("SABLE_DAY_", 8:11) ~ "Peak obesity",
      sable_idx %in% paste0("SABLE_DAY_", 12:15) ~ "BW loss",
      sable_idx %in% paste0("SABLE_DAY_", 16:19) ~ "BW maintenance",
      sable_idx %in% paste0("SABLE_DAY_", 20:23) ~ "BW regain",
      TRUE ~ NA_character_
    ),
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")
  ) %>%
  
  # Keep only TEE values
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>%
  
  # Compute ZT time and identify complete days
  group_by(ID, SABLE) %>%
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_day = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>%
  ungroup() %>%
  
  # Keep only full days (ZT0–ZT23)
  group_by(ID, SABLE, complete_day) %>%
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>%
  ungroup() %>%
  filter(is_complete_day == 1) %>%
  
  # Remove problematic mice
  filter(!ID %in% c(3709, 3712, 3715, 3717, 3718, 3723, 3724, 3725)) %>%
  
  # --- STEP 1: Compute mean kcal/hr and kcal/day per complete day ---
  group_by(ID, SABLE, complete_day) %>%
  summarise(
    tee_hr = mean(value, na.rm = TRUE),
    tee_day = tee_hr * 24,  # convert kcal/hr → kcal/day
    .groups = "drop"
  ) %>%
  
  # --- STEP 2: Average across complete days to get one value per ID × SABLE ---
  group_by(ID, SABLE) %>%
  summarise(
    avg_tee_hr = mean(tee_hr, na.rm = TRUE),
    avg_tee_day = mean(tee_day, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # --- Add GROUP and DRUG info ---
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    ),
    SABLE = factor(SABLE, levels = c("Baseline", "Peak obesity", "BW loss", "BW maintenance", "BW regain"))
  )

#echoMRI
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join lean mass info into TEE dataset
sable_daily_TEE_adj <- sable_daily_TEE %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean, Fat, Weight),
    by = c("ID", "SABLE")
  )


#statistical model
model_daily_lean <- lmer(avg_tee_day ~ SABLE * GROUP + Lean + (1 | ID), data = sable_daily_TEE_adj)
summary(model_daily_lean)

n_distinct(sable_daily_TEE_adj$ID) #good we have 16 animals


emm_daily_lean <- emmeans(model_daily_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_daily_lean_df <- as.data.frame(emm_daily_lean)



# Pairwise contrasts within each GROUP
contrasts_by_group_daily_lean <- contrast(emm_daily_lean, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_group_daily_lean_df <- as.data.frame(contrasts_by_group_daily_lean)

# Filter for restricted group and baseline vs BW maintenance
restricted_contrast_group_daily_lean <- contrasts_group_daily_lean_df %>%
  filter(GROUP == "Restricted", contrast == "Baseline - BW regain")
restricted_contrast_group_daily_lean

# Pairwise contrasts within each GROUP
contrasts_by_SABLE_daily_lean <- contrast(emm_daily_lean, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_daily_lean_df <- as.data.frame(contrasts_by_SABLE_daily_lean)
contrasts_SABLE_daily_lean_df

# time point comparisons
Sig_SABLE_daily_lean <- contrasts_SABLE_daily_lean_df %>%
  filter(p.value <= 0.05) 
Sig_SABLE_daily_lean

#Plot of TEE (adj. for lean mass) vs Sable phase--plot of emmeans
ggplot() +
  geom_jitter(data = sable_daily_TEE_adj, 
              aes(x = SABLE, y = avg_tee_day, color = GROUP),
              width = 0.2, alpha = 0.4, size = 2) +
  geom_point(data = emm_daily_lean_df,
             aes(x = SABLE, y = emmean, color = GROUP),
             position = position_dodge(0.2), size = 3) +
  geom_line(data = emm_daily_lean_df,
            aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
            position = position_dodge(0.2), linewidth = 1) +
  geom_errorbar(data = emm_daily_lean_df,
                aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                width = 0.1, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(y = "TEE (adjusted for Lean mass)", x = "SABLE phase",
       color = "Group",
       title = "Ne method: TEE across SABLE phases (adjusted for Lean mass)") +
  theme(legend.position = "top")

# Collapse adjusted data the same way you did before for non lean mass correction
TEE_adj_plotdata <- emm_TEE_df %>%
  mutate(
    PlotGroup = case_when(
      SABLE %in% c("baseline", "peak obesity") ~ "all",          # collapse all
      SABLE %in% c("BW loss", "BW maintenance") ~ GROUP,         # separate by GROUP
      SABLE == "BW regain" ~ paste(GROUP, "RTIOXA_47", sep = "_") # keep consistent labeling
    )
  )

