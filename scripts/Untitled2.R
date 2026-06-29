#This script is a work in progress for the Kevin Hall graph

#Resting metabolic rate
#RMR = 30 minute period with the lowest average measure for TEE (kcal/hr)
#Extrapolate the average TEE (kcal/hr) to kcal/day by multiplying by 24
#Learning how to do this procedure:
#1. Find RMR for ID=3706 at the BW regain time point
#2. Find RMR for all mice at the BW regain time point
#3. Find RMR for all mice at all SABLE time points
#4. Do multiple linear regression to calculate predictions for 
#RMR (kcal/day) adjusted for lean mass
#5. Determine if there is an EE gap based on RMR
#6. Graph adjusted RMR 

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
                     "SABLE_DAY_7") ~ "baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
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
  
  group_by(ID, SABLE) %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_RMR_data <- sable_RMR_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )
#filter (ID == "3706") %>%
#filter (SABLE == "BW regain")

#### Code including all IDs at all SABLE time points ####
#This step requires a df that has minute data for TEE rather than pre-calculated average daily TEE

# Compute sliding 30-minute averages for each mouse *and* period
lowest_windows_all <- sable_RMR_data %>%
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(avg_30min_value = slide_dbl(
    .x = value,
    .f = mean,
    .before = 29,       # previous 29 rows + current = 30 minutes
    .complete = TRUE),
    window_end_time   = DateTime,
    window_start_time = DateTime - minutes(29)) %>%
  filter(!is.na(avg_30min_value)) %>%
  ungroup()

# For each ID and SABLE, find the lowest 30-minute average
lowest_windows_summary <- lowest_windows_all %>%
  group_by(ID, SABLE) %>%
  slice_min(avg_30min_value, n = 1) %>%
  ungroup() %>%
  select(ID, SABLE, window_start_time, window_end_time, avg_30min_value)

# View summary
lowest_windows_summary

#### Change RMR units from kcal_hr to kcal_day to match tee units####
lowest_windows_summary <- lowest_windows_summary %>%
  rename(RMR_kcal_hr = avg_30min_value) %>%
  mutate(RMR_kcal_day = RMR_kcal_hr*24) %>%
  group_by(ID, SABLE)


##### Process sable_dwn into sable_TEE_data to get Avg daily TEE (tee) for each mouse at each time point####
# build the summarized dataset 
#version with creation of tee 
#Join it with echo data to create sable_tee_adj 
#Then join sable_tee_adj with lowest_windows_summary (i.e. df with RMR) to create sable_TEE_adj_RMR
#use this compiled code to do tee-BMR and for linear regression models and graphing)
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
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
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>%
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )


#### Attach echoMRI_data to sable_TEE_data --> name new df as sable_TEE_adj ####
#echo info is from NZO_Figure7 - TEE_correctedbyLean (rev. LM).R on (LM accessed on 10-16-25)

#Process echoMRI info for NZO mice
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      Date == as.Date("2025-02-20") ~ "peak obesity",
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
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join Lean, Fat, and Weight info into TEE dataset
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean, Weight, Fat),
    by = c("ID", "SABLE")
  )

##### Combine lowest_windows_summary with sable_TEE_adj ####
sable_TEE_adj_RMR <- sable_TEE_adj %>%
  left_join(
    lowest_windows_summary %>% select(ID, SABLE, window_start_time, window_end_time, RMR_kcal_day),
    by = c("ID", "SABLE")
  ) %>%
  group_by(ID, SABLE) %>%
  mutate(TEE_minus_RMR = tee - RMR_kcal_day)

#### Multiple linear regression model for RMR (FFM, FM) ####
#statistical model
plot_prep <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted") %>%
  filter(SABLE == "baseline")
model_RMR_regline <- lm(RMR_kcal_day ~ Lean + Fat, data = plot_prep)
summary(model_RMR_regline)

#n_distinct(sable_TEE_adj_RMR$ID) #good we have 16 animals

#emm_RMR <- emmeans(model_RMR_lean, ~ SABLE * GROUP, cov.reduce = mean)
#emm_RMR_df <- as.data.frame(emm_RMR)

# Pairwise contrasts within each GROUP
#contrasts_by_group <- contrast(emm_RMR, method = "pairwise", by = "GROUP")

# Convert to a data frame
#contrasts_df <- as.data.frame(contrasts_by_group)

# Filter for restricted group and baseline vs BW maintenance
#restricted_contrast <- contrasts_df %>%
 # filter(GROUP == "restricted", contrast == "baseline - BW maintenance")

#restricted_contrast

# ==========================================================
# STEP 1: Fit the linear regression model
# ==========================================================
# Predicting Resting Metabolic Rate (RMR) from Fat and Lean
model <- lm(RMR_kcal_day ~ Fat + Lean + (1 | ID), data = plot_prep)
summary(model)

# ==========================================================
# STEP 2: Get tidy outputs of the regression results
# ==========================================================
library(broom)

# Coefficients table (estimate, std error, p-value)
tidy(model)

# Model fit statistics (R², AIC, etc.)
glance(model)

# ==========================================================
# STEP 3: Add predicted RMR values to your dataset
# ==========================================================
library(dplyr)

plot_prep2 <- plot_prep %>%
  mutate(predicted_RMR = predict(model, newdata = plot_prep))

# ==========================================================
# STEP 4: Plot predicted vs actual RMR to check model fit
# ==========================================================
library(ggplot2)

ggplot(plot_prep, aes(x = predicted_RMR, y = RMR)) +
  geom_point(color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Predicted vs Actual Resting Metabolic Rate",
    x = "Predicted RMR",
    y = "Actual RMR"
  ) +
  theme_minimal()

# ==========================================================
# STEP 5: Extract regression equation coefficients
# ==========================================================
coef(model)

# Example interpretation (values will depend on your data):
# Predicted RMR = Intercept + (b1 * Fat) + (b2 * Lean)

# ==========================================================
# STEP 6: Plot adjusted (predicted) RMR vs Lean, colored by Fat
# ==========================================================
ggplot(plot_prep, aes(x = Lean, y = predicted_RMR, color = Fat)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Adjusted Resting Metabolic Rate vs Lean Mass",
    x = "Lean Mass",
    y = "Adjusted (Predicted) RMR",
    color = "Fat Mass"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#statistical model
model1 <- lmer(RMR_kcal_day ~ Fat + Lean + (1 | ID), data = plot_prep)
summary(model1)




emm_TEE <- emmeans(model_TEE_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_df <- as.data.frame(emm_TEE)
#------------------------------------------------------
#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  
  # remove background grid lines only
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  # keep axis lines
  axis.line = element_line(color = "black")
)

####Graph predicted RMR adjusted for lean mass####
ggplot() +
  geom_jitter(data = sable_TEE_adj_RMR, 
              aes(x = SABLE, y = RMR_kcal_day, color = GROUP),
              width = 0.2, alpha = 0.4, size = 2) +
  geom_point(data = emm_RMR_df,
             aes(x = SABLE, y = emmean, color = GROUP),
             position = position_dodge(0.2), size = 4) +
  geom_line(data = emm_RMR_df,
            aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
            position = position_dodge(0.2), linewidth = 1.5) + 
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  geom_errorbar(data = emm_RMR_df,
                aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                width = 0.15, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(y = "Adjusted REE (kcal/day)", x = "Time point",
       color = "Group",
       title = "NZO resting energy expenditure adjusted for lean mass") +
  #format.plot +
  theme(legend.position = "top", plot.title = element_text(hjust=0.5), 
        axis.text = element_text( 
          color="black", 
          size=12))