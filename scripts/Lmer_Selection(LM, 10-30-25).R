#LM (revised 10-29-25)

#Model to identify EE gap in NZO mice
    #First, repeat the calculations that identify RMR for each mouse and create 
    #a df with RMR, body composition, and TEE for each mouse ID
#Purpoe of this script is to select a linear mixed model

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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
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



#### Create a mixed model to check for EE gap ####

####Determine if I should use interaction term for Sable and GROUP or just include 
    #both in the mixed effects model
#make both models
model_int <- lmer(RMR_kcal_day ~ GROUP * SABLE + Lean + (1 | ID), data = sable_TEE_adj_RMR)
model_add <- lmer(RMR_kcal_day ~ GROUP + SABLE + Lean + (1 | ID), data = sable_TEE_adj_RMR)

#Compare the two models using 
#1. likelihood ratio test (tells you the relative utility of each model)
#2. AIC (tells you if each model is good)
anova(model_add, model_int)
#The p-value for the interaction term is 2.675e-8, which is <0.0001 suggesting
  #that including interaction significantly improves the model

####Compare models visually (measured, predicted, emmeans)####

#Step 1a: create predictions for each model
sable_TEE_adj_RMR_check <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  mutate(
    pred_int = predict(model_int, re.form = NA), #makes new column
    pred_add = predict(model_add, re.form = NA) #new column
  )

pred_df <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  mutate(
    pred_int = predict(model_int, re.form = NA), #makes new column
    pred_add = predict(model_add, re.form = NA) #makes new column
  )

#Step 1b: compute estimated marginal means
emm_int <- emmeans(model_int, ~ SABLE * GROUP)
emm_add <- emmeans(model_add, ~ SABLE * GROUP)

emm_int_df <- as.data.frame(emm_int) %>%
  mutate(model = "Interaction (GROUP*SABLE)")

emm_add_df <- as.data.frame(emm_add) %>%
  mutate(model = "Additive (GROUP + SABLE)")

emm_combined <- bind_rows(emm_int_df, emm_add_df)

#Step 2: Combine everything for visualization
pred_long <- pred_df %>%
  pivot_longer(
    cols = c(pred_int, pred_add),
    names_to = "model",
    values_to = "pred_RMR"
  ) %>%
  mutate(model = recode(model,
                        pred_int = "Interaction (GROUP*SABLE)",
                        pred_add = "Additive (GROUP + SABLE)"))

#Step 3: make graphs to visualize the two models
#Graph of GROUP + SABLE (additive mixed model) 
    #-->measured RMR (each ID), prediceted RMR (each ID), emmeans RMR (each GROUP and SABLE)
#Graph of GROUP*SABLE (interaction mixed model) 
#-->measured RMR (each ID), prediceted RMR (each ID), emmeans RMR (each GROUP and SABLE)
ggplot() +
  # Observed points
  geom_point(
    data = pred_long,
    aes(x = SABLE, y = RMR_kcal_day, color = GROUP),
    alpha = 0.3) +
  # Individual predicted lines (faint)
  geom_line(
    data = pred_long,
    aes(x = SABLE, y = pred_RMR, color = GROUP, group = ID),
    alpha = 0.3) +
  # Group mean predicted lines (thick)
  geom_line(
    data = emm_combined,
    aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
    size = 1.2) +
  # Error bars for EMM 95% CI
  # Note that this uses CI rather than SEM for error bars. I could change this
  geom_errorbar(
    data = emm_combined,
    aes(x = SABLE, ymin = lower.CL, ymax = upper.CL, color = GROUP),
    width = 0.2) +
  facet_wrap(~model) +
  theme_classic() +
  labs(
    title = "Predicted RMR: Additive vs Interaction Mixed Models",
    x = "SABLE timepoint",
    y = "RMR (kcal/day)",
    color = "Diet group")

#### Multiple linear regression model for RMR (adjusted for lean mass) ####
    #Based on the chunk of code above I conclude that I should use the model
    #with interaction between SABLE*GROUP (rather than additive)

#Build statistical model (to get RMR that is adjusted for Lean, multiple measures, etc.)
model_RMR_lean <- lmer(RMR_kcal_day ~ SABLE * GROUP + Lean + (1 | ID), data = sable_TEE_adj_RMR)
summary(model_RMR_lean)

n_distinct(sable_TEE_adj_RMR$ID) #Confirm that we have 16 animals

emm_RMR <- emmeans(model_RMR_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_df <- as.data.frame(emm_RMR)

# Pairwise contrasts within each GROUP
contrasts_by_group <- contrast(emm_RMR, method = "pairwise", by = "GROUP")

# Convert to a data frame
contrasts_df <- as.data.frame(contrasts_by_group)

#Significant contrasts
contrasts_df <- contrasts_df%>%
  mutate(Significant = p.value<0.05)
contrasts_df

# Filter for GROUP == restricted and SABLE=baseline vs SABLE=BW maintenance
Restricted_contrast <- contrasts_df %>%
  filter(GROUP == "Restricted", contrast == "Baseline - BW maintenance")

Restricted_contrast

# Filter for GROUP == ad libitum and SABLE=baseline vs SABLE=BW maintenance
adlib_contrast <- contrasts_df %>%
  filter(GROUP == "ad lib", contrast == "Baseline - BW maintenance")

adlib_contrast

