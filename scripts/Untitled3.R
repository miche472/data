#This script is a work in progress for RMR for C57 mice

#I copy and pasted the version of NZO RMR that I later realized had bugs, 
#I will try using that as a starting point for C57 RMR
  #I won't close this tab until I do that though

#Extrapolate resting metabolic rate (RMR) for C57 mice at each SABLE time point
#For each ID find 30 minute block with the lowest average TEE (kcal/hr)
#Extrapolate the average TEE during this 30min block (kcal/hr) to kcal/day by multiplying by 24

#Learning how to find 30min block with lowest avg TEE:
#1. Find RMR for ID=7861 at the BW regain time point
#2. Find RMR for all mice at the BW regain time point

#### 3. Goal: Calculate RMR for all mice at all SABLE time points ####
#4. Do multiple linear regression to calculate predictions for 
#RMR (kcal/day) adjusted for lean mass
#5. Determine if there is an EE gap based on RMR
#6. Graph RMR adjusted for lean mass vs SABLE time point

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

#TEE####
# build the summarized dataset
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT %in% c(2)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Off", "On")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3","SABLE_DAY_4") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7","SABLE_DAY_8") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12") ~ "BW maintenance", 
    sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain"
  )) %>%
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE, SEX) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days, SEX) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE, SEX) %>% 
  #summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
  filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE, SEX) %>% 
  #summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )) %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
    ) %>%
  #filter (ID == "7861") %>%
  #filter (SABLE == "BW regain")
  
  #### Make sable_TEE_adj by combining echomri data with sable_TEE_data ####
#Took this code from C57_Figure7 - TEE_correctedbyLean.R on 10-17-25
#Will combine sable_TEE_adj with the RMR df at the end of the code 

#echoMRI
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 2) %>% # Just C57 males and females
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "Restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG, SEX, DIET_FORMULA) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      Date == as.Date("2025-03-07") ~ "Peak obesity",
      Date == as.Date("2025-04-21") ~ "BW loss",
      Date == as.Date("2025-06-05") ~ "BW maintenance",
      Date %in% as.Date(c("2025-09-11", "2025-09-10","2025-09-05","2025-09-04",
                          "2025-09-02","2025-09-01","2025-08-28","2025-08-27")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) #%>%  # <-- optional (removes EchoMRI measurements that don't correspond to a sable time point)

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join Lean, Fat, and Weight info into TEE dataset
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean, Weight, Fat, SEX, DIET_FORMULA),
    by = c("ID", "SABLE")
  )

#### Find 30min with the lowest average TEE for one ID at one SABLE time point #### 
lowest_window <- sable_TEE_adj %>%
  arrange(DateTime) %>%
  mutate(
    # 30-minute moving average: current + previous 29 minutes
    avg_30min = slide_dbl(
      .x = value,
      .f = mean,
      .before = 29,      # include 30 total minutes
      .complete = TRUE)) %>%
  filter(!is.na(avg_30min)) %>%
  slice_min(avg_30min, n = 1) %>%
  mutate(
    window_end   = DateTime,
    window_start = DateTime - minutes(29)) %>%
  select(window_start, window_end, avg_30min)
lowest_window

#### Table of every 30 min time windows for ID=7861 during BW regain####
#In the code that makes the data frome "sable_TEE-data", filter for ID=7861 and SABLE=BW regain
all_windows <- sable_TEE_data %>%
  arrange(DateTime) %>%
  mutate(avg_30min = slide_dbl(
    .x = value,
    .f = mean,
    .before = 29,
    .complete = TRUE),
    window_end   = DateTime,
    window_start = DateTime - minutes(29)) %>%
  filter(!is.na(avg_30min)) %>%
  select(window_start, window_end, avg_30min) %>%
  arrange(avg_30min)

head(all_windows, 10) # Show the lowest few windows

#### Graphs all 30 minute windows for 3706 during BW regain####

# Step 1: Compute all 30-minute sliding averages
all_windows <- sable_TEE_data %>%
  arrange(DateTime) %>%
  mutate(avg_30min = slide_dbl(
    .x = value,
    .f = mean,
    .before = 29,
    .complete = TRUE),
    window_end   = DateTime,
    window_start = DateTime - minutes(29)) %>%
  filter(!is.na(avg_30min)) %>%
  select(window_start, window_end, avg_30min)

# Step 2: Identify the lowest-average 30-minute window
lowest_window <- all_windows %>% 
  slice_min(avg_30min, n = 1)

# Step 3: Plot with shaded area and labels
ggplot(sable_TEE_data, aes(x = DateTime, y = value)) +
  geom_line(color = "gray40") +
  # Highlight the lowest window
  geom_rect(
    data = lowest_window,
    aes(xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf),
    fill = "skyblue", alpha = 0.3, inherit.aes = FALSE) +
  # Optional: show the mean value as a point at the center of the window
  geom_point(
    data = lowest_window,
    aes(x = window_end - minutes(15), y = lowest_window$avg_30min),
    color = "blue", size = 3) +
  # Add annotation label
  geom_text(data = lowest_window, aes(x = window_end - minutes(15), y = lowest_window$avg_30min, 
                                      label = paste0(
                                        "Lowest 30-min avg = ", round(avg_30min, 2),
                                        "\nStart: ", format(window_start, "%H:%M"),
                                        "\nEnd: ", format(window_end, "%H:%M"))),vjust = -1, hjust = 0.5, size = 3.8, color = "black") +
  labs(title = "Value Over Time — Lowest 30-Minute Period Highlighted", x = "Time", y = "Value") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())

#### Adjust the code to include all mice at SABLE= BW regain ####
#Build df = sable_TEE_data without filtering for ID= 7861
# Compute sliding 30-minute averages for each mouse ID
lowest_windows_all <- sable_TEE_data %>%
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(avg_30min_value = slide_dbl(
    .x = value,
    .f = mean,
    .before = 29,       # previous 29 rows + current = 30 minutes
    .complete = TRUE),
    window_end_time   = DateTime,
    window_start_time = DateTime - minutes(29)) %>%
  filter(!is.na(avg_30min_value)) %>%
  ungroup()

# For each ID, find the lowest-average 30-minute period
lowest_windows_summary <- lowest_windows_all %>%
  group_by(ID) %>%
  slice_min(avg_30min_value, n = 1) %>%
  ungroup() %>%
  select(ID, window_start_time, window_end_time, avg_30min_value)

lowest_windows_summary # View results

#### Adjust code to include all IDs at all SABLE time points ####
#Create df= sable_TEE_data without filtering SABLE or ID
# Compute sliding 30-minute averages for each mouse *and* period
lowest_windows_all <- sable_TEE_data %>%
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
RMR_data <- lowest_windows_summary %>%
  rename(RMR_kcal_hr = avg_30min_value) %>%
  mutate(RMR_kcal_day = RMR_kcal_hr*24) %>%
  group_by(ID, SABLE)

#### Combine RMR_data with sable_TEE_adj ####
sable_TEE_RMR_adj <- sable_TEE_adj %>%
  left_join(
    RMR_data %>% select(ID, SABLE, window_start_time, window_end_time, RMR_kcal_day),
    by = c("ID", "SABLE")
  ) %>%
  group_by(ID, SABLE) %>%
  mutate(TEE_delt_RMR = tee - RMR_kcal_day)

#### Multiple linear regression model for RMR ####
#statistical model
model_RMR_lean <- lmer(RMR_kcal_day ~ SABLE * GROUP + Lean + (1 | ID), data = sable_TEE_RMR_adj)
summary(model_RMR_lean)

n_distinct(RMR_data$ID)

n_distinct(sable_TEE_RMR_adj$ID) #good we have 16 animals

emm_RMR <- emmeans(model_RMR_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_RMR_df <- as.data.frame(emm_RMR)

# Pairwise contrasts within each GROUP
contrasts_by_group <- contrast(emm_RMR, method = "pairwise", by = "GROUP")

# Convert to a data frame
contrasts_df <- as.data.frame(contrasts_by_group)

# Filter for restricted group and baseline vs BW maintenance
restricted_contrast <- contrasts_df %>%
  filter(GROUP == "restricted", contrast == "baseline - BW maintenance")

restricted_contrast

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
  geom_jitter(data = sable_TEE_RMR_adj, 
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
  labs(y = "Adjusted RMR (kcal/day)", x = "Time point",
       color = "Group",
       title = "NZO basal metabolic rate adjusted for lean mass") +
  format.plot +
  theme(legend.position = "top", plot.title = element_text(hjust=0.5), 
        axis.text = element_text( 
          color="black", 
          size=12))