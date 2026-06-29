#Goal of script: Graph and do stats for NZO adiposity index

# Revised: 11-5-25

#Need to run script to make sable_TEE_RMR_adj to use as source data for adiposity index etc.
#Located at the bottom of this script if it is not already in the Environment

# Adiposity index: ANOVA, T-test, linear mixed model and emmeans, bar graph
# Use sable_tee_adj_RMR (uses echoMRI_data adiposity index and only the mice with Sable measures)

####Libraries####
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

####----- Graphs of adiposity index at the 5 time points (n=16 NZO mice) -----####
#****Add bars to show statistical differences in PowerPoint

#Format plot
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#3498DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

#####----- Barplot - NZO adiposity index (measured values from echoMRI_data/sable_TEE_adj_RMR) ####
AI_barplot_measured <- sable_TEE_adj_RMR %>%
  ggplot(aes(x = SABLE, y = adiposity_index, fill = GROUP)) + 
  stat_summary( # mean bars
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  stat_summary( # error bars (mean ± SE)
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(  # individual data points
    aes(color = GROUP), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  labs(title = "NZO adiposity index", 
       y = "Adiposity index (Fat/Lean)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  scale_y_continuous( # set y-axis breaks every 0.2 units
    breaks = seq(0, max(sable_TEE_adj_RMR$adiposity_index, na.rm = TRUE), by = 0.2)) +
  format.plot +
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color="black", size=14),
    plot.title = element_text(hjust = 0.5, face = "bold"))
AI_barplot_measured

##### Scatterplot -NZO adiposity index (measured values from echoMRI_data/sable_TEE_adj_RMR)####

AI_scatterplot_measured <- sable_TEE_adj_RMR %>%
  ggplot(aes(x = SABLE, y = adiposity_index, color = GROUP)) + 
  #geom_point(
  #position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
  #alpha = 0.7,size = 2) +
  geom_line(aes(group = ID), alpha = 0.3, size = 0.8) + # connect points from the same ID
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = GROUP),
    #position = position_dodge(width = 0.6), shape = 18,
    size = 3) +
  stat_summary( # line connecting group means
    fun = mean,
    geom = "line",
    #position = position_dodge(width = 0.6),
    aes(group = GROUP),
    size = 1) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    aes(group = GROUP),
    #position = position_dodge(width = 0.6),
    width = 0.2) +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal() +
  labs(
    title = "NZO adiposity index", 
    y = "Adiposity index (Fat/Lean)", 
    x = "Time point", 
    color = "Treatment group") +
  scale_y_continuous(
    limits = c(0, NA),   # start at 0, upper limit adjusts automatically
    expand = c(0, 0),    # optional: removes extra padding above/below
    breaks = seq(0, max(sable_TEE_adj_RMR$adiposity_index, na.rm = TRUE), by = 0.2)) +
  #scale_y_continuous(breaks = seq(0, max(sable_TEE_adj_RMR$adiposity_index, na.rm = TRUE), by = 0.2)) +
  format.plot +
  theme(legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color="black", size=14),
    plot.title = element_text(hjust = 0.5, face = "bold"))
AI_scatterplot_measured

####-----T-test, repeated measures ANOVA, post hoc for ANOVA -----####
#T-test: pairwise comparison between GROUP (ad lib and restricted) at each STATUS
ttest_results_AI <- sable_TEE_adj_RMR %>%
  group_by(SABLE) %>%
  t_test(adiposity_index ~ GROUP, var.equal = TRUE) %>%   # or var.equal = FALSE if not assumed
  adjust_pvalue(method = "bonferroni") %>%   # optional multiple-comparison correction
  add_significance("p.adj")                  # adds stars based on adjusted p-values
ttest_results_AI

#Repeated measures ANOVA: comparison of adiposity index between groups across the 5 time points
# Perform repeated-measures ANOVA separately for each GROUP
anova_results_AI <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  group_by(GROUP) %>%
  anova_test(dv = adiposity_index, wid = ID, within = SABLE)
anova_results_AI

#Post-hoc test for ANOVA
pairwise_results_AI <- sable_TEE_adj_RMR %>%
  group_by(GROUP) %>%
  pairwise_t_test(adiposity_index ~ SABLE, paired = TRUE, p.adjust.method = "bonferroni") %>%
  filter(p.adj.signif <= 0.05)
pairwise_results_AI

####----- Linear mixed model, emmeans, contrasts -----####
# Build linear mixed model
model_AI <- lmer(adiposity_index ~ SABLE * GROUP + (1 | ID), data = sable_TEE_adj_RMR)
summary(model_AI)

#Confirm the number of mice in the data frame is correct...16 mice
n_distinct(sable_TEE_adj_RMR$ID) 

# Calculate estimated marginal means (emmeans) 
emm_AI <- emmeans(model_AI, ~ SABLE * GROUP, cov.reduce = mean)
emm_AI_df <- as.data.frame(emm_AI)

# Pairwise contrasts within each GROUP
contrasts_by_group_AI <- contrast(emm_AI, method = "pairwise", by = "GROUP")
# Convert to a data frame
contrasts_AI_df <- as.data.frame(contrasts_by_group_AI)
contrasts_AI_df
#Filter for significant contrasts
Sig_contrasts_AI <- contrasts_AI_df %>%
  filter(p.value <=0.05)
Sig_contrasts_AI

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_AI <- contrast(emm_AI, method = "pairwise", by = "SABLE")
# Convert to a data frame
contrasts_SABLE_AI_df <- as.data.frame(contrasts_by_SABLE_AI)
contrasts_SABLE_AI_df

#####----- Barplot of estimates from mixed model & emmeans -----####
# Define custom colors
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

barplot_emm_AI <- emm_AI_df %>%
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
    title = "NZO adiposity index (emmeans)",
    y = "Adiposity index (Fat/Lean)",
    x = "Time point",
    fill = "Treatment group") +
  scale_y_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2)) +
  #scale_y_continuous(
    #breaks = seq(0, max(emm_AI_df$emmean, na.rm = TRUE), by = 0.2)) +
  format.plot +
  theme(legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text = element_text(color="black", size=14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))
barplot_emm_AI

#####-----Scatter plot of estimates from linear mixed model & emmeans -----####
ggplot() +
  geom_point(data = emm_AI_df,
             aes(x = SABLE, y = emmean, color = GROUP),
             position = position_dodge(0.2), size = 4) +
  geom_line(data = emm_AI_df,
            aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
            position = position_dodge(0.2), linewidth = 1.5) + 
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  geom_errorbar(data = emm_AI_df,
                aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                width = 0.15, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(title = "NZO adiposity index (emmeans)",
       y = "Adiposity index (Fat/Lean)", 
       x = "Time point",
       color = "Treatment group") +
  scale_y_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2)) +
  #scale_y_continuous(
    #limits = c(0, NA),   # start at 0, upper limit adjusts automatically
    #expand = c(0, 0),    # optional: removes extra padding above/below
    #breaks = seq(0, max(emm_AI_df$emmean + emm_AI_df$SE, na.rm = TRUE), by = 0.2)) +
  #scale_y_continuous(breaks = seq(0, max(emm_AI_df$emmean, na.rm = TRUE), by = 0.2)) +
  format.plot +
  theme(legend.position = "top", 
        axis.ticks.y = element_line(),
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.text = element_text(color="black", size=14),
        axis.text.x = element_text(angle = 45, hjust = 1))

####------- Generate source df: sable_TEE_adj_RMR -------####
#functions#
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

# Resting metabolic rate (RMR): Identify 30min with lowest avg. TEE #
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
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"))

sable_RMR_data <- sable_RMR_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))

# Code including all IDs at all SABLE time points #
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

# Change RMR units from kcal_hr to kcal_day to match tee units#
lowest_windows_summary <- lowest_windows_summary %>%
  rename(RMR_kcal_hr = avg_30min_value) %>%
  mutate(RMR_kcal_day = RMR_kcal_hr*24) %>%
  group_by(ID, SABLE)


# Process sable_dwn into sable_TEE_data to get Avg daily TEE (tee) for each mouse at each time point#
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"))

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))


# Attach echoMRI_data to sable_TEE_data --> name new df as sable_TEE_adj #
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
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
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
      TRUE ~ NA_character_)) %>% 
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
    echoMRI_data %>% select(ID, SABLE, Lean, Weight, Fat, adiposity_index),
    by = c("ID", "SABLE"))

# Combine lowest_windows_summary with sable_TEE_adj #
sable_TEE_adj_RMR <- sable_TEE_adj %>%
  left_join(
    lowest_windows_summary %>% select(ID, SABLE, window_start_time, window_end_time, RMR_kcal_day),
    by = c("ID", "SABLE")) %>%
  group_by(ID, SABLE) %>%
  mutate(TEE_minus_RMR = tee - RMR_kcal_day)
