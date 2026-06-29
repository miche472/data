# APS poster --- tentative figure 2

#Graphs for poster figure with change in BW, daily FI, and energy efficiency
      #Used "Calculate BW, FI, and duration of stages (Rev. 3-27-26).R" 
      #and modified aesthetics of the graphs for APS poster

#Statistical analyses corresponding to graphs above
#Code modified from: "Change in BW during weight cycle (start 3-19-26).R"

#Started: 4-3-26
#Revised: 4-3-26

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(patchwork)

#Graph formatting ####
#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black")) # keep axis lines

# Define custom colors
custom_colors <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
custom_colors2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")

#------------------------------------#
#Create dfs with data for BW, FI, duration, and energy efficiency ####

#Calculate cumulative FI and BW change during each transition period

#Create df1 directly in this script (originally created in "FI & BW (Started 2-16-26).R")
#Prepare BW.csv
BW_to_join <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>%
  rename(COMMENTS_BW = COMMENTS)

#Prepare FI.csv
FI_to_join <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% #remove time when fed chow
  filter(corrected_intake_gr < 20 & corrected_intake_gr >= 0) %>% #removes 1-29-25 measurements 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>%
  rename(COMMENTS_FI = COMMENTS)

#Create FI_BW_joined
#Join prepared BW and FI data frames
FI_BW_joined <- FI_to_join %>%
  left_join(BW_to_join %>% 
      select(ID, DATE, BW, COMMENTS_BW),
    by = c("ID", "DATE"))

# Create df1 
#Adds variables: GROUP, DRUG, STATE, day_rel, FI_rel, FI_cum to joined BW & FI data
df1 <- FI_BW_joined %>%
  ungroup() %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(day_rel = DATE - first(DATE),
         FI_cum_INTAKE_kcal =cumsum((INTAKE_GR*3.82)),
         STATE = case_when(
           ##Baseline: First day of LFD/first day of obesity development
           #Date is the first date after start of LFD for which there was a BW and FI measurement
           ID %in% c(3706, 3707, 3708, 3709, 3710, 3711) & DATE == as.Date("2024-11-20") ~ "Baseline", 
           ID %in% c(3713, 3714, 3717, 3718, 3719) & DATE == as.Date("2024-11-27") ~ "Baseline",
           ID %in% c(3716, 3720, 3721, 3722, 3723, 3724, 3725, 3726) & DATE == as.Date("2024-12-04") ~ "Baseline",
           ID %in% c(3727, 3728, 3729) & DATE == as.Date("2024-12-11") ~ "Baseline",
           ##Peak obesity: End peak obesity period (last day of Peak obesity sable)
           #First day of calorie restriction for Weight Cycled mice. 
           #All Weight Cycled mice started restriction on the same day, so I also used this date for Control mice
           ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 
                     3725, 3728, 3729) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3716, 3717, 3718,
                     3719, 3726) & DATE == as.Date("2025-02-24") ~ "Peak obesity",
           ID %in% c(3727) & DATE == as.Date("2025-03-10") ~ "Peak obesity",
           ##BW loss: End of Sable recording for acute BW loss (i.e. start of BW maintenance period)
           #Sable day 15
           ID %in% c(3724, 3727, 3728, 3729) & DATE == as.Date("2025-04-11") ~ "BW loss", 
           ID %in% c(3707, 3709, 3711, 3723, 3725) & DATE == as.Date("2025-04-15") ~ "BW loss",
           ID %in% c(3706, 3713, 3720, 3721, 3722) & DATE == as.Date("2025-04-19") ~ "BW loss",
           ID %in% c(3708, 3710, 3714, 3716, 3718, 3726) & DATE == as.Date("2025-04-23") ~ "BW loss",
           ID %in% c(3717, 3719) & DATE == as.Date("2025-04-27") ~ "BW loss",
           ##BW maintenance: End of Sable recording for BW maintenance
           #First day of injections for all mice. (ad libitum LFD was restored on the same day as start of injections)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           ID %in% c(3713, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           ID %in% c(3714) & DATE == as.Date("2025-06-19") ~ "BW maintenance",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-22") ~ "BW maintenance",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           ##BW regain: Final day of experiment --> day of sac (End of regain)
           ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(STATE)) %>%   #Remove measurements that aren't from sign post dates (i.e. STATE)
  mutate(STATE = factor(STATE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance", "BW regain"))) 

##FI by stage: Use df1 to calculate cumulative FI for each stage ####
# Summarize cumulative FI per ID and STATE
FI_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(FI_cum_end = max(FI_cum_INTAKE_kcal, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = FI_cum_end) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `Peak obesity` - Baseline,
    kcal_peak_to_loss = `BW loss` - `Peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
FI_stage_long <- FI_stage_summary %>%
  select(ID, DRUG, GROUP,
         kcal_baseline_to_peak,
         kcal_peak_to_loss,
         kcal_loss_to_maint,
         kcal_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("kcal_"),
    names_to = "Transition",
    values_to = "kcal") %>%
  mutate(
    Transition = factor(
      Transition,
      levels = c("kcal_baseline_to_peak",
                 "kcal_peak_to_loss",
                 "kcal_loss_to_maint",
                 "kcal_maint_to_regain"),
      labels = c("BW gain",         
                 "BW loss",                 
                 "BW maintenance",          
                 "BW regain")))             

##Duration of each stage: use df1 to calculate duration (days) for each transition period ####
# Summarize duration (days) per ID and STATE
Days_stage_summary <- df1 %>%
  group_by(ID, GROUP, DRUG, STATE) %>%
  summarise(max_day_rel = max(day_rel, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,
    values_from = max_day_rel) %>%
  # Calculate days between stages
  mutate(
    Days_baseline_to_peak = `Peak obesity` - Baseline,
    Days_peak_to_loss = `BW loss` - `Peak obesity`,
    Days_loss_to_maint = `BW maintenance` - `BW loss`,
    Days_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
Days_stage_long <- Days_stage_summary %>%
  select(ID, DRUG, GROUP,
         Days_baseline_to_peak,
         Days_peak_to_loss,
         Days_loss_to_maint,
         Days_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("Days_"),
    names_to = "Transition",
    values_to = "Duration_day") %>%
  mutate(
    Transition = factor(
      Transition,
      levels = c("Days_baseline_to_peak",
                 "Days_peak_to_loss",
                 "Days_loss_to_maint",
                 "Days_maint_to_regain"),
      labels = c("BW gain",    
                 "BW loss",            
                 "BW maintenance",     
                 "BW regain")))        

# df2: Join duration (days) and cummulative FI during each phase
FI_duration_joined <- FI_stage_long %>%
  left_join(
    Days_stage_long %>% 
      select(Transition, ID, Duration_day),
    by = c("ID", "Transition"))

df2 <- FI_duration_joined %>%
  group_by(GROUP, Transition, ID) %>%
  mutate(Daily_kcal = kcal/as.numeric(Duration_day))

#Verify that there are 22 mice in each of the 4 stages 
df2 %>% 
  group_by(Transition) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO in each stage

##ΔBW in each stage: use df1 to calculate delta_BW (g) for each stage ####
# Summarize Δ BW (g) per ID and STATE
BW_stage_summary <- df1 %>% #df1 has 5 values for each mouse (one per time point)
  group_by(ID, GROUP, STATE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,         #time period (5 total)
    values_from = BW_end) %>%   #BW at end of each time period
  # Calculate change in BW between time points for each ID
  mutate(
    BWg_baseline_to_peak = `Peak obesity` - Baseline,
    BWg_peak_to_loss = `BW loss` - `Peak obesity`,
    BWg_loss_to_maint = `BW maintenance` - `BW loss`,
    BWg_maint_to_regain = `BW regain` - `BW maintenance`)

# Convert to long format
BW_stage_long <- BW_stage_summary %>%
  select(ID, GROUP,
         BWg_baseline_to_peak,
         BWg_peak_to_loss,
         BWg_loss_to_maint,
         BWg_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("BWg_"),
    names_to = "Transition",       #name of new column which lists all transition periods for all mice
    values_to = "delta_BW_g") %>%  #name of new column with delta BW values for each transition period
  mutate(
    Transition = factor(
      Transition,
      levels = c("BWg_baseline_to_peak",
                 "BWg_peak_to_loss",
                 "BWg_loss_to_maint",
                 "BWg_maint_to_regain"),
      labels = c("BW gain",    
                 "BW loss",            
                 "BW maintenance",     
                 "BW regain")))        

#Combine cumulative FI, duration, ΔBW for each phase -> Add BW_stage_long to df2 
#df2 has cumulative FI and # of days in each transition period
FI_duration_BW_joined <- df2 %>%
  left_join(BW_stage_long %>% 
      select(Transition, ID, delta_BW_g),
    by = c("ID", "Transition"))

# Calculate change in BW per kcal consumed (phase_BW_FI) for each ID during each stage
df4 <- FI_duration_BW_joined %>%
  group_by(GROUP, Transition, ID) %>%
  mutate(phase_BW_FI = delta_BW_g/kcal)

#For df4: Verify that there are 22 mice in each of the 4 transition periods 
df4 %>% 
  group_by(Transition) %>%
  summarise(n_ID = n_distinct(ID)) #we have 22 NZO in each transition period

##% change in BW during each stage ####
#Percent change relative to the "start" of each STATE rather than to Baseline BW
#to do this, modify the code that creates BW_stage_summary.
#Modification = calculate % change rather than change in grams between phases
#ΔBW during each transition: use df1 to calculate delta_BW (g) for each transition period
# Summarize Δ BW (g) per ID and STATE
BW_stage_summary_percent <- df1 %>% #df1 has 5 values for each mouse (one per time point)
  group_by(ID, GROUP, STATE) %>%
  summarise(BW_end = max(BW, na.rm = TRUE), .groups = "drop") %>% #grouped by ID and STATE, so "max BW" is equivalent to "BW" 
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATE,         #time period (5 total)
    values_from = BW_end) %>%   #BW at end of each time period
  # Calculate change in BW between time points for each ID
  group_by(ID) %>%
  mutate(
    BWp_baseline_to_peak = 100*((`Peak obesity` - Baseline)/Baseline), 
    BWp_peak_to_loss = 100*((`BW loss` - `Peak obesity`)/`Peak obesity`),
    BWp_loss_to_maint = 100*((`BW maintenance` - `BW loss`)/`BW loss`),
    BWp_maint_to_regain = 100*((`BW regain` - `BW maintenance`)/`BW maintenance`))

# Convert to long format
BW_stage_long_percent <- BW_stage_summary_percent %>%
  select(ID, GROUP,
         BWp_baseline_to_peak,
         BWp_peak_to_loss,
         BWp_loss_to_maint,
         BWp_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("BWp_"),
    names_to = "Transition",       #name of new column which lists all transition periods for all mice
    values_to = "delta_BW_perce") %>%  #name of new column with delta BW values for each transition period
  mutate(
    Transition = factor(
      Transition,
      levels = c("BWp_baseline_to_peak",
                 "BWp_peak_to_loss",
                 "BWp_loss_to_maint",
                 "BWp_maint_to_regain"),
      labels = c("BW gain",    #formerly Baseline → Peak obesity
                 "BW loss",            #formerly Peak obesity → BW loss
                 "BW maintenance",     #formerly BW loss → BW maintenance
                 "BW regain")))        #formerly BW maintenance → BW regain

# df5: Add BW_stage_long to df4. df4 has cumulative FI, # of days, and ΔBW (g) in each transition period
df5 <- df4 %>%
  left_join(
    BW_stage_long_percent %>% 
      select(Transition, ID, delta_BW_perce),
    by = c("ID", "Transition"))

#------------------------------------#
#------------------------------------#
#Fig 2A: Change in BW (grams) ####
##Stats (BW) ####

###Linear mixed model for change in BW (g) ####


##Graphs (BW) ####
###Graph BW (emmeans) ####

### Graph BW (raw) ####
plot_delta_BW_g_stage <- ggplot(df4, aes(x = Transition, y = delta_BW_g, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot+
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  labs(
    title="Change in body weight (g)",
    y = "Change in body weight (g)")
plot_delta_BW_g_stage

#Export plot to folder called "APS_figures" 
ggsave(plot_delta_BW_g_stage,
       filename="Raw_delta_BW_g_stage.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#------------------------------------#
#------------------------------------#
#Fig 2B: Daily FI ####
##Stats (FI) ####

###Linear mixed model for daily FI ####
#Build multiple linear regression model
model_FI_daily_stage <- lmer(Daily_kcal ~ Transition * GROUP + (1 | ID), data = df4)
summary(model_FI_daily_stage)

# Calculate estimated marginal means 
emm_FI_daily_stage <- emmeans(model_FI_daily_stage, ~ Transition * GROUP, cov.reduce = mean)
emm_FI_daily_stage_df <- as.data.frame(emm_FI_daily_stage)

# Pairwise contrasts within each GROUP
contrasts_by_group_FI_daily_stage <- contrast(emm_FI_daily_stage, method = "pairwise", by = "GROUP")
contrasts_by_group_FI_daily_stage_df <- as.data.frame(contrasts_by_group_FI_daily_stage)

# Pairwise contrasts within each stage of weight cycling
contrasts_by_SABLE_FI_daily_stage <- contrast(emm_FI_daily_stage, method = "pairwise", by = "Transition")
contrasts_SABLE_FI_daily_stage_df <- as.data.frame(contrasts_by_SABLE_FI_daily_stage)

##Graphs (FI) ####
###Graph FI (emmeans) ####
#Predicted daily FI by stage (emmeans for MLR) 

plot_emm_FI_daily_stage <- ggplot(emm_FI_daily_stage_df, aes(x = Transition, y = emmean, fill = GROUP)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.73) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.75, color="#454441") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", angle=45, hjust=1, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
labs(title = "Daily energy intake (kcal/day)",
  y = "Energy intake (kcal/day)")
plot_emm_FI_daily_stage

#Export plot to folder called "APS_figures" 
ggsave(plot_emm_FI_daily_stage,
       filename="emm_FI_daily_stage_plot.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")

#Conclusion: During the weight regain stage, FI was not significantly different 
#for weight cycled and control mice (emmeans, p=0.91).
#Weight cycled mice had comparable daily FI during obesity development and weight regain (emmeans, p=0.12)


### Graph FI (raw) ####
plot_daily_FI_stage <- ggplot(df4, aes(x = Transition, y = Daily_kcal, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  theme(#legend.position = "right",
    legend.position = "none",
    plot.title = element_text(size=20, hjust = 0.5, face="bold"),
    #legend.title = element_text(size = 12, face="bold"),
    #legend.text = element_text(size = 12),
    legend.title=element_blank(),
    legend.text=element_blank(),
    axis.text.x = element_text(size= 20, color="black", hjust=1, vjust=1, angle=45),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size = 20),
    panel.border = element_blank()) +
  format.plot +
  labs(title="Daily energy intake (kcal/day)",
       y = "Energy intake (kcal/day)") 
plot_daily_FI_stage

#Export plot to folder called "APS_figures" 
ggsave(plot_daily_FI_stage,
       filename="Raw_daily_FI_stage_plot.png", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/APS_figures")
  
#------------------------------------#
#------------------------------------#
#Fig 2C: Energy efficiency ####
## Stats (Energy efficiency) ####
### Linear mixed model for energy efficiency ####


##Graphs (Energy efficiency) ####
###Graph Energy efficiency (emmeans) ####


### Graph Energy efficiency (raw) ####
#(use!) Graph Energy Efficiency for ONLY Obesity development & Weight regain ####
df4_efficient <- df4 %>%
  filter(Transition %in% c("Obesity development", "Weight regain"))

ggplot(df4_efficient, aes(x = Transition, y = phase_BW_FI, fill = GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), 
           width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8),
                width = 0.3, 
                color="black") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors2) +
  #theme_minimal() + #removes the lines that make the graph a box
  format.plot+
  theme(legend.position = "right",
        plot.title = element_text(size=16, hjust = 0.5, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size= 13, color="black"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  labs(title="Figure 2C: Energy efficiency",
       #x = "Stage of weight cycle",
       y = "ΔBW (g) per kcal consumed",
       fill = "Treatment group",
       color = "Treatment group")




