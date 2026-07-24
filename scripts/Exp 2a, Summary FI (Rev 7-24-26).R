# Summary Food intake for Exp 2A

#Started: 7-23-26
#Revised: 7-24-26

#Graph and analyze average cumulative intake by week and average daily intake by week
  #One week before start of BW loss
  #During 10nmol/kg
  #During 20nmol/kg
  #Week 1 of regain
  #week 2 of regain

#Left off: daily food intake by week (7/24/26)

library(mmand)
library(pacman)
library(this.path)
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
library(lme4)
library(grid)

#Format plot (LM version)
format.plot_LM2 <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"), 
  panel.border = element_blank(),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black"),
  plot.title = element_text(size=17, hjust = 0.5, face="bold", vjust=2),
  legend.title=element_text(size=15, face="bold"),
  legend.text=element_text(size=13),
  axis.title.x = element_text(face="bold", size= 15),
  axis.text.x = element_text(size= 13, angle=25, vjust=0.5, hjust=0.7),
  axis.title.y = element_text(face="bold", size= 15),
  axis.text.y = element_text(size = 13))
# Define custom colors
custom_colors1_GLP_Exp2 <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")
custom_colors2_GLP_Exp2 <- c("Tirzepatide" = "#1E90FF", "Vehicle" = "#8B8989")
custom_colors3_GLP_Exp2 <- c("Tirzepatide" = "#104E8B", "Vehicle" = "#403d3c")


# Update BW.csv & FI_LM.csv ####
#Create df ####
# bodyweight and food intake

# libs 
pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand
)

cohort_csv_files <- tibble(
  filepath = list.files("../data", full.names = TRUE)) %>% 
  filter(
    grepl("COHORT_[0-9]+[0-9]*.csv", filepath)) #now we can used cohort > 10
cohort_csv_files

cohort_open_files <- cohort_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map_dfr(
    ., function(X){
      read_csv(X$filepath) %>% 
        select(ID, FOOD_WEIGHT_START_G, FOOD_WEIGHT_END_G, DATE, DIET, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS) %>% 
        mutate(
          INTAKE_GR = (FOOD_WEIGHT_START_G - FOOD_WEIGHT_END_G),
          DATE = lubridate::mdy(DATE)
        ) %>% 
        select(ID, INTAKE_GR, DATE, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS) %>% 
        rename(
          BW = BODY_WEIGHT_G
        ) %>% 
        mutate(BW=as.numeric(BW), ID=as.factor(ID))})

# load food description
food_desc <- read_csv("../data/food_description.csv")

# load metadata
#I changed the META.csv file on my local computer to include COHORT 19...
#NOTE: META.csv --> cohort 19 will be wiped away the next time I pull from origin ####
metadata <- read_csv("../data/META.csv") %>% 
  select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# output food-intake file
FI_LM <- cohort_open_files %>%
  select(ID, DIET_FORMULA, INTAKE_GR, DATE, COMMENTS) %>%
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    delta_alt = {
      intake_idx <- !is.na(INTAKE_GR) #creates a logical vector (TRUE/FALSE) where rows that have INTAKE_GR=NA --> FALSE and rows with a value for INTAKE_GR -->TRUE
      intake_dates <- DATE[intake_idx] #Unconfirmed: only keeps rows for which intake_idx is TRUE
      
      # compute differences only on valid intake rows
      diffs <- c(NA, as.numeric(diff(intake_dates)))
      
      # create full-length vector and fill only intake rows
      out <- rep(NA_real_, n())
      out[intake_idx] <- diffs
      out
    }
  ) %>%
  mutate(delta_measurement = DATE - lag(DATE)) %>% #just use to remove first observation for each mouse
  drop_na(delta_measurement) %>% #just use to remove first observation for each mouse
  mutate(corrected_intake_gr = INTAKE_GR / as.numeric(delta_alt)) %>%
  left_join(., food_desc, by = "DIET_FORMULA") %>%
  mutate(corrected_intake_kcal = corrected_intake_gr * KCAL_G) %>%
  left_join(., metadata, by = "ID")  %>%
  select(-delta_measurement)


# output bodyweight file
BW <- cohort_open_files %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, BW, DATE,COMMENTS) %>% 
  drop_na(BW) %>% 
  left_join(., metadata, by = "ID")

write_csv(x = FI_LM, "../data/FI_LM.csv")
write_csv(x = BW, "../data/BW.csv")

#Read in BW and FI ####
#Read in BW and filter for cohort 19 (Spring 2026 NZO mice)
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19)

#Read in FI and filter for cohort 19 (Spring 2026 NZO mice)
FI_LM_COHORT19 <- read_csv("~/Documents/GitHub/data/data/FI_LM.csv") %>%
  filter(COHORT == 19)

#Create df with BW and FI
BW_FI_19 <- BW_COHORT19 %>% #Join FI and BW
  left_join(
    FI_LM_COHORT19 %>% 
      select(ID, INTAKE_GR, DATE, delta_alt, corrected_intake_gr, corrected_intake_kcal, KCAL_G),
    by = c("ID", "DATE")) %>%
  mutate(ID = as.factor(ID)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(day_rel = DATE - first(DATE),
         day_rel = as.numeric(day_rel))
replace_na(list(#INTAKE_GR=0, 
  #delta_alt=0, 
  #corrected_intake_gr=0,
  #corrected_intake_kcal=0, 
  KCAL_G=3.82))

# Make df for analysis
Exp2_tracker <- BW_FI_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0)) %>%
  mutate( 
    STAGE = case_when(ID %in% c(3742, 3743, 3744, 3745, 3748, 3749, 3750, 3751) & DATE< "2026-06-03" ~ "BW gain",
                      ID %in% c(3742, 3743, 3749, 3750, 3751) & DATE> "2026-06-02" & DATE<"2026-06-17" ~ "BW loss",
                      ID %in% c(3744, 3745, 3748) & DATE> "2026-06-02" & DATE<"2026-06-20" ~ "BW loss",
                      ID %in% c(3742, 3743, 3750, 3751) & DATE>"2026-06-16" ~ "BW regain",
                      
                      ID %in% c(3746, 3747, 3752, 3753) & DATE< "2026-06-08" ~ "BW gain",
                      
                      ID == 3747 & DATE>= "2026-06-08" & DATE<"2026-06-18" ~ "BW loss",
                      ID == 3747 & DATE> "2026-06-17" ~ "BW regain",
                      
                      ID == 3753 & DATE>= "2026-06-08" & DATE<"2026-06-22" ~ "BW loss",
                      ID == 3753 & DATE> "2026-06-21" ~ "BW regain",
                      
                      
                      ID %in% c(3746, 3752) 
                      & DATE> "2026-06-07"& DATE<"2026-06-23" ~ "BW loss"),
    
    
    DRUG = case_when(ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
                     ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide"),
    # All days since treatment started (including days w/o injection)
    Treatment_day =case_when(COMMENTS=="INJECT_DAY_1_DOSE_ONE"~0,
                             COMMENTS=="SABLE_DAY_4_INJECT_DAY_1_DOSE_ONE"~0, #For 3747,3752,3753
                             
                             COMMENTS=="INJECT_DAY_2_DOSE_ONE"~1,
                             COMMENTS=="INJECT_DAY_3_DOSE_ONE"~2,
                             COMMENTS=="INJECT_DAY_4_DOSE_ONE"~3,
                             COMMENTS=="INJECT_DAY_5_DOSE_ONE"~4,
                             COMMENTS=="INJECT_DAY_6_DOSE_ONE"~5,
                             COMMENTS=="INJECT_DAY_7_DOSE_ONE"~6,
                             
                             COMMENTS=="INJECT_DAY_1_DOSE_TWO"~7,
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_1_DOSE_TWO"~7, #For 3746,3747
                             
                             COMMENTS=="INJECT_DAY_2_DOSE_TWO"~8,
                             COMMENTS=="SABLE_DAY_2_INJECT_DAY_2_DOSE_TWO"~8, #For 3747
                             
                             COMMENTS=="INJECT_DAY_3_DOSE_TWO"~9,
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_3_DOSE_TWO"~9, #For 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_3_INJECT_DAY_3_DOSE_TWO"~9, #For 3747
                             
                             COMMENTS=="INJECT_DAY_4_DOSE_TWO"~10,
                             COMMENTS=="SABLE_DAY_2_INJECT_DAY_4_DOSE_TWO"~10, #For 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_4_REGAIN_DAY_1"~10, #For 3747
                             
                             COMMENTS=="INJECT_DAY_5_DOSE_TWO"~11,
                             COMMENTS=="SABLE_DAY_3_INJECT_DAY_5_DOSE_TWO"~11, #For 3742,3743,3750,3751,3753
                             
                             COMMENTS=="INJECT_DAY_6_DOSE_TWO"~12,
                             COMMENTS=="SABLE_DAY_4_INJECT_DAY_6_DOSE_TWO"~12, #for 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_6_DOSE_TWO"~12, #For 3746,3752
                             
                             COMMENTS=="INJECT_DAY_7_DOSE_TWO"~13,
                             COMMENTS=="SABLE_DAY_5_INJECT_DAY_7_DOSE_TWO" ~13, #for 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_2_INJECT_DAY_7_DOSE_TWO"~13, #For 3746,3752
                             
                             COMMENTS=="SABLE_DAY_6_REGAIN_DAY_1"~14,
                             COMMENTS=="SABLE_DAY_1" & ID %in% c(3744, 3745, 3748) & DATE>"2026-06-02" & DATE < "2026-06-24"~14, #For 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_3" & ID %in% c(3746,3752) & DATE>"2026-06-02" & DATE < "2026-06-24"~14,
                             ID==3749 & DATE == "2026-06-17"~14,
                             
                             COMMENTS=="SABLE_DAY_2" & ID %in% c(3744, 3745, 3748) &  DATE>"2026-06-02" & DATE <"2026-06-24"~15,
                             ID == 3749 & DATE == "2026-06-18"~15,
                             COMMENTS=="SABLE_DAY_4" & ID %in% c(3746,3752) & DATE>"2026-06-02" & DATE <"2026-06-24"~15,
                             
                             COMMENTS=="SABLE_DAY_3" & ID %in% c(3744, 3745, 3748) &  DATE>"2026-06-02" & DATE <"2026-06-24"~16,
                             ID == 3749 & DATE == "2026-06-19"~16,
                             COMMENTS=="REGAIN_DAY_3" & ID %in% c(3742,3743,3750,3751,3753)~16, 
                             
                             COMMENTS=="SABLE_DAY_4" & ID %in% c(3744, 3745, 3748) & DATE>"2026-06-02" & DATE<"2026-06-24"~17,
                             ID==3749 & DATE == "2026-06-20"~17,
                             COMMENTS=="REGAIN_DAY_4" & ID %in% c(3742,3743,3750,3751,3753)~17,
                             
                             #ID==3749 & DATE == "2026-06-21"~18,
                             COMMENTS=="REGAIN_DAY_5" & ID %in% c(3742,3743,3750,3751,3753)~18,
                             
                             #ID == 3749 & DATE == "2026-06-22"~19,
                             COMMENTS=="REGAIN_DAY_6" & ID %in% c(3742,3743,3750,3751,3753) ~19,
                             COMMENTS=="REGAIN_DAY_7" & ID %in% c(3742,3743,3750,3751,3753) ~20,
                             COMMENTS=="REGAIN_DAY_8" & ID %in% c(3742,3743,3750,3751,3753) ~21,
                             COMMENTS=="REGAIN_DAY_9" & ID %in% c(3742,3743,3750,3751,3753) ~22,
                             COMMENTS=="REGAIN_DAY_10" & ID %in% c(3742,3743,3750,3751,3753) ~23,
                             COMMENTS=="REGAIN_DAY_11" & ID %in% c(3742,3743,3750,3751,3753) ~24,
                             COMMENTS=="REGAIN_DAY_12" & ID %in% c(3742,3743,3750,3751,3753) ~25,
                             COMMENTS=="REGAIN_DAY_13" & ID %in% c(3742,3743,3750,3751,3753) ~26), 
    Dose_step = case_when(
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-06-02" & DATE < "2026-6-10" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-6-09" & DATE < "2026-06-17" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      
      ID %in% c(3747, 3753) & STAGE == "BW loss" & 
        DATE > "2026-06-05" & DATE < "2026-6-14" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID %in% c(3747, 3753) & STAGE == "BW loss" & 
        DATE > "2026-6-13" & DATE < "2026-06-21" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      
      ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) 
      & STAGE == "BW loss" & DRUG == "Vehicle" ~ "Vehicle")) %>%
  mutate(DRUG = as.factor(DRUG)) %>%
  drop_na(STAGE) #This is risky, but needed to remove the period when the exp 2a vehicle mice were repurposed for BW loss in exp 2b



#Before injections: Cumulative intake for a 7 day "window" before COMMENTS=="INJECT_DAY_1_DOSE_ONE"
#During initial dose: Cumulative intake for a 7 days window after COMMENTS=="INJECT_DAY_1_DOSE_ONE" 
      #(not including this day)...captures days when Dose_step=="10nmol/kg"
#During escalated dose: Cumulative intake for a 7 days window after COMMENTS=="INJECT_DAY_1_DOSE_TWO" 
      #(not including this day)...captures days when Dose_step=="20nmol/kg"
#Cumulative intake for 7 day window after injections ceased
#Cumulative intake for 7 day window after "regain day 7"

Exp2_tracker_week <- Exp2_tracker %>%
  group_by(ID) %>%
  mutate(
    #first_inj = DATE[COMMENTS %in% c("INJECT_DAY_1_DOSE_ONE", "SABLE_DAY_4_INJECT_DAY_1_DOSE_ONE")][1],
    first_inj = DATE[COMMENTS %in% c("INJECT_DAY_2_DOSE_ONE")][1], # captures the effect of the first dosage for 7 days
    days_from_first = as.integer(DATE - first_inj),
    Relative_week = floor(days_from_first / 7)) %>%  # Relative week
  ungroup()

weekly_intake <- Exp2_tracker_week %>%
  group_by(ID, DRUG, Relative_week) %>%
  summarise(
    weekly_food_g = sum(INTAKE_GR, na.rm = TRUE),
    weekly_food_kcal = sum(INTAKE_GR * 3.82, na.rm = TRUE),
    n_days = n(),
    .groups = "drop")

daily_intake_by_wk <- Exp2_tracker_week %>%
  group_by(ID, DRUG, Relative_week) %>%
  summarise (daily_food_g_mean = mean(corrected_intake_gr),
             daily_food_g_sum_7 = (sum(INTAKE_GR, na.rm = TRUE)/7),
             diff_method = daily_food_g_mean-daily_food_g_sum_7)

#Filter out where necessary
try_filter <- Exp2_tracker_week %>%
  filter(!(Relative_week == 2 & DRUG=="Vehicle")) %>% #these mice started BW loss in exp 2a wk 2 (were rolled into Exp 2b)
  filter(!(Relative_week %in% c(-10, -9))) %>% #Data isn't clean this early on
  filter(!(ID=="3747")) #Filter out ID=="3747" because this mouse didn't undergo 2 weeks of weight loss

#Graph FI during each week ####

#Prepare df
Weekly_intake_plot <- try_filter %>%
  group_by(ID, DRUG, Relative_week) %>%
  summarise(
    weekly_food_g = sum(INTAKE_GR, na.rm = TRUE),
    weekly_food_kcal = sum(INTAKE_GR * 3.82, na.rm = TRUE),
    n_days = n(),
    .groups = "drop") %>%
  mutate(Relative_week = as.numeric(Relative_week)) %>%
  filter(Relative_week > -2) 

#Line graph
ggplot(Weekly_intake_plot, aes(x=Relative_week , y=weekly_food_g, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values = custom_colors1_GLP_Exp2) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  geom_vline(xintercept=0, linetype="dotted", linewidth=1, color="#62748E")+ #start BW loss
  geom_vline(xintercept=1, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  geom_vline(xintercept=2, linetype="dotted", linewidth=1, color="#62748E")+ #end BW loss
  labs(x="Relative week",
       y= "Cumulative intake (g)",
       title= "Cumulative food intake (g/wk)",
       color="Treatment", fill="Treatment")


#Bar graph- cumulative FI (grams per week)
plot_weekly_FI_exp2a <- ggplot(Weekly_intake_plot, aes(x=Relative_week , y=weekly_food_g, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") + #454441
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  labs(x="Week",
       y= "Cumulative intake (g)",
       title= "Cumulative food intake (grams/week)",
       color="Treatment", fill="Treatment")
plot_weekly_FI_exp2a

##LMM: cumulative FI ####
For_stats_weekly_intake <- Weekly_intake_plot %>%
  mutate(Relative_week = as.factor(Relative_week)) 
  

#Build multiple linear regression model for RMR not adjusted for BW or lean #
model_cumulative_FI_gr <- lmer(weekly_food_g ~ Relative_week * DRUG + (1 | ID), data = For_stats_weekly_intake)
summary(model_cumulative_FI_gr)

#Calculate estimated marginal means #
emm_cumulative_FI_gr <- emmeans(model_cumulative_FI_gr, ~ Relative_week * DRUG, cov.reduce = mean)
emm_cumulative_FI_gr_df <- as.data.frame(emm_cumulative_FI_gr)

# Pairwise contrasts within each GROUP
contrasts_by_Treatment_cumulative_FI_gr <- contrast(emm_cumulative_FI_gr, method = "pairwise", by = "DRUG")
contrasts_by_Treatment_cumulative_FI_gr_df <- as.data.frame(contrasts_by_Treatment_cumulative_FI_gr)

# Pairwise contrasts within each Relative_week (time point)
contrasts_by_week_cumulative_FI_gr <- contrast(emm_cumulative_FI_gr, method = "pairwise", by = "Relative_week")
contrasts_by_week_cumulative_FI_gr <- as.data.frame(contrasts_by_week_cumulative_FI_gr)

# Daily food intake by week ####
## Bar graph ####
## LMM - daily FI by week ####


