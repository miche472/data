# Experiment 2A: BW ####

#female NZO mice at peak obesity (Total of n=12 )
#one week of 10nmol/kg tirzepatide in saline (vehicle n=6, saline)
#one week of 20nmol/kg tirzepatide in saline (vehicle n=6, saline)

#Objectives ####
## 1.) Change in BW (grams) ####
  #First injection -> end of first week of injections
  #End of first dose to end of second dose

## 2.) Percent change in BW ####
  #First injection -> end of first week of injections
  #End of first dose to end of second dose
  #week 1 of regain
  #Week 2 of regain...until initial BW is restored

## 3.) How many days did it take mice to regain their original BW (on average)? ####

## 4.) Energy efficiency during each of the treatment days ####

## 5.) Rate of change in BW ####
    #Rate of change in BW pre BW loss, during BW loss week 1, BW loss week 2, and during regain

#left off: graphs, need to check accuracy of FI for 3747 on 6/19/26

#Libraries ####
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

#Format plot (LM version 2) ####
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
                             COMMENTS=="REGAIN_DAY_13" & ID %in% c(3742,3743,3750,3751,3753) ~26,
                             COMMENTS=="REGAIN_DAY_14" & ID %in% c(3742,3743,3750,3751,3753) ~27,
                             COMMENTS=="REGAIN_DAY_15" & ID %in% c(3742,3743,3750,3751,3753) ~28,
                             COMMENTS=="REGAIN_DAY_16" & ID %in% c(3742,3743,3750,3751,3753) ~29,
                             COMMENTS=="REGAIN_DAY_17" & ID %in% c(3742,3743,3750,3751,3753) ~30,
                             COMMENTS=="REGAIN_DAY_18" & ID %in% c(3742,3743,3750,3751,3753) ~31,
                             COMMENTS=="REGAIN_DAY_19" & ID %in% c(3742,3743,3750,3751,3753) ~32,
                             COMMENTS=="REGAIN_DAY_20" & ID %in% c(3742,3743,3750,3751,3753) ~33,
                             COMMENTS=="REGAIN_DAY_21" & ID %in% c(3742,3743,3750,3751,3753) ~34,
                             COMMENTS=="REGAIN_DAY_22" & ID %in% c(3742,3743,3750,3751,3753) ~35,
                             COMMENTS=="REGAIN_DAY_23" & ID %in% c(3742,3743,3750,3751,3753) ~36,
                             COMMENTS=="REGAIN_DAY_24" & ID %in% c(3742,3743,3750,3751,3753) ~37,
                             COMMENTS=="REGAIN_DAY_25" & ID %in% c(3742,3743,3750,3751,3753) ~38,
                             COMMENTS=="REGAIN_DAY_26" & ID %in% c(3742,3743,3750,3751,3753) ~39,
                             COMMENTS=="REGAIN_DAY_27" & ID %in% c(3742,3743,3750,3751,3753) ~40,
                             COMMENTS=="REGAIN_DAY_28" & ID %in% c(3742,3743,3750,3751,3753) ~41,
                             COMMENTS=="REGAIN_DAY_29" & ID %in% c(3742,3743,3750,3751,3753) ~42), 
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
#---
#---


# start delete ####
# Modify df for weeks ####
Exp2a_tracker_BW <- Exp2_tracker %>%
  group_by(ID) %>%
  mutate(
    first_inj_10nmol = DATE[COMMENTS %in% c("INJECT_DAY_2_DOSE_ONE")][1], # using day 2 vs day 1 captures the effect of the first dosage for 7 days
    days_from_first_10nmol = as.integer(DATE - first_inj_10nmol),
    BW_pct_change_wk1 = BW[days_from_first_10nmol %in% c(8)]-BW[days_from_first_10nmol %in% c(1)],
    
    first_inj_20nmol = DATE[COMMENTS %in% c("INJECT_DAY_2_DOSE_TWO")][1],
    days_from_first_20nmol = as.integer(DATE - first_inj_20nmol))
#end delete ####


#---
#---
# 1.) Change in BW (grams) ####
BW_change_gr_exp2a <- Exp2_tracker %>%
   filter(!(ID=="3747")) %>% #Filter out ID=="3747" because this mouse didn't undergo 2 weeks of weight loss
  group_by(ID) %>%
  mutate(
    BW_gr_loss_0_7 = case_when(Treatment_day == 7 ~ BW - BW[match(0, Treatment_day)], TRUE ~ NA_real_), #After 1st wk treatment
    BW_gr_loss_7_14 = case_when(Treatment_day == 14 ~ BW - BW[match(7, Treatment_day)], TRUE ~ NA_real_), #After 2nd week treatment
    BW_gr_loss_0_14 = case_when(Treatment_day == 14 ~ BW - BW[match(0, Treatment_day)], TRUE ~ NA_real_), # Total treatment period
    BW_gr_regain_wk1 = case_when(Treatment_day == 21 ~ BW - BW[match(14, Treatment_day)], TRUE ~ NA_real_),  #First week regain
    BW_gr1_regain_wk2 = case_when(Treatment_day == 28 ~ BW - BW[match(21, Treatment_day)], TRUE ~ NA_real_)) %>%
  ungroup()


#---
#---
# 2.) Percent change in BW ####
  #First injection -> end of first week of injections
  #End of first dose to end of second dose
  #week 1 of regain
  #Week 2 of regain...until initial BW is restored

##Make df for percent change in BW ####
BW_change_pct_exp2a <- Exp2_tracker %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>% #only IDs undergoing injections
  filter(Treatment_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW))) %>% #Calculate percent change in BW since BW loss started
  filter(!(Treatment_day > 14 & DRUG=="Vehicle")) %>% #these mice started BW loss in exp 2a wk 2 (were rolled into Exp 2b)
  filter(!(ID=="3747")) #Filter out ID=="3747" because this mouse didn't undergo 2 weeks of weight loss

##Line graph ####
ggplot(BW_change_pct_exp2a, aes(x=Treatment_day , y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
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
  geom_hline(yintercept=0, linetype="solid", linewidth=1, color="#62748E")+ #start BW loss
  geom_vline(xintercept=6, linetype="dotted", linewidth=1, color="#62748E")+ #start BW loss
  geom_vline(xintercept=13, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  labs(x="Day relative to first injection",
       y= "% change in BW",
       title= "Percent change in BW",
       color="Treatment", fill="Treatment")


##Bar graph ####
Plot_BW_pct_change_exp2a <- ggplot(BW_change_pct_exp2a, aes(x=Treatment_day , y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") + #454441
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  labs(x="Day relative to first injection",
       y= "% change in BW",
       title= "Percent change in BW",
       color="Treatment", fill="Treatment")
Plot_BW_pct_change_exp2a


#---
#---
## 3.) How many days did it take mice to regain their original BW (on average)? ####
#Takes less than 5 days
###Method 1 ####
exp2a_when_regain <- Exp2_tracker %>%
  filter(DRUG=="Tirzepatide") %>% 
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
         BW_pct_change = 100*((BW - first(BW)) / first(BW))) %>%
   filter(Treatment_day >14) %>%
  filter(BW_pct_change >-5 & BW_pct_change <5)

###Method 2 ####
exp2a_when_regain_try <- exp2a_when_regain %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  filter(Treatment_day >14) %>%
mutate(delta_BW= BW-first(BW)) %>%
filter(delta_BW >=0)

#---
#---
## 4.) Energy efficiency during each of the treatment days ####
# definied as change in grams of BW per kcal consumed
Exp2a_efficiency <- Exp2_tracker %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
mutate(delta_BW= BW-first(BW),
       efficiency = delta_BW/corrected_intake_kcal) %>%
  ungroup() %>%
  group_by(DATE, DRUG) %>%
  summarise(Avg_efficiency = mean(efficiency),
            SD_efficiency = sd(efficiency),
            SE_efficiency = SD_efficiency/n())

## 5.) Rate of change in BW ####
    #Rate of change in BW pre BW loss, during BW loss week 1, BW loss week 2, and during regain


