# Incretin receptor agonist optimization experiment 2

#started: 5-26-26
#revised: 

# ----Goals/Objectives ----####
#Objective 1: Track BW and FI during Sable acclimation & BW loss
#Objective 2: Assign mice to treatment groups (TZP or vehicle)
#Objective 3: Monitor body composition for n=12 mice
#Objective 4: Assess energy expenditure during peak obesity reading
    #to verify proper functioning of Sable system

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
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb","Survodutide" = "#c93618", "Vehicle" = "#403d3c")
custom_colors_GLPB <- c("Tirzepatide" = "lightblue","Survodutide" = "lightgreen", "Vehicle" = "darkgray")


# libraries
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

#------------------------------------------------------------#
#Objective 1

# Update BW.csv & FI_LM.csv ####
#Create df ####
# bodyweight and food intake
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
#I changed the META.csv file on my local computer to include COHORT 19...this 
#will be wiped away the next time I pull from origin
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
      out}) %>%
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

##Read in BW and FI
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
replace_na(list(KCAL_G=3.82))

#-----------------------------------------------------------------#
#-------------------------Goal 1: sable acclimation-------------------------#####
#-----------------------------------------------------------------#
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
    # NOTE: dates for the two rounds of sable recording will be different
    STAGE = case_when(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                      & DATE< "2026-06-02" ~ "BW gain",
                      ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                      & DATE> "2026-06-01" ~ "BW loss",
                      
                      ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                      & DATE< "2026-06-02" ~ "BW gain",
                      ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                      & DATE> "2026-06-01" ~ "BW loss"),
    
    DRUG = case_when(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) ~ "Vehicle",
                     ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) ~ "Tirzepatide"),
    
    Treatment_day =case_when(COMMENTS=="INJECT_DAY_ONE"~0, # All days since treatment started (including days w/o injection)
                             COMMENTS=="INJECT_DAY_TWO"~1),
    
    Dose_step = case_when(
                          ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)
                          & STAGE == "BW loss" & DATE< "2026-5-03"& DRUG=="Tirzepatide" ~ "10nmol/kg",
                          ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                          & STAGE == "BW loss" & DATE> "2026-5-02" & DATE< "2026-05-17" 
                          & DRUG == "Tirzepatide" ~ "20nmol/kg",
                          ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)
                          & STAGE == "BW loss" & DATE >"2026-5-16" & DATE <"2026-5-22"
                          & DRUG == "Tirzepatide" ~ "30nmol/kg",
                          
                          ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753) 
                          & STAGE == "BW loss" & DRUG == "Vehicle" ~ "Vehicle")) %>%
  mutate(DRUG = as.factor(DRUG))
 

#---
##Plot BW & FI  ####
# BW throughout study
ggplot(Exp2_tracker, aes(x = DATE, y = BW, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="Day",
       y= "Body weight (grams)",
       title= "Body weight")

#---
#FI throughout study
Exp2_tracker_plot_FI <- Exp2_tracker %>%
  filter(corrected_intake_kcal>0)

ggplot(Exp2_tracker_plot_FI, aes(x = DATE, y = corrected_intake_kcal, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="Day",
       y= "Daily food intake (kcal)",
       title= "Daily food intake (kcal)")

#------------------------------------------#

##Stage: Acclimation (sable peak obesity) ####
Exp2_acc_plot <- Exp2_tracker %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  filter(DATE> "2026-05-20") %>%
  mutate(
        day_rel = DATE-first(DATE),
        BW_pct_change = 100*((BW - first(BW)) / first(BW)),
        BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR),
        BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal))

### BW ####
#acclimation
ggplot(Exp2_acc_plot, aes(x = day_rel, y = BW, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="Day",
       y= "Body weight (grams)",
       title= "BW: sable acclimation")

#### BW % (acclimation)
ggplot(Exp2_acc_plot, aes(x = day_rel, y = BW_pct_change, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="Day",
       y= "ΔBody weight (%) ",
       title= "Percent ΔBW during sable acclimation")

###FI ####
#Acclimation
Exp2_acc_plot_FI <- Exp2_acc_plot %>%
  filter(corrected_intake_kcal>0)

ggplot(Exp2_acc_plot_FI, aes(x = DATE, y = corrected_intake_kcal, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="Day",
       y= "Daily food intake (kcal)",
       title= "Daily food intake (kcal)")


#------------------------------------------#
##Stage: BW loss ####
Exp2_BWloss_plot <- Exp2_tracker %>%
  filter(Treatment_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    BW_pct_change = 100*((BW - first(BW)) / first(BW)),
    BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR),
    BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal))

### BW  ####
#during BW loss
#### Percent ΔBW during BW loss
ggplot(Exp2_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  #stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1), color=DRUG), vjust = -1, size = 4, 
  #fill = "white", fontface="bold", label.padding = unit(0.15, "lines"), show.legend=FALSE) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  labs(x="Day",
       y= "Δ Body weight (%) ",
       title= "Percent ΔBW during BW loss",
       color="Treatment", fill="Treatment")


### FI ####
#Daily FI (kcal) during BW loss


#-----------------------------------------------------------------#
#-------------------------Goal 2: Assign groups-------------------------#####
#-----------------------------------------------------------------#
#Either 6 vehicle/6 TZP or 4 vehicle/8 TZP

###Method: Use adiposity index to assign mice to groups ####
#Need to save the echoMRI data files with the correct formatting and title
#to the echoMRI folder in the GitHub data folder. Then I need to re-run
#data_proc so that an updated echoMRI.csv is written.
#The code below will pull data from echoMRI.csv

#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_data_2 <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  filter(Date == "2026-05-21") %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index) %>%
  ungroup()

#---
## Two groups (n=6) ####
#Create pattern for group assignment
pattern_2B <- c("Group1", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group1") 

echoMRI_data_2B <- echoMRI_data_2 %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern_2B, length.out = n())) %>%
  ungroup() %>%
  mutate(treatment = as.factor(treatment))

#T-test to see if the mean adiposity index is significantly different for Restrict and Ad libitum groups
pacman::p_load(
  tidyverse,
  ggplot2,
  ggpubr,
  furrr)
library(tidyverse)
library(ggpubr)
library(rstatix)

#Get summary stats
echoMRI_data_2B %>%
  group_by(treatment) %>%
  get_summary_stats(adiposity_index, type = "mean_sd") 

##Test pre-requisites for T Test and perform T test
#Test for outliers  
echoMRI_data_2B %>%
  group_by(treatment) %>%
  identify_outliers(adiposity_index)

#Test for normality
echoMRI_data_2B %>%
  group_by(treatment) %>%
  shapiro_test(adiposity_index)

#Equal variance (p value is 0.329 suggesting that there is not a significant difference in the 
                #variance within each group)
echoMRI_data_2B %>% levene_test(adiposity_index ~ treatment)

#Student t test (assumes equal variance)
stat.test_2B <- echoMRI_data_2B %>%
  t_test(adiposity_index ~ treatment, var.equal = TRUE) %>%
  add_significance()
stat.test_2B
#p= 0.871 --> mean adiposity index is not sig diff between 2 treatment groups

#Graph adiposity index for treatment groups 1 and 2 
ggplot(echoMRI_data_2B, aes(x = treatment, y = adiposity_index, fill = treatment)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "Adiposity Index (fat mass/lean mass)")

#For Weight:
#Get summary stats
echoMRI_data_2B %>%
  group_by(treatment) %>%
  get_summary_stats(Weight, type = "mean_sd") 

##Test pre-requisites for T Test and perform T test
#Test for outliers  
echoMRI_data_2B %>%
  group_by(treatment) %>%
  identify_outliers(Weight)

#Test for normality
echoMRI_data_2B %>%
  group_by(treatment) %>%
  shapiro_test(Weight)

#qqplot by group
ggqqplot(echoMRI_data_2B, x = "adiposity_index", facet.by = "treatment")

#Equal variance (p value is 0.135 suggesting that there is not a significant difference in the 
#variance within each group)
echoMRI_data_2B %>% levene_test(Weight ~ treatment)

#Student t test (assumes equal variance)
stat.test_2B_BW <- echoMRI_data_2B %>%
  t_test(Weight ~ treatment, var.equal = TRUE) %>%
  add_significance()
stat.test_2B_BW
#p= 0.847 --> mean adiposity index is not sig diff between 2 treatment groups

#Graph Weight for treatment groups 1 and 2 
ggplot(echoMRI_data_2B, aes(x = treatment, y = Weight, fill = treatment)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "Body weight (g)")


#---
##Three groups (n=4, n=4, n=4 TZP) ####
#Create pattern for group assignment
pattern_2C <- c("Group1", 
                "Group2", 
                "Group3", 
                "Group1", 
                "Group2", 
                "Group3", 
                "Group3", 
                "Group2", 
                "Group1", 
                "Group3", 
                "Group2", 
                "Group1") 

echoMRI_data_2C <- echoMRI_data_2 %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern_2C, length.out = n())) %>%
  ungroup() %>%
  mutate(treatment = as.factor(treatment))

#Test for outliers  
echoMRI_data_2C %>%
  group_by(treatment) %>%
  identify_outliers(adiposity_index)

#Test for normality
echoMRI_data_2C %>%
  group_by(treatment) %>%
  shapiro_test(adiposity_index)

#qqplot by group
ggqqplot(echoMRI_data_2C, x = "adiposity_index", facet.by = "treatment")

#Test for equal variance
echoMRI_data_2C %>% levene_test(adiposity_index ~ treatment)

# Run the ANOVA
anova_model_2C <- aov(adiposity_index ~ treatment, data = echoMRI_data_2C)

# View the ANOVA table
summary(anova_model_2C)

#Graph adiposity index for treatment groups 1, 2, 3 
ggplot(echoMRI_data_2C, aes(x = treatment, y = adiposity_index, fill = treatment)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "Adiposity Index (fat mass/lean mass)")

#---
## LOOK HERE! ####
#IF I combine groups 1 and 3 would that create two group Group1/3 (n=8) and Group2 (n=4)
#with equal variance, normal distribution, no outliers, and means that aren't significantly different?

echoMRI_data_2C_combined <- echoMRI_data_2C %>%
  mutate(treatment2 = if_else(treatment == "Group2", "GroupA", "GroupB")) %>%
  ungroup() %>%
  mutate(treatment2 = as.factor(treatment2))

#Test for outliers  
echoMRI_data_2C_combined %>%
  group_by(treatment2) %>%
  identify_outliers(adiposity_index)

#Test for normality
echoMRI_data_2C_combined %>%
  group_by(treatment2) %>%
  shapiro_test(adiposity_index)

#qqplot by group
ggqqplot(echoMRI_data_2C_combined, x = "adiposity_index", facet.by = "treatment2")

#Test for equal variance
echoMRI_data_2C_combined %>% levene_test(adiposity_index ~ treatment2)

#Student t test (assumes equal variance)
stat.test_2C_combined <- echoMRI_data_2C_combined %>%
  t_test(adiposity_index ~ treatment2, var.equal = TRUE) %>%
  add_significance()
stat.test_2C_combined

#Graph adiposity index for treatment groups A (1 and 3) and group B (2)
ggplot(echoMRI_data_2C_combined, aes(x = treatment2, y = adiposity_index, fill = treatment2)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "Adiposity Index (fat mass/lean mass)")

#Graph BW for treatment groups A (1 and 3) and group B (2)
ggplot(echoMRI_data_2C_combined, aes(x = treatment2, y = Weight, fill = treatment2)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "BW (g)")

#---
## Two groups (n=8 & n=4) ####
pattern_2D <- c("Group1", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group1", 
                "Group1", 
                "Group1", 
                "Group1", 
                "Group2", 
                "Group1", 
                "Group2", 
                "Group1") 

echoMRI_data_2D <- echoMRI_data_2 %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern_2D, length.out = n())) %>%
  ungroup() %>%
  mutate(treatment = as.factor(treatment))

#Test for outliers  
echoMRI_data_2D %>%
  group_by(treatment) %>%
  identify_outliers(adiposity_index)

#Test for normality
echoMRI_data_2D %>%
  group_by(treatment) %>%
  shapiro_test(adiposity_index)

#qqplot by group
ggqqplot(echoMRI_data_2D, x = "adiposity_index", facet.by = "treatment")

#Test for equal variance
echoMRI_data_2D %>% levene_test(adiposity_index ~ treatment)

ggplot(echoMRI_data_2D, aes(x = treatment, y = adiposity_index, fill = treatment)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.position = "right") + 
  format.plot +
  labs(x=NULL) +labs(y= "Adiposity Index (fat mass/lean mass)")

#-----------------------------------------------------------------#
#-------------------------Goal 3: Body composition-------------------------#####
#-----------------------------------------------------------------#

#Echo dates (don't have data for all 12 mice on all dates)
# 2026-04-29 -> 3747, 3748, 3749, 3750, 3751, 3752, 3753
# 2026-05-07 -> 3742, 3743, 3744, 3745, 3746
# 2026-05-22 -> 3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753

# These mice were called "Naive" in the context of incretin agonist experiment 1
#Need to save the echoMRI data files with the correct formatting and title
#to the echoMRI folder in the GitHub data folder. Then I need to re-run
#data_proc so that an updated echoMRI.csv is written.

#Excerpt from data_proc.R to create the echo data frame and write the csv
# load metadata
metadata <- read_csv("../data/META.csv") %>% 
  select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# echo MRI
echomri_csv_files <- tibble(
  filepath = list.files("../data/echoMRI", full.names = TRUE)) %>% 
  filter(grepl("*.xlsx", filepath)) 
echomri_csv_files

echomri_open_files <- echomri_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map(., function(X){
    readxl::read_xlsx(X$filepath) %>% 
      select(Label, Fat, Lean, Weight, TimeDateDura) %>% 
      rename(ID = Label) %>% 
      separate_wider_delim(TimeDateDura, delim = ";", names = c("Date", "A", "B")) %>% 
      select(-A, -B) %>% 
      separate_wider_delim(Date, delim = " ", names = c("hms", "month", "day", "year")) %>% 
      mutate(day = gsub(",", "", day),
             Date = paste(year, month, day, sep = "-"),
             Date = lubridate::ymd(Date),
             ID =  as.factor(ID)) %>% 
      select(-hms, -month, -day, -year)}) %>% 
  bind_rows() %>% 
  left_join(., metadata, by = "ID")
echomri_open_files

# compare adiposity index = fat / lean
echomri_data <- echomri_open_files %>% 
  mutate(adiposity_index = Fat / Lean) %>% 
  group_by(ID) %>% 
  mutate(n_measurement = as.numeric(as.factor(Date)))
echomri_data

write_csv(x = echomri_data, "../data/echomri.csv")

##The code below will pull data from echoMRI.csv
#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_data_19 <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) 

echoMRI_data_exp2 <- echoMRI_data_19 %>%
  filter(ID %in% c(3742,3743,3744,3745,3746,3747,3748,3749,3750, 3751,3752,3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(
         #DRUG = case_when(ID %in% c() ~ "Vehicle",
                          #ID %in% c() ~ "Tirzepatide"),
         #DRUG = as.factor(DRUG),
         stage = case_when(n_measurement %in% c(1,2,3) ~ "BW gain"))


#Just BW gain
echoMRI_calc_exp2 <-echoMRI_data_exp2 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  filter(stage=="BW gain") %>%
  mutate(delta_AI= adiposity_index - first(adiposity_index),
         delta_Weight = Weight - first(Weight),
         delta_Lean = Lean - first(Lean),
         delta_Fat = Fat-first(Fat),
         Lean_per_Weight = delta_Lean/delta_Weight,
         Fat_per_Weight = delta_Fat/delta_Weight,
         AI_pct = (delta_AI/first(adiposity_index))*100,
         Lean_pct = ((delta_Lean/first(Lean))*100), #lean mass gained for every g of BW gained
         Fat_pct = (delta_Fat/first(Fat))*100,
         Weight_pct = (delta_Weight/first(Weight))*100)
         






#-----------------------------------------------------------------#
#-------------------------Goal 3: Body composition-------------------------#####
#-----------------------------------------------------------------#


















