#

#Version 1: 5-29-2025 LM
# Context: NZO mice were previously assigned to two diet_group (ad lib or restricted) based upon adiposity index
#Goal: Assign n=22 female NZO mice (cohorts 3-5) to rtioxa47 or vehicle injection groups during BW regain period. 
  #Use adiposity index as criteria for group assignment.
  #Generate 4 groups with intergroup homogeneity between ad lib/rtioxa47 & ad lib/vehicle 
  #and between restricted/rtioxa47 & restricted/vehicle

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()

#----ADIPOSITY INDEX----

#open echo_mri_data
echomri_data <- read_csv("../data/echomri.csv")

echomri_data_2<- echomri_data %>%
  filter(COHORT > 2 & COHORT < 6) %>%
  filter(Date == "2025-05-27") %>%
  select(ID, adiposity_index) %>%
  mutate(diet_group = if_else(ID %in% c(
    3708, 3714, 3720, 3721, 3710,
    3722, 3723, 3724, 3725, 3727, 3728, 3729
  ), "RESTRICTED", "AD_LIB"))

#Create a new column called "treatment" to indicate rtioxa47 or vehicle injection
#Note: for now I used "Group1" and "Group2" in place of "rtioxa47" and "vehicle" because I don't know how to replace 
  #Group1 and Group2 with "rtioxa47" and "vehicle" without introducing bias
  pattern <- c("Group2", "Group2", "Group1", "Group1")
  
  echomri_data_4 <- echomri_data_2 %>%
    group_by(diet_group) %>%
    arrange(desc(adiposity_index), .by_group = TRUE) %>%
  mutate(treatment = c("Group1", rep(pattern, length.out = n() - 1)))
#-----------------------------------------
#Verify that the assignment method produced homogeneity between treatment groups, within diet groups
  # 1) Create two data frames each with data from just one diet_group (AD_LIBITUM or RESTRICTED)
  # 2) Conduct t tests -> AD_LIB/Group1 vs. AD_LIB/Group2; RESTRICTED/Group1 vs RESTRICTED/Group2
  # 3) Visualize AI for sorted groups -> AD_LIB/Group1 vs. AD_LIB/Group2; RESTRICTED/Group1 vs RESTRICTED/Group2
  
#Create data frame with just AD_LIB mice (n=10)
AD_LIB_only <- echomri_data_4 %>%
    ungroup() %>%
    filter(diet_group== "AD_LIB")
  
#T test: Within the AD_LIB diet_group check if the mean AI differs between Group1 and Group2 (p value)
  t_test_result_AD_LIB <- lm(adiposity_index ~ treatment, data = AD_LIB_only)
  # View the result
  print(t_test_result_AD_LIB)
  # Summary
  summary(t_test_result_AD_LIB)
  
#Create data frame with just RESTRICTED NZO mice (n=12)
  RESTRICTED_only <- echomri_data_4 %>%
    ungroup() %>%
    filter(diet_group== "RESTRICTED")
  
#T test: Within the RESTRICTED diet_group check if the mean AI differs between Group1 and Group2 (p value)
  t_test_result_RESTRICTED <- lm(adiposity_index ~ treatment, data = RESTRICTED_only)
  # View the result
  print(t_test_result_RESTRICTED)
  # Summary
  summary(t_test_result_RESTRICTED)
  
#Visualize distribution of assigned treatment groups for AD LIBITUM diet_group 
  echomri_data_3_adlib <-AD_LIB_only %>%
    ggplot(aes(x = treatment, y = adiposity_index, color = ID)) +
    geom_point(size = 3, alpha = 0.8) + # Individual points
    geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +    # Labels for each point
    stat_summary(
      fun.data = mean_se,
      geom = "pointrange",
      size = 1.5,
      shape = 21,
      color = "black",
      fill = "red",
      position = position_dodge(width = 0.2)     # Mean and SEM as big points
    ) +
    theme_minimal() +     # Aesthetic improvements
    labs(
      title = "Adiposity Index by Treatment Group (Ad Libitum)",
      x = "Treatment",
      y = "Adiposity Index (fat/lean mass)",
      color = "Animal ID"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  echomri_data_3_adlib
  
#Visualize distribution of assigned treatment groups for RESTRICTED diet_group 
  echomri_data_3_restricted <-RESTRICTED_only %>%
    ggplot(aes(x = treatment, y = adiposity_index, color = ID)) +
    geom_point(size = 3, alpha = 0.8) + # Individual points
    geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +    # Labels for each point
    stat_summary(
      fun.data = mean_se,
      geom = "pointrange",
      size = 1.5,
      shape = 21,
      color = "black",
      fill = "red",
      position = position_dodge(width = 0.2)     # Mean and SEM as big points
    ) +
    theme_minimal() +     # Aesthetic improvements
    labs(
      title = "Adiposity Index by Treatment Group (Food restricted)",
      x = "Treatment",
      y = "Adiposity Index (fat mass/lean mass)",
      color = "Animal ID"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  echomri_data_3_restricted
