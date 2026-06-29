####Objective is to graph distribution of glucose measurements in NZO mice####
####date is March 25 2025

# libs ----

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand
)

#Read in cohorts 3 data

Cohort3 <- read_csv("../data/COHORT_3.csv") 

Cohort3T <- Cohort3 %>%
  drop_na(FASTED_GLU_mg_dL) %>%
  group_by(ID)

####This works! It deletes rows where both fasted and basal glucose are NA####  
Cohort3D <- Cohort3 %>%
  ungroup() %>%
  filter(FASTED_GLU_mg_dL>0 | BASAL_GLU_mg_dL>0) %>%
  mutate(glucose_time = if_else(BASAL_GLU_mg_dL = NULL, "Basal", "Fasted")) %>%
  #mutate(glucose_time = if_else(FASTED_GLU_mg_dL > 0, "Fasted", "Basal"))

####Combine Basal and Fasted glucose into one column (this also deletes each of the columns...can I retain original?)####
Cohort3_glucose <- Cohort3D %>%
  unite(Glucose_mg_dL, BASAL_GLU_mg_dL, FASTED_GLU_mg_dL, sep = ", ", na.rm = TRUE)

#Graph blood glucose distribution
  ggplot(Cohort3_glucose, aes(x = DATE, y = Glucose_mg_dL, group = ID)) + 
    geom_point() + geom_jitter(position=position_jitter(0.2)) +
    theme(legend.position = "right") + labs(x=NULL) +labs(y= "Blood glucose (mg/dL)")
  