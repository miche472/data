#Cohort 3-5 NZO mice (Nov. 2024-July 2025)

#MILLIPLIEX (read on LUMINEX plate reader) data -- run 5/7/2025 ####
#Assay conducted by Lauren Michels and Audrey Daugherty in Kotz lab
#Panel of analytes related to metabolism
#Samples taken from NZO mice at peak of obesity and at the end of BW loss
#BW loss achieved via caloric restriction of Research Diets D12450Ki

#Note that we should have used protease inhibitors for GLP-1, 
      #active Ghrelin, glucagon, & amylin. Therefore levels of these analytes 
      #may not be accurate

#started script on 10-1-26
#revised on 10-5-26

##Objectives ####
###1. Evaluate change over time for the same mice ####
###2. Evaluate differences between ad libitum and caloric restriction mice after BW loss ####
#Use linear mixed models to account for multiple measurements from the same mice

#When relevant, remove wells with ND, BLOQ, or ND_BLOQ
#Legend for "Messages":
  #ND = Response is below limit of detection (LoD);     
  #EXT = Using extrapolated values;     
  #BLOQ = Response is below the lower limit of quantification (LLoQ);    
  #Defects = Well(s) have defects

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

#Libraries
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
library(hms)

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#Format plot (LM version 2) 
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
custom_colors <- c("Ad libitum" = "#FAAC41","Restricted" = "#3498DB")
custom_colors2 <- c("Ad libitum" = "#E67E22","Restricted" = "#1d5e8a")

#Read in MILLIPLEX data
MILLIPLEX_readin <- read_csv("~/Documents/GitHub/data/data/MILLIPLEX_L.csv") 

#Format MILLIPLEX data
MILLIPLEX <- MILLIPLEX_readin %>%
mutate(ID = as.factor(ID),
       DATE = lubridate::mdy(DATE),
       DATE_blood = DATE,
       SABLE = factor(SABLE, levels = c("Peak obesity", "BW loss"))
       )


# LMM: Hardcode for single analyte ####

## GIP ####
# gastric inhibitory polypeptide -> glucose-dependent insulinotropic polypeptide
df_GIP <- MILLIPLEX %>%
  ungroup() %>%
  filter(Analyte == "GIP total") %>%
  group_by(ID, SABLE, GROUP) %>%
  drop_na(Avg_Result)

#Build linear mixed model for GIP in serum (derived from retro orbital blood)
model_GIP <- lmer(Avg_Result ~ SABLE * GROUP + (1 | ID), data = df_GIP)
summary(model_GIP)

#Calculate estimated marginal means 
emm_GIP <- emmeans(model_GIP, ~ SABLE * GROUP, cov.reduce = mean)
emm_GIP_df <- as.data.frame(emm_GIP)

# Pairwise contrasts within each GROUP
contrasts_by_group_GIP <- contrast(emm_GIP, method = "pairwise", by = "GROUP")
contrasts_by_group_GIP_df <- as.data.frame(contrasts_by_group_GIP)

# Pairwise contrasts within each SABLE (time point)
contrasts_by_SABLE_GIP <- contrast(emm_GIP, method = "pairwise", by = "SABLE")
contrasts_SABLE_GIP_df <- as.data.frame(contrasts_by_SABLE_GIP)

#Bar plot
Plot_GIP <- ggplot(df_GIP, aes(x=SABLE, y=Avg_Result, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) + 
  scale_color_manual(values = custom_colors2) +
  format.plot_LM2 +
  labs(x="Stage",
       y= "pg/mL",
       title= "GIP",
       color="Treatment", fill="Treatment", group="Treatment")
Plot_GIP



#Box plot
boxPlot_GIP <- ggplot(df_GIP, aes(x = SABLE, y = Avg_Result, 
                      group = interaction(SABLE, GROUP), fill = GROUP, color = GROUP)) +
  geom_boxplot(position = position_dodge(width = 0.8),
    width = 0.65, alpha = 0.4, linewidth = 0.8) +
  geom_point(aes(group = GROUP), position = position_jitterdodge(dodge.width = 0.8,
      jitter.width = 0.2), alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = ID, group = GROUP),
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
    size = 3, color = "black", max.overlaps = Inf, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) + 
  scale_color_manual(values = custom_colors2) +
  format.plot_LM2 +
  labs(x = "Stage",y = "pg/mL", title = "GIP",
    color = "Treatment",
    fill = "Treatment")
boxPlot_GIP


#------------------------------#.
#------------------------------#.
## Leptin ####
df_Leptin_echo <- MILLIPLEX_echo %>%
  ungroup() %>%
  filter(Analyte == "Leptin") %>%
  group_by(ID, SABLE, GROUP) %>%
  drop_na(Avg_Result)

### Leptin: LMM, emmeans, contrasts ####
model_Leptin <- lmer(Avg_Result ~ SABLE * GROUP + (1 | ID), data = df_Leptin_echo)
emm_Leptin <- emmeans(model_Leptin, ~ SABLE * GROUP, cov.reduce = mean) 
emm_Leptin_df <- as.data.frame(emm_Leptin)
contrasts_by_group_Leptin_df <- as.data.frame(contrast(emm_Leptin, 
                                                       method = "pairwise", by = "GROUP"))
contrasts_by_SABLE_Leptin_df <- as.data.frame(contrast(emm_Leptin, 
                                                       method = "pairwise", by = "SABLE"))


###Leptin+Lean: LMM, emmeans, contrasts ####
model_Leptin_lean <- lmer(Avg_Result ~ SABLE * GROUP + Lean + (1 | ID), data = df_Leptin_echo)
emm_Leptin_lean <- emmeans(model_Leptin_lean, ~ SABLE * GROUP, cov.reduce = mean)
emm_Leptin_df_lean <- as.data.frame(emm_Leptin_lean)
contrasts_by_group_Leptin_df_lean <- as.data.frame(contrast(emm_Leptin_lean, 
                                                            method = "pairwise", by = "GROUP"))
contrasts_by_SABLE_Leptin_df_lean <- as.data.frame(contrast(emm_Leptin_lean, 
                                                            method = "pairwise", by = "SABLE"))


###Leptin+Fat: LMM, emmeans, contrasts ####
model_Leptin_fat <- lmer(Avg_Result ~ SABLE * GROUP + Fat + (1 | ID), data = df_Leptin_echo)
emm_Leptin_fat <- emmeans(model_Leptin_fat, ~ SABLE * GROUP, cov.reduce = mean)
emm_Leptin_df_fat <- as.data.frame(emm_Leptin_fat)
contrasts_by_group_Leptin_df_fat <- as.data.frame(contrast(emm_Leptin_fat, 
                                                            method = "pairwise", by = "GROUP"))
contrasts_by_SABLE_Leptin_df_fat <- as.data.frame(contrast(emm_Leptin_fat, 
                                                            method = "pairwise", by = "SABLE"))


###Leptin+BW: LMM, emmeans, contrasts ####
model_Leptin_BW <- lmer(Avg_Result ~ SABLE * GROUP + BW + (1 | ID), data = df_Leptin_echo)
emm_Leptin_BW <- emmeans(model_Leptin_BW, ~ SABLE * GROUP, cov.reduce = mean)
emm_Leptin_df_BW <- as.data.frame(emm_Leptin_BW)
contrasts_by_group_Leptin_df_BW <- as.data.frame(contrast(emm_Leptin_BW, 
                                                          method = "pairwise", by = "GROUP"))
contrasts_by_SABLE_Leptin_df_BW <- as.data.frame(contrast(emm_Leptin_BW, 
                                                          method = "pairwise", by = "SABLE"))

#Bar plot
Plot_Leptin <- ggplot(df_Leptin_echo, aes(x=SABLE, y=Avg_Result, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="black") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) + 
  scale_color_manual(values = custom_colors2) +
  format.plot_LM2 +
  labs(x="Stage",
       y= "pg/mL",
       title= "Leptin",
       color="Treatment", fill="Treatment", group="Treatment")
Plot_Leptin


#Box plot
boxPlot_Leptin <- ggplot(df_Leptin_echo, aes(x = SABLE, y = Avg_Result, 
                      group = interaction(SABLE, GROUP), fill = GROUP, color = GROUP)) +
  geom_boxplot(position = position_dodge(width = 0.8),
    width = 0.65, alpha = 0.4, linewidth = 0.8) +
  geom_point(aes(group = GROUP), position = position_jitterdodge(dodge.width = 0.8,
      jitter.width = 0.2), alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = ID, group = GROUP),
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
    size = 3, color = "black", max.overlaps = Inf, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = custom_colors) + 
  scale_color_manual(values = custom_colors2) +
  format.plot_LM2 +
  labs(x = "Stage",y = "pg/mL", title = "Leptin",
    color = "Treatment",
    fill = "Treatment")
boxPlot_Leptin

#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.


#Check missing data ####
#How many missing Avg Result values for each Analyte?
  #Message = BLOQ -> Avg_Result = value 
  #Message = BLOQ_ND, EXT_BLOQ, EXT_ND, or EXT -> Avg_Result = NA

MILLIPLEX %>%
  group_by(Analyte) %>%
  summarise(
    n_mice = n_distinct(ID[!is.na(Avg_Result)]),
    n_observations = sum(!is.na(Avg_Result))
  )

#Number of IDs with Message = BLOQ for each analyte for at least one SABLE time point
MILLIPLEX %>%
  filter(Message == "BLOQ") %>%
  group_by(Analyte) %>%
  summarise(
    n_ID = n_distinct(ID)
  )

#Total number of Message = BLOQ for a given analyte 
MILLIPLEX %>%
  filter(Message == "BLOQ") %>%
  group_by(Analyte) %>%
  summarise(
    n_ID = n())


#--------------------------------------------------------#.
# LMM: All analytes     ####
#--------------------------------------------------------#.

##1. Function for LMM and emmeans for each analyte ####
# creating a function for LMM and emmeans circumvents the need to replicate the 
#analysis I did with GIP for each analyte (i.e. more parsimonious code)

analyze_analyte <- function(analyte_name) {
  
  df <- MILLIPLEX %>%
    filter(Analyte == analyte_name) %>%
    drop_na(Avg_Result) 
  
  # Run linear mixed model
  model <- tryCatch(lmer(Avg_Result ~ SABLE * GROUP + (1 | ID), data = df),
    error = function(e) NULL)
  
  # If model failed, return failure --> model may fail for an analyte if there are insufficient observations
  if (is.null(model)) {
    return(list(success = FALSE, model = NULL, emm = NULL,
      contrasts_SABLE = NULL, contrasts_GROUP = NULL))
  }
  
  # Estimated marginal means
  emm <- emmeans(model, ~ SABLE * GROUP)
  
  # Pairwise comparisons within each SABLE
  contrasts_SABLE <- contrast(emm, method = "pairwise", by = "SABLE")
  
  # Pairwise comparisons within each GROUP
  contrasts_GROUP <- contrast(emm, method = "pairwise", by = "GROUP")
  
  list(success = TRUE, model = model, emm = emm,
    contrasts_SABLE = contrasts_SABLE, contrasts_GROUP = contrasts_GROUP)
}


## 2. Run the function for each analyte ####
analytes <- unique(MILLIPLEX$Analyte)
results <- lapply(analytes, analyze_analyte)
names(results) <- analytes


## 3. Check which analytes successfully ran -> create df ####
model_status <- data.frame(
  Analyte = analytes,
  Success = sapply(results, function(x) x$success))


## 4. Combine SABLE pairwise comparisons for all analytes into one df ####
contrasts_SABLE_all <- bind_rows(
  lapply(names(results), function(analyte) {
    if (results[[analyte]]$success) {
      as.data.frame(results[[analyte]]$contrasts_SABLE) %>%
        mutate(Analyte = analyte)
      }}))


## 5. Combine GROUP pairwise comparisons for all analytes into one df ####
contrasts_GROUP_all <- bind_rows(
  lapply(names(results), function(analyte) {
    if (results[[analyte]]$success) {
      as.data.frame(results[[analyte]]$contrasts_GROUP) %>%
        mutate(Analyte = analyte)
    }}))


#--------------------------------------------------------#.
# LMM: All analytes (BLOQ removed)        ####
#--------------------------------------------------------#.

## 1. Remove Message = BLOQ rows ####
MILLIPLEX_BLOQ <- MILLIPLEX %>%
  filter(Message != "BLOQ" | is.na(Message))
  #ungroup() %>%
  #group_by(Analyte, SABLE) %>%
  #summarise(n_IDs = n())


## 2. Function to analyze each analyte ####
analyze_analyte_BLOQ <- function(analyte_name) {
  
  df_BLOQ <- MILLIPLEX_BLOQ %>%
    filter(Analyte == analyte_name) %>%
    drop_na(Avg_Result) 
  
  # Run linear mixed model
  model_BLOQ <- tryCatch(lmer(Avg_Result ~ SABLE * GROUP + (1 | ID), data = df_BLOQ),
    error = function(e) NULL)
  
  # If model failed, return failure --> model may fail for an analyte if there are insufficient observations
  if (is.null(model_BLOQ)) {
    return(list(success = FALSE, model_BLOQ = NULL, emm_BLOQ = NULL,
      contrasts_SABLE_BLOQ = NULL, contrasts_GROUP_BLOQ = NULL))
  }
  
  # Estimated marginal means
  emm_BLOQ <- emmeans(model_BLOQ, ~ SABLE * GROUP)
  
  # Pairwise comparisons within each SABLE
  contrasts_SABLE_BLOQ <- contrast(emm_BLOQ, method = "pairwise", by = "SABLE")
  
  # Pairwise comparisons within each GROUP
  contrasts_GROUP_BLOQ <- contrast(emm_BLOQ, method = "pairwise", by = "GROUP")
  
  list(success = TRUE, model_BLOQ = model_BLOQ, emm_BLOQ = emm_BLOQ,
    contrasts_SABLE_BLOQ = contrasts_SABLE_BLOQ, contrasts_GROUP_BLOQ = contrasts_GROUP_BLOQ)
}


## 3. Run the function for each analyte ####
analytes_BLOQ <- unique(MILLIPLEX_BLOQ$Analyte)
results_BLOQ <- lapply(analytes_BLOQ, analyze_analyte_BLOQ)
names(results_BLOQ) <- analytes_BLOQ


## 4. Check which analytes successfully ran -> create df ####
model_status_BLOQ <- data.frame(
  Analyte = analytes_BLOQ,
  Success = sapply(results_BLOQ, function(x) x$success))


## 5. Combine SABLE pairwise comparisons for all analytes into one df ####
contrasts_SABLE_all_BLOQ <- bind_rows(
  lapply(names(results_BLOQ), function(analyte) {
    if (results_BLOQ[[analyte]]$success) {
      as.data.frame(results_BLOQ[[analyte]]$contrasts_SABLE_BLOQ) %>%
        mutate(Analyte = analyte)
      }}))


## 6. Combine GROUP pairwise comparisons for all analytes into one df ####
contrasts_GROUP_all_BLOQ <- bind_rows(
  lapply(names(results_BLOQ), function(analyte) {
    if (results_BLOQ[[analyte]]$success) {
      as.data.frame(results_BLOQ[[analyte]]$contrasts_GROUP_BLOQ) %>%
        mutate(Analyte = analyte)
    }}))


#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.
#Correlation tests ####

##Correlation between fat and leptin level? ####
fat_leptin_cor <- df_Leptin_echo %>%
  group_by(GROUP, SABLE) %>%
  summarise(
    n = sum(complete.cases(Fat, Avg_Result)),
    r = cor(Fat, Avg_Result, method = "pearson", use = "complete.obs"),
    p = cor.test(Fat, Avg_Result, method = "pearson")$p.value,
    .groups = "drop"
  )

fat_leptin_cor

ggplot(df_Leptin_echo, aes(x = Fat, y = Avg_Result, color = interaction(GROUP, SABLE))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Fat mass",
    y = "Leptin",
    color = "GROUP / SABLE"
  ) +
  theme_classic()

#It almost seems like leptin becomes uncoupled from fat mass after weight loss...why could this be?

##Correlation between change in fat and change in leptin level? ####
df_change <- df_Leptin_echo %>%
  group_by(ID) %>%
  arrange(SABLE) %>%
  summarise(
    GROUP = first(GROUP),
    Fat_change = last(Fat) - first(Fat),
    Leptin_change = last(Avg_Result) - first(Avg_Result),
    .groups = "drop"
  )

change_cor <- df_change %>%
  group_by(GROUP) %>%
  summarise(
    n = sum(complete.cases(Fat_change, Leptin_change)),
    r = cor(Fat_change, Leptin_change, method = "pearson", use = "complete.obs"),
    p = cor.test(Fat_change, Leptin_change, method = "pearson")$p.value,
    .groups = "drop"
  )

change_cor
#I think that p and r values produced by this correlation analysis suggest that for ad lib mice a greater 
#change in fat mass is associated with a greater change in leptin. 
#In contrast, for mice that lost weight the change in fat mass is not significantly correlated with the 
#change in leptin level

ggplot(df_change, aes(x = Fat_change, y = Leptin_change, color = GROUP)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Change in fat mass", y = "Change in leptin", color = "GROUP") +
  theme_classic()

#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.
# BW & FI: Cohorts 3-5 ####
#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.

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

#Read in BW and FI 
#Read in BW and filter for cohort 19 (Spring 2026 NZO mice)
BW_345 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT %in% c(3,4,5))

#Read in FI and filter for cohorts 3-5 (NZO mice: Nov 2024 - July 2025)
FI_LM_345 <- read_csv("~/Documents/GitHub/data/data/FI_LM.csv") %>%
  filter(COHORT %in% c(3,4,5))

#Create df with BW and FI
BW_FI_345 <- BW_345 %>% #Join FI and BW
  left_join(
    FI_LM_345 %>% 
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



#----------------------------------------------------------------------------#.
#----------------------------------------------------------------------------#.
# EchoMRI: Cohorts 3-5 ####
#----------------------------------------------------------------------------#.
#----------------------------------------------------------------------------#.

#Read in echomri.csv
echoMRI_readin <- read_csv("~/Documents/GitHub/data/data/echomri.csv") 

echoMRI_345 <- echoMRI_readin %>%
  filter(COHORT %in% c(3,4,5)) %>%
  filter(!(ID %in% c(3710,3712,3715,3720,3727))) %>% #Mice that died (3712,3715) or were not included in assay
  mutate(SABLE = case_when(
                  Date=="2025-02-20" & ID >=3706 & ID <=3729 ~ "Peak obesity",
                  Date=="2025-03-28" & ID >=3706 & ID <=3729 ~ "BW loss")) %>%
  drop_na(SABLE) %>%
  group_by(ID, Date) %>%
  mutate(BW=Fat+Lean) %>%
  select(-(Weight)) %>%
  mutate(ID = as.factor(ID),
         SABLE = factor(SABLE, levels = c("Peak obesity", "BW loss")),
         Date_echo = Date) %>%
  ungroup()

#Join body composition data to MILLIPLEX data
       

MILLIPLEX_echo <- MILLIPLEX %>%
  left_join(
    echoMRI_345 %>% 
      select(ID, Date_echo, Lean, Fat, BW, SABLE), 
    by = c("ID", "SABLE"))

