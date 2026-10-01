#Cohort 3-5 NZO mice (Nov. 2024-July 2025)

#MILLIPLIEX (read on LUMINEX plate reader) data -- run 5/7/2025 ####
#Assay conducted by Lauren Michels and Audrey Daugherty in Kotz lab
#Panel of analytes related to metabolism
#Samples taken from NZO mice at peak of obesity and at the end of BW loss
#BW loss achieved via caloric restriction of Research Diets D12450Ki

#started script on 10-1-26
#revised on __

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

# Define custom colors
custom_colors <- c("Ad libitum" = "#FAAC41","Restricted" = "#3498DB")
custom_colors2 <- c("Ad libitum" = "#E67E22","Restricted" = "#1d5e8a")

#Read in MILLIPLEX data
MILLIPLEX_readin <- read_csv("~/Documents/GitHub/data/data/MILLIPLEX_L.csv") 

#Format MILLIPLEX data
MILLIPLEX <- MILLIPLEX_readin %>%
mutate(ID = as.factor(ID),
       DATE = lubridate::mdy(DATE),
       SABLE = factor(SABLE, levels = c("Peak obesity", "BW loss"))
       )


# Analyze Milliplex results ####

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




#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.

# Create BW and FI data frame for cohorts 3, 4, 5 ####
## COHORT_19 BW and FI ####

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

#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.
#----------------------------------------------------------------------------.

