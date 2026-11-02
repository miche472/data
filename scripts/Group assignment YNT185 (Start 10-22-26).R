#Group assignment for OXA dose response study (Oct/Nov 2026) ####

#we have 8 NZO female mice that will be used to pilot the use of voluntary
#oral consumption of YNT-185*2HCl in a strawberry gummy
#Mice will all receive one vehicle gummy and one gummy with each of 3 doses 
#There is a 48hr washout period between doses
#Mice will be recording in Sable for the duration of the study

#This script is based on "Randomization OXA dose (Rev 8-27-26).R"
#Started: 10-22-26
#Revised:

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand
)

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

#-----------------------------------------------------.#
#Calculate doses based on dose for a mouse of average BW 
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


# Identify average BW ####
Group_avg_BW <- BW_FI_19 %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3741)) %>%
  group_by(ID) %>%
  arrange(DATE) %>%
    slice_tail(n = 1) %>%
  ungroup() %>%
  summarise(avg_BW_G = mean(BW)) %>%
  mutate(avg_BW_KG = avg_BW_G/1000)
Group_avg_BW

#Calculate doses for a mouse of average BW
#30mg/kg, 60mg/kg, 100mg/kg

#If 100mg per kg, then "x" mg per Group_avg_BW
Dose_mass_G <- Group_avg_BW %>%
  mutate(Mass_100mgkg_G = avg_BW_KG *100,
         Mass_60mgkg_G = avg_BW_KG *60,
         Mass_30mgkg_G = avg_BW_KG *30)

#-----------------------------------------------------.#
#-----------------------------------------------------.#
#-----------------------------------------------------.#

set.seed(123)

#-----------------------------------------------------.#
## 1. Establish variables for IDs and treatments ####
#-----------------------------------------------------.#

# Placeholder IDs and Doses
ID <- c("S1","S2","S3","S4","S5","S6","S7","S8")

DOSE <- c("Dose A",
          "Dose B",
          "Dose C",
          "Dose D")

# Actual mouse IDs and dose values
actual_ID <- c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3741)

actual_DOSE <- c(
  "Vehicle",
  "1.8mg",
  "3.5mg",
  "5.8mg")

#-----------------------------------------------------.#
## 2. Build standard 4 × 4 Latin square ####
#-----------------------------------------------------.#
latin_square <- tribble(
  ~INJECT_DAY_1, ~INJECT_DAY_2, ~INJECT_DAY_3, ~INJECT_DAY_4,
  DOSE[1],  DOSE[2],  DOSE[3],  DOSE[4],
  DOSE[2],  DOSE[3],  DOSE[4],  DOSE[1],
  DOSE[3],  DOSE[4],  DOSE[1],  DOSE[2],
  DOSE[4],  DOSE[1],  DOSE[2],  DOSE[3])


# Randomize order of the three sequences
latin_square <- latin_square %>%
  slice_sample(n = 4)

#-----------------------------------------------------.#
## 3. Replicated latin square: 8 mice & 4 treatments ####
#-----------------------------------------------------.#
# Shuffle the eight sequences and assign IDs
design <- bind_rows(latin_square, latin_square) %>%
  slice_sample(n = 8) %>%
  mutate(ID = sample(ID), .before = 1)
design

#-----------------------------------------------------.#
## 4. Replace placeholder variable names ####
#-----------------------------------------------------.#
# Replace placeholder ID and DOSE with study specific values
id_key <- setNames(sample(actual_ID), ID)
dose_key <- setNames(actual_DOSE, DOSE)

design_final <- design %>%
  mutate(
    ID = id_key[ID],
    across(starts_with("INJECT"), ~ dose_key[.x]))
design_final

write_csv(x = design_final, "../data/design_final_10-22-26.csv")
#write_csv(x = design_final, "../data/design_final_10-28-26.csv")

#-----------------------------------------------------.#
## 5. Change to long format #### 
#-----------------------------------------------------.#
#Allows for later joining with other data frames
# Convert to long format --> One row per mouse × dose combination

design_long <- design_final %>%
  pivot_longer(
    cols = starts_with("INJECT_DAY"),
    names_to = "INJECT_DAY",
    values_to = "DOSE") %>%
  mutate(PERIOD = readr::parse_number(INJECT_DAY)) %>%
  arrange(ID, INJECT_DAY)
design_long

write_csv(x = design_final, "../data/design_long_10-22-26.csv")
#write_csv(x = design_final, "../data/design_long_10-28-26.csv")