#Orexin A icv dose response (July/August 2026)

#Started: 7-30-26
#Revised: 7-31-26 (see section below with simplified approach)

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

#Libraries####
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
  mmand
)

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

#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#Functions####
#This method gives integer values for zt_time and accounts for the 6:30 to 18:30 cycle
#It is more explicit in its format/logic compared to original zt_time function
zt_time <- function(DateTime){
  time <- hour(DateTime) + minute(DateTime)/60
  zt <- time - 18.5
  zt <- if_else(zt < 0, zt + 24, zt)
  floor(zt)}

# 0.5-2.5hrs after injection...cumulative
# 0.5-4.5hrs after injection...cumulative 
# 0.5-24.5hrs after injection...cumulative

#TEE, RMR and SPA
#Calculate hourly basis --> use photoperiod code from Exp 2a Photo EE & RQ

#Locomotion: time spent moving and distance traveled
#Look at resulting formatting and decide how to append dose information. Could hardcode, but there may be a better approach

# For group assignment I should read in the .csv file that I created in the group assignment script
# use the one from 7-17 for all mice except 3738. Use the one from 7-30 for 3738

Group_assignment <- read_csv("../data/design_final_7-17.csv") 
Group_assignment_3738 <- read_csv("../data/design_final_7-30.csv") %>%
  filter(ID==3738)
Group_assignment_OXA <- bind_rows(Group_assignment, Group_assignment_3738) %>%
   mutate(ID=as.factor(ID))

Group_assignment_OXA_long <- Group_assignment_OXA %>%
  pivot_longer(
    cols = starts_with("INJECT_DAY"),
    names_to = "INJECT_DAY",
    values_to = "DOSE") %>%
  mutate(PERIOD = readr::parse_number(INJECT_DAY)) %>%
  arrange(ID, INJECT_DAY)
Group_assignment_OXA_long

#Associate calendar date with each INJECT_DAY
Group_assignment_OXA_long <- Group_assignment_OXA_long %>%
mutate(
  date= case_when(
  INJECT_DAY == "INJECT_DAY_1" ~ ymd("2026-07-18"),
  INJECT_DAY == "INJECT_DAY_2" ~ ymd("2026-07-21"),
  INJECT_DAY == "INJECT_DAY_3" ~ ymd("2026-07-24"),
  INJECT_DAY == "INJECT_DAY_4" ~ ymd("2026-07-27"),
  INJECT_DAY == "INJECT_DAY_5" ~ ymd("2026-07-30"))) %>%
mutate(DRUG= case_when(date>"2026-07-04" & date <"2026-08-03" ~ "Orexin A")) %>%
  ungroup()

#Filter Sable data
sable_dwn_OXA <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()

#Join Sable data with group assignment data
OXA_Sable <- sable_dwn_OXA %>% #Join FI and BW
  left_join(
    Group_assignment_OXA_long %>% 
      select(ID, INJECT_DAY, DOSE, PERIOD, date, DRUG),
    by = c("ID", "date"))

#Add exact injection timestamp using META_INJECTIONS.csv
injection_time <- read_csv("../data/META_INJECTIONS.csv") %>%
  select(!("TEMPERATURE_F")) %>%
  drop_na(VALID) %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME)) 

injection_time_PERIOD <- injection_time %>%
   mutate(ID=as.factor(ID)) %>%
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = dense_rank(as.Date(INJECTION_DateTime))) %>%
  ungroup()
 

#Combine OXA_Sable with injection_time_PERIOD and create a new column indicating the injection date/time for each Period (ie injection #)
 OXA_SABLE_injections <- OXA_Sable %>% 
  left_join(
    injection_time_PERIOD %>% 
      select(ID, INJECTION_TIME, VALID, INJECTION_DateTime, PERIOD),
    by = c("ID", "PERIOD"))
 
 
 
 
 #Possible alternative approach --> simplified####
 #Instead of using the data frame created by the group assignation script, extract dose, drug, period, date of inejction
 #etc. from the injection_metadata (this requires adding a column called DOSE to the metadata)
 
 #Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
 #Filter Sable data
sable_dwn_OXA <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup()
 
 #Bring in drug, dose, order of injection, DateTime of injections
injection_time <- read_csv("../data/META_INJECTIONS_LM.csv") %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  #arrange(ID, INJECTION_DateTime) %>%
  #group_by(ID) %>%
  #mutate(PERIOD = dense_rank(as.Date(INJECTION_DateTime))) %>%
  #ungroup()
 
 
#Use a rolling join to bring together sable data and injection data.
#Idea is to assign the date/time of the injection after which each row was recorded
 OXA_Sable_joined <- sable_dwn_OXA %>%
  left_join(
    injection_time,
    by = "ID" #it is expected that R gives a message about a many to many relationship. If I want to silence this I can add "relationship = "many-to-many"
  ) %>%
  filter(DateTime >= INJECTION_DateTime) %>%
  group_by(ID, DateTime) %>%
  slice_max(INJECTION_DateTime, n = 1) %>%
  ungroup()
 
 
 #Check: confirm that every row got assigned an injection
 OXA_Sable_joined %>%
  summarise(missing_injection = sum(is.na(INJECTION_DateTime)), total_rows = n())
 #All rows got an injection time, as intended
 
#Calculate time after injection (in minutes and in hours)
 #Make a column for whether an observation is within 24hrs of an injection --> this is a complete day within my experimental paradigm
OXA_Sable_joined_2 <- OXA_Sable_joined %>%
  mutate(minutes_post = as.numeric(difftime(DateTime,INJECTION_DateTime,units = "mins"))) %>%
  mutate(hrs_post = minutes_post/60) %>%
  mutate(In_2hr = if_else((minutes_post<=120), 1, 0),
        In_4hr = if_else((minutes_post<=240), 1, 0),
        In_24hr = if_else((minutes_post<=1440), 1, 0))
 
 
 
 
