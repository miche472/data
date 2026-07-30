#Orexin A icv dose response (July/August 2026)

#Started: 7-30-26
#Revised: 

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

# FOr group assignment I should read in the .csv file that i created in the group assignment script
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
mutate(date= case_when(
  INJECT_DAY == "INJECT_DAY_1" ~ "2026-07-18",
  INJECT_DAY == "INJECT_DAY_2" ~ "2026-07-21",
  INJECT_DAY == "INJECT_DAY_3" ~ "2026-07-24",
  INJECT_DAY == "INJECT_DAY_4" ~ "2026-07-27",
  INJECT_DAY == "INJECT_DAY_5" ~ "2026-07-30"),
  DRUG= case_when(DATE>="2026-07-18" & DATE <="2026-07-30" ~ "Orexin A")) %>%
  mutate(date=as.date(date))

#Join Sable data with 
sable_dwn_Coh_19 <- sable_dwn %>%
  filter(COHORT==19)

OXA_Sable <- sable_dwn_Coh_19 %>% #Join FI and BW
  left_join(
    Group_assignment_OXA_long %>% 
      select(ID, INJECT_DAY, DOSE, PERIOD, date, DRUG),
    by = c("ID", "date"))
