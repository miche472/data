#Orexin A icv dose response (July/August 2026)

#Started: 7-30-26
#Revised: 8-3-26 (simplified approach from 7-31-26)

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


 
 
 
 
#Calculate time after injection (in minutes and in hours)
 #Make a column for whether an observation is within 24hrs of an injection --> this is a complete day within my experimental paradigm
OXA_Sable_joined_2 <- OXA_Sable_joined %>%
  mutate(minutes_post = as.numeric(difftime(DateTime,INJECTION_DateTime,units = "mins"))) %>%
  mutate(hrs_post = minutes_post/60) %>%
  mutate(In_2hr = if_else((minutes_post<=120), 1, 0),
        In_4hr = if_else((minutes_post<=240), 1, 0),
        In_24hr = if_else((minutes_post<=1440), 1, 0))
 
 


#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
 
 #Filter Sable data
sable_dwn_OXA <- sable_dwn %>%
  filter(COHORT==19) %>%
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup() %>%
  #Need to get rid of sable data for 3731 after the invalid injection. 
  group_by(ID, DateTime) %>%
  arrange_by(DateTime) %>%
  filter(!(ID= "3731" & DateTime >= "2026-07-30 17:33:42" & DateTime < "2026-08-03 11:33:42")) %>%
  ungroup()
  
#Read in meta data for injection time
read_injection_time <- read_csv("../data/META_INJECTIONS_LM.csv")

#Process injection time data
injection_time <- read_injection_time %>%
  mutate(INJECTION_DateTime = lubridate::mdy_hm(INJECTION_TIME),
         ID=as.factor(ID)) %>%
  filter(VALID == "VALID") %>% #removes INVALID injection for 3731
  arrange(ID, INJECTION_DateTime) %>%
  group_by(ID) %>%
  mutate(PERIOD = row_number()) %>%
  ungroup()

#Identify the exact injection date/time which corresponds to each chunk of recording (i.e. Period of recording)
#This data frame has one row for each DateTime +ID...it does not yet account for the fact that multiple parameters (10)
#were measured at once for each mouse
 injection_assignment <- sable_dwn_OXA %>%
  distinct(ID, DateTime) %>%
  left_join(
    injection_time,
    by = join_by(
      ID,
      DateTime >= INJECTION_DateTime
    ),
    relationship = "many-to-many"
  ) %>%
  group_by(ID, DateTime) %>%
  slice_max(
    INJECTION_DateTime,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
 
#Join 
OXA_Sable_joined <- sable_dwn_OXA %>%
  left_join(
    injection_assignment,
    by = c("ID", "DateTime")
  )

