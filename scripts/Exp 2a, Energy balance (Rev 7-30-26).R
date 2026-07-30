#combine FI, BW and EE for experiment 2 mice

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

# code to create BW and FI df called "BW_FI_19" ####
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


# Code to calculate Energy expenditure data (df called Minute_sum_EE_hr_alldays) ####
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

#Get just locomotion for sable_dwn
filter_loc1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(zt_time = zt_time(DateTime),
    is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
    complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, SABLE, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                             min(DateTime),
                                             units="hours")),
    zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
    zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()


# Get just TEE from sable_dwn
filter_TEE1_revised <-sable_dwn %>%
  filter(COHORT ==19) %>%  
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3747,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10") ~ "BW loss",
    ID %in% c(3744,3745,3746,3748,3749,3752) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
                       "SABLE_DAY_8") ~ "BW loss")) %>%
  filter(grepl("kcal_hr_*", parameter)) %>%
  ungroup() %>% 
  arrange(ID, SABLE, DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(zt_time = zt_time(DateTime),
         is_zt_init = replace_na(as.numeric(zt_time == 0 & lag(zt_time, default = -1) != 0), 0),
         complete_days = cumsum(is_zt_init)) %>%
  filter(complete_days > 0) %>%
  group_by(ID, SABLE, complete_days) %>%
  mutate(recording_duration = as.numeric(difftime(max(DateTime),
                                                  min(DateTime),
                                                  units="hours")),
         zt_hours = n_distinct(zt_time)) %>%
  filter(recording_duration >= 22.5, #complete day must have a total of at least 22.5 hrs of data
         zt_hours >= 22) %>% #complete day must have data from at least 22 of the zt hours
  ungroup()


#---
#In df filter_loc2 use mutate to make a column called AllMeters_ using the value column data
#In df filter_TEE2 use mutate to make a column called kcal_hr_ using the value column data

filter_loc2 <- filter_loc1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_TEE2 <- filter_TEE1_revised %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)

#---
#Join filter_loc2 and filter_TEE2 into a df called filter_loc_TEE2
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_loc_TEE2 <- filter_loc2 %>%
  left_join(
    filter_TEE2 %>% 
      select(Kcal_Hr, ID, DateTime, sable_idx), # Should this be changed now that I have defined days differently? maybe use complete_days and SABLE rather than sable_idx?
    by = c("ID", "DateTime", "sable_idx")) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss"))) %>%
  filter(!ID %in% c(3748, 3749, 3751)) #issues with recording


#On some days mouse A was recording in cage I in the morning and mouse B starting recording in cage I
#in the evening. META.csv only indicates dates not hours, so I need to extricate the data from these two mice
#in the code
new_filter_loc_TEE2 <- filter_loc_TEE2 %>%
  group_by(ID, SABLE) %>%
  filter(max(complete_days) <= 2 |complete_days != max(complete_days)) %>% 
  ungroup()
#Now I should have accurate data for days when two mice record in the same cage number

#Part 1 of 2: TRY associating a date with complete days to facilitate pairing of manual BW/FI with Sable data####
Complete_days <- new_filter_loc_TEE2 %>%
  group_by(ID, SABLE, complete_days) %>%
  summarise(
    start_time = min(DateTime),
    end_time   = max(DateTime),
    start_date = as.Date(start_time),
    end_date   = as.Date(end_time),
    .groups = "drop"
  )
#end ####


#Calculate EE for EACH complete day #### 
Minute_sum_EE_hr_alldays <- new_filter_loc_TEE2 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, complete_days) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, complete_days, hr, lights) %>% #added lights as an additional grouping factor
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    # Total Energy Expenditure (TEE) --> minute by minute summation
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    # Resting metabolic rate for each hour calculated using 10th percentile method for RMR:
    #Note: this is the rate of RMR rate within each hr not for the entire 24hrs
    #Calculating this will allow for calculation of NEAT
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    # Total RMR energy across observed time in the hour
    RMR_kcal = RMR_rate * (n_obs / 60),
    # NEAT= TEE-RMR when the mouse is moving
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE), .groups = "drop") %>%
  
  # Verify that TEE = RMR + NEAT --> TEEvsRMR_NEAT should be close to zero
  mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal)) %>%
  #Re-attach DRUG
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Peak obesity","BW loss")))

#Part 2 of 2: TRY associating a date with complete days to facilitate pairing of manual BW/FI with Sable data####
#Join df= Complete_days to df= Minute_sum_EE_hr_alldays
Complete_days_Minute_sum_EE_hr <- Minute_sum_EE_hr_alldays %>%
  left_join(
    Complete_days %>% 
      select(ID, SABLE, complete_days, end_date), 
    by = c("ID", "SABLE", "complete_days"))


# end ####

# Start delete ? ####
#Calculate daily sum (add all hours together)
Daily_EE_alldays <- Minute_sum_EE_hr_alldays %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG, complete_days) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)))
# end delete ? ####

#---
#---
#---
# Join df for EE --> Minute_sum_EE_hr_alldays AND df for BW/FI --> BW_FI_19
#Want to keep all EE values

#change name of DATE in df=BW_FI_19 to end_date
BW_FI_19_date <- BW_FI_19 %>%
  group_by(ID, DATE) %>%
  mutate(end_date=DATE) %>%
  ungroup()

BW_FI_19_EE <- Complete_days_Minute_sum_EE_hr %>%
  left_join(
    BW_FI_19_date %>% 
      select(ID, end_date, corrected_intake_gr, corrected_intake_kcal, BW), 
    by = c("ID", "end_date")) 
  #filter(SABLE=="BW loss") The NA values go away when I only look at BW loss mice...I didn't measure BW and FI every day at peak obesity but I did during BW loss

#Now calculate EE and RMR

#Calculate daily sum (add all hours together)
Daily_EE_alldays_BW_FI <- BW_FI_19_EE %>%
  ungroup() %>%
  group_by(ID, SABLE, DRUG, complete_days) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day)),
            corrected_intake_gr = first(corrected_intake_gr),
            corrected_intake_kcal = first(corrected_intake_kcal),
            BW = first(BW),
            end_date = first(end_date))

#I will need to find a different way to attach FI during peak obesity, for now 
#just look at energy balance for the two drug groups during BW loss

# Energy balance = energy intake-energy expenditure

Energy_balance <- Daily_EE_alldays_BW_FI %>%
  ungroup() %>%
  group_by(ID, SABLE, complete_days) %>%
  mutate(E_balance_kcal = corrected_intake_kcal - TEE_kcal_day) 

#Graph energy balance by drug group ####
#graph EE vs BW and FI vs BW

BW_loss_balance <- Energy_balance %>%
  filter(SABLE == "BW loss") %>%
  filter(complete_days %in% c(1,2))

ggplot(BW_loss_balance, aes(x=BW , y=TEE_kcal_day, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values = custom_colors1_GLP_Exp2) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  geom_hline(yintercept=0, linetype="solid", linewidth=1, color="#62748E")+ #start BW loss
  labs(x="Body weight",
       y= "Daily FI (kcal)",
       title= "FI vs BW",
       color="Treatment", fill="Treatment")

#Plots of energy balance: Exp 2A, SABLE= BW loss ####

#Prepare df for plotting
BW_loss_balance <- Energy_balance %>%
  filter(SABLE == "BW loss") %>%
  #filter(complete_days %in% c(1,2)) %>%
  drop_na(E_balance_kcal)

##Bar plot of energy balance: Exp 2A, SABLE= BW loss ####
ggplot(BW_loss_balance, aes(x=complete_days, y=E_balance_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
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
  geom_hline(yintercept=0, linetype="solid", linewidth=1, color="#62748E")+ 
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  labs(x="Day",
       y= "Energy balance",
       title= "Energy balance",
       color="Treatment", fill="Treatment")

##Box plot of energy balance: Exp 2A, SABLE= BW loss ####
ggplot(BW_loss_balance, aes(x = factor(complete_days), y= E_balance_kcal,fill = DRUG,color = DRUG)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.2),alpha = 0.6,size = 2) +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 0,linetype = "solid",linewidth = 1,color = "#62748E") +
  scale_fill_manual(values = custom_colors2_GLP_Exp2) +
  scale_color_manual(values = custom_colors3_GLP_Exp2) +
  format.plot_LM2 +
  labs(x = "Day", y = "Energy balance (kcal)",
    title = "Energy balance during BW loss",
    color = "Treatment",
    fill = "Treatment")