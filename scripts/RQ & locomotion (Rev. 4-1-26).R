#libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)

#RQ and locomotion (PedMeters)

#Process - PedMeters ####
##------- Get only the PedMeters parameter --> value column has Ped_meters data
filter_ped1 <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("PedMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) 

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
filter_ped2 <- filter_ped1 %>%
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
    #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
    #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
    #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
    #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_ped3 <- filter_ped2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2))

#Check number of mice in filter_ped3
filter_ped3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

# Process - AllMeters ####
##------- Get only the AllMeters parameter --> value column has All_meters data
filter_Alm1 <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) 

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
filter_Alm2 <- filter_Alm1 %>%
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
    #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
    #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
    #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
    #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_Alm3 <- filter_Alm2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2))

#Check number of mice in filter_ped3
filter_Alm3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

# Process - RQ ####
#------- Get only the RQ parameter --> value column has RQ data
filter_RQ1 <-sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain")) %>% 
  filter(grepl("RQ_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))) %>% 
  ungroup() %>% 
  group_by(ID) 

#Get complete day 1 and 2 by explicitly defining the zt_time and SABLE_DAY
filter_RQ2 <- filter_RQ1 %>% 
  #Baseline
  mutate(LM_complete_day = case_when(
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_1" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_2" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_2" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_3" ~ 2,
    #Peak obesity
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_8" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_9" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_9" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_10" ~ 2,
    #BW loss
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_12" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_13" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_13" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_14" ~ 2,
    #BW maintenance
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_16" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_17" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_17" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_18" ~ 2,
    #BW regain
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_20" ~ 1,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_21" ~ 1,
    zt_time %in% c(0, 1, 2, 3) & sable_idx=="SABLE_DAY_21" ~ 2,
    zt_time>3 & zt_time<24 & sable_idx=="SABLE_DAY_22" ~ 2)) %>%
  
  filter(!ID %in% c(3715,3712)) %>% #remove mice that died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #technical issues with Sable cages
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Get only the complete day needed (either day 1 or 2)
filter_RQ3 <- filter_RQ2 %>%
  mutate(use_day = case_when(
    SABLE == "Baseline" & LM_complete_day==1 ~ 1,
    SABLE == "Peak obesity" & LM_complete_day==1 ~ 1,
    SABLE == "BW loss" & ID %in% c(3708, 3710, 3714, 3716, 3726) & LM_complete_day==2 ~ 2,
    SABLE == "BW loss" & ID %in% c(3706, 3707, 3711, 3713, 3719, 3720, 3721, 3722, 3727, 3728, 3729) & LM_complete_day==1 ~ 1,
    SABLE == "BW maintenance" & LM_complete_day==2 ~ 2,
    SABLE == "BW regain" & LM_complete_day==1 ~ 1)) %>%
  filter(use_day %in% c(1,2)) #remove observations not from complete day 1 or 2

#Check number of mice in filter_RQ3
filter_RQ3 %>% 
  group_by(SABLE) %>%
  summarise(n_ID = n_distinct(ID)) #Good, n=16 for all SABLE

#---
#Prep RQ, ped meters, and all meters for joining

filter_ped <- filter_ped3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Ped_meters = value) %>%
  rename(parameter_PedMeters = parameter) %>%
  rename(fix_value_PedMeters = fix_value) 

filter_All_meter <- filter_Alm3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) 

filter_RQ <- filter_RQ3 %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(RQ = value) %>%
  rename(parameter_RQ = parameter) %>%
  rename(fix_value_RQ = fix_value)

#---
#Join filter_ped, filter_All_meter, and filter_RQ into a df called filter_ped_RQ_alm
# Add kcal_hr to filter_locom (by ID, DateTime, sable_idx)
filter_ped_alm_RQ <- filter_RQ %>%
  left_join(
    filter_ped %>% 
      select(Ped_meters, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx")) %>%
  left_join(
    filter_All_meter %>% 
      select(All_meters, ID, DateTime, sable_idx),
    by = c("ID", "DateTime", "sable_idx"))

#-----------------------------------------#
#Locomotion calculations ####
##Ambulation (PedMeters) ####
#Compute ambulatory movement

#Graph Ped_meters to identify if there are recording anomalies
#BW loss
BWloss_filter_ped_alm_RQ <- filter_ped_alm_RQ %>%
  filter(SABLE=="BW loss") %>%
  filter(ID==3729)

ggplot(BWloss_filter_ped_alm_RQ,
       aes(x = DateTime, y = Ped_meters,
           color = ID, fill = ID)) +
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
  labs(x="DateTime",
       y= "Ped_meters",
       title= "BW loss:3729")

#After graphing, I think the biggest issue is that the cumulative counter
#resets when Sable is stopped and started again. 
#Identify reset points as times when Ped_meters > lag(Ped_meters)
Obesity_filter_ped_alm_RQ <- filter_ped_alm_RQ %>%
  filter(SABLE=="Peak obesity")

df_fixed <- Obesity_filter_ped_alm_RQ %>%
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(
    reset = Ped_meters < lag(Ped_meters, default = first(Ped_meters)),
    offset = cumsum(if_else(reset, lag(Ped_meters), 0, missing = 0)),
    Ped_continuous = Ped_meters + offset,
    distance = Ped_continuous - lag(Ped_continuous)) %>%
  
ggplot(df_fixed, aes(DateTime)) +
  geom_line(aes(y = Ped_meters), color = "red") +
  geom_line(aes(y = Ped_continuous), color = "blue") +
  facet_wrap(~ID)

#Now try doing this for all Sable time points in one df
df_fixed_all <- filter_ped_alm_RQ %>%
  arrange(DateTime) %>%
  group_by(ID, SABLE) %>%
  mutate(
    reset = Ped_meters < lag(Ped_meters, default = first(Ped_meters)),
    offset = cumsum(if_else(reset, lag(Ped_meters), 0, missing = 0)),
    Ped_continuous = Ped_meters + offset,
    distance = Ped_continuous - lag(Ped_continuous)) %>%
  drop_na(distance) %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(total_distance_m=sum(distance), .groups="drop")

#Can use this as a check for df_fixed_all--must not run the summarise step though
df_fixed_all_plot <- df_fixed_all %>%
  filter(SABLE=="Peak obesity")

ggplot(df_fixed_all_plot, aes(DateTime)) +
  geom_line(aes(y = Ped_meters), color = "red") +
  geom_line(aes(y = Ped_continuous), color = "blue") +
  facet_wrap(~ID)
    
#Alternative method to account for restarting of recording
#Got this from CS's "NZO_figure8_locomotion_collapsed.R"
Ambulation_df <- filter_ped_alm_RQ %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  arrange(DateTime) %>%
  mutate(distance_m = Ped_meters-lag(Ped_meters)) %>%
  filter(distance_m >=0) %>%
  summarise(total_distance_m = sum(distance_m), .groups="drop") %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))


#Box plot
ggplot(Ambulation_df, aes(x = SABLE, y = total_distance_m, fill = GROUP)) + 
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")
  
##All locomotion: ambulation and fine motor (AllMeters) ####
All_movement_df <- filter_ped_alm_RQ %>%
  ungroup() %>%
  group_by(ID, SABLE, LM_complete_day, hr) %>%
  arrange(DateTime) %>%
  mutate(distance_m = All_meters-lag(All_meters)) %>%
  filter(distance_m >=0) %>%
  summarise(total_distance_m = sum(distance_m), .groups="drop") %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
      ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47")) %>%
  mutate(SABLE = factor(SABLE, levels = c("Baseline","Peak obesity","BW loss","BW maintenance","BW regain")))

#Bar plot: AllMeters
ggplot(All_movement_df,
       aes(x = hr, y = total_distance_m,
           color = GROUP, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.3, color="black") +
  #geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             #alpha = 0.6, size = 2) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="DateTime",
       y= "Total - All locomotion (m)",
       title= "BW loss, All locomotion") +
  facet_wrap(~SABLE)

#Box plot: AllMeters
ggplot(All_movement_df, aes(x = SABLE, y = total_distance_m, fill = GROUP)) + 
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white") +
  geom_jitter(position = position_jitter(0.2))

#Linear mixed model 


#-----------------------------------------#
#RQ calculations####
#Calculate RQ on an hourly basis
df_RQ_hour <- filter_ped_RQ %>%
  filter(RQ>0.5 & RQ<1.2) %>%
  group_by(ID, SABLE, hr, GROUP) %>%
  arrange(DateTime) %>%
  summarise(RQ_hour = mean(RQ))

#---Graph raw average RQ during each hour for each GROUP
ggplot(df_RQ_hour, aes(x = hr, y = RQ_hour, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  facet_wrap(~SABLE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title="RQ--Raw",
    x = "Hour",
    y = "RQ",
    fill = "Treatment group")

#Calculate RQ by photo period
df_RQ_photo <- filter_ped_RQ %>%
  filter(RQ>0.5 & RQ<1.2) %>%
  group_by(ID, SABLE, lights, GROUP) %>%
  arrange(DateTime) %>%
  summarise(RQ_photo = mean(RQ))

#---Graph raw average RQ during each photo period for each GROUP
ggplot(df_RQ_photo, aes(x = lights, y = RQ_photo, fill = GROUP)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)) +
  format.plot+
  facet_wrap(~SABLE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title="RQ, photo period--Raw",
    x = "Photo period",
    y = "RQ",
    fill = "Treatment group")

# RQ: Left off at version 2, belowhere...####
#trying to replicate what is in the consensus paper

#Graph using the approach described (smoothing/rolling average and ribbon with s.e.m.)
#Start with just BW loss time point

##Version 1: RQ ####
# 1. Create continuous hour variable
#I think this is necessary so that I don't have an arbitrary hourly bin that would interfere with the 
#calculation of rolling average -> existing variable, hr, is not continuous

filter_ped_BWloss <- filter_ped_RQ %>%
  filter(SABLE == "BW loss") %>%
  arrange(ID, DateTime)

filter_ped_BWloss <- filter_ped_RQ %>%
  arrange(ID, DateTime) %>%
  group_by(ID) %>%
  mutate(
    hour = as.numeric(difftime(DateTime, min(DateTime), units = "hours"))
  ) %>%
  ungroup()
# 2. Smooth within each mouse
# Use k = 30 for true 30-min smoothing with 1-min data
df_smooth_BWloss <- filter_ped_BWloss %>%
  group_by(ID, GROUP) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(
    RQ_smooth = zoo::rollmean(RQ, k = 30, fill = NA, align = "center")
  ) %>%
  ungroup()

# 3. Compute mouse-level means at each time point
df_mouse_BWloss <- df_smooth_BWloss %>%
  group_by(GROUP, zt_time, ID) %>%
  summarise(
    RQ_mouse = mean(RQ_smooth, na.rm = TRUE),
    .groups = "drop")

# 4. Compute group mean ± SEM
df_summary_BWloss <- df_mouse_BWloss %>%
  group_by(GROUP, zt_time) %>%
  summarise(
    mean_RQ_GROUP = mean(RQ_mouse, na.rm = TRUE),
    sem_RQ_GROUP  = sd(RQ_mouse, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# 5. Plot
ggplot(df_summary_BWloss, aes(x = zt_time, y = mean_RQ_GROUP, color = GROUP, fill = GROUP)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_RQ_GROUP - sem_RQ_GROUP,
                  ymax = mean_RQ_GROUP + sem_RQ_GROUP),
              alpha = 0.25,
              color = NA) +
  scale_x_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24)) +
  labs(
    x = "Time (zt_time, hours)",
    y = "Respiratory exchange ratio (RER)",
    title = "24-hour RER (30-min Smoothed)",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

summary(df_smooth_BWloss$zt_time)


##Version 2: RQ ####
#Actually I think if I use this method i can avoid arbitrary binning while still maintaining clock/circadian time

# 1. Smooth within each mouse
df_smooth <- filter_ped_RQ %>%
  filter(SABLE == "BW loss") %>%
  arrange(ID, DateTime) %>%
  group_by(ID, GROUP) %>%
  mutate(
    RQ_smooth = zoo::rollmean(RQ, k = 30, fill = NA, align = "center")) %>%
  ungroup()

#2. Use continuous time (zt_time): zt_time is currently integers --> make fractional
df_smooth <- df_smooth %>%
  mutate(
    zt_time_cont = zt_time + lubridate::minute(DateTime)/60)

#3. Compute mouse-level values at EACH time point
#this step replaces binning (verify that ChatGPT is correct in this statement)
df_mouse <- df_smooth %>%
  group_by(GROUP, ID, zt_time_cont) %>%
  summarise(
    RQ_mouse = mean(RQ_smooth, na.rm = TRUE),
    .groups = "drop")

# 4. Compute group mean ± SEM
df_summary <- df_mouse %>%
  group_by(GROUP, zt_time_cont) %>%
  summarise(
    mean_RQ = mean(RQ_mouse, na.rm = TRUE),
    sem_RQ  = sd(RQ_mouse, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

#5 back convert zt_time to clock time
df_summary2 <- df_summary %>%
  mutate(
    clock_time = if_else(
      zt_time_cont < 4,
      zt_time_cont + 20,
      zt_time_cont - 4
    )
  )
#5. Plot (continuous, smooth curve)
ggplot(df_summary2,
       aes(x = clock_time, y = mean_RQ,
           color = GROUP, fill = GROUP)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = mean_RQ - sem_RQ,
                  ymax = mean_RQ + sem_RQ),
              alpha = 0.25,
              color = NA) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 4),
    labels = function(x) sprintf("%02d:00", x),
    limits = c(0, 24)
  ) +
  labs(
    x = "Clock Time",
    y = "Respiratory Quotient (RQ)",
    title = "RQ during BW Loss (30-min Rolling Average)"
  ) +
  theme_classic(base_size = 14)
