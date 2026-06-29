#Incretin Agonist optimization study script

#Started:5-7-26
#Revised:5-8-26

#Objectives:
  #1. Create a df that includes naive mice as a 2nd age-matched control group
      #in addition to the vehicle control group 
  #2. Monitor body comp. (echo data) throughout the weight loss study

#Format plot
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black")) # keep axis lines
# Define custom colors
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb","Survodutide" = "#c93618", "Vehicle" = "#403d3c", "Naive" ="darkgray")
custom_colors_GLPB <- c("Tirzepatide" = "lightblue","Survodutide" = "lightgreen", "Vehicle" = "darkgray")


# libraries
library(mmand)
library(pacman)
library(this.path)
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

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)
#---
#Objective 1: 2nd control group ####
#Need to create a variable for naive mice called matched_days --> pretend
#like the day after their first echoMRI is their first day of treatment 
#For mice undergoing weight loss Matched_days = Treatment_day
#For naive mice 

naive_control_1 <- BW_FI_19 %>%
  ##LOOK HERE! change date directly below to "2026-05-06" ####
  filter((ID %in% c(3742, 3743, 3744, 3745, 3746)), DATE>"2026-05-06") 
         
naive_control_2 <- BW_FI_19 %>%
  filter((ID %in% c(3747,3748,3749,3750,3751,3752,3753)), DATE>"2026-04-29") 

naive_control_3 <- bind_rows(naive_control_1, naive_control_2)

naive_control_4 <- naive_control_3 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(matched_day = DATE - first(DATE),
         matched_day = as.numeric(matched_day)) %>%
  ungroup()

#In GLP1 tracker make matched_days = Treatment_day
#Then join that df with naive_control_4 by ID and DATE
#If ID is 3742-3753 then DRUG=naive otherwise DRUG=DRUG

All_mice_1 <- GLP1_tracker %>%
  left_join(naive_control_4 %>% 
      select(ID, DATE, matched_day),
    by = c("ID", "DATE"))  %>% 
  mutate(matched_day = if_else(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741), Treatment_day, matched_day),
         DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                          ID %in% c(3730, 3731, 3738, 3741) ~ "Tirzepatide",
                          ID %in% c(3732, 3733, 3737, 3739) ~ "Survodutide",
                          ID %in% c(3742, 3743, 3744, 3745, 3746, 3747,3748,3749,3750,3751,3752,3753) ~ "Naive")) %>%
  select(-SEX, -AIM, -DIET_FORMULA, -STRAIN) %>%
  filter(matched_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW))) %>%
  mutate(BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR)) %>%
  mutate(BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal))


##Plot ΔBW (%) -- SVD, TZP, Vehicle, Naive ####
ggplot(All_mice_1, aes(x=matched_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  #format.plot+
  #facet_wrap(~DRUG)+
  geom_hline(yintercept=0)+
 geom_vline(xintercept=16, linetype="dashed", color="darkgray")+ #increased doses
 geom_vline(xintercept=20, linetype="dashed", color="darkgray")+ #started doing injections every other day 
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "ΔBW (%) during BW loss",
       color="Treatment", fill="Treatment")

  
##Daily FI (kcal) -- SVD, TZP, Vehicle, Naive ####
ggplot(All_mice_1, aes(x=matched_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.25) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  #format.plot+
  geom_vline(xintercept=16, linetype="dashed", color="darkgray")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", color="darkgray")+ #started doing injections every other day
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Food intake (kcal/day)",
       color="Treatment", fill="Treatment")

#---
#Objective 2: Body comp ####
#Need to save the echoMRI data files with the correct formatting and title
  #to the echoMRI folder in the GitHub data folder. Then I need to re-run
  #data_proc so that an updated echoMRI.csv is written.


#Excerpt from data_proc.R to create the echo data frame and write the csv
# load metadata

metadata <- read_csv("../data/META.csv") %>% 
  select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# echo MRI ----
echomri_csv_files <- tibble(
  filepath = list.files("../data/echoMRI", full.names = TRUE)) %>% 
  filter(grepl("*.xlsx", filepath)) 
echomri_csv_files

echomri_open_files <- echomri_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map(., function(X){
    readxl::read_xlsx(X$filepath) %>% 
      select(Label, Fat, Lean, Weight, TimeDateDura) %>% 
      rename(ID = Label) %>% 
      separate_wider_delim(TimeDateDura, delim = ";", names = c("Date", "A", "B")) %>% 
      select(-A, -B) %>% 
      separate_wider_delim(Date, delim = " ", names = c("hms", "month", "day", "year")) %>% 
      mutate(day = gsub(",", "", day),
             Date = paste(year, month, day, sep = "-"),
             Date = lubridate::ymd(Date),
             ID =  as.factor(ID)) %>% 
      select(-hms, -month, -day, -year)}) %>% 
  bind_rows() %>% 
  left_join(., metadata, by = "ID")
echomri_open_files

# compare adiposity index = fat / lean ----
echomri_data <- echomri_open_files %>% 
  mutate(adiposity_index = Fat / Lean) %>% 
  group_by(ID) %>% 
  mutate(n_measurement = as.numeric(as.factor(Date)))
echomri_data

write_csv(x = echomri_data, "../data/echomri.csv")

##The code below will pull data from echoMRI.csv
#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_data_19 <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) %>%
  
Body_compt_BWloss <- echoMRI_data_19 %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                   ID %in% c(3730, 3731, 3738, 3741) ~ "Tirzepatide",
                   ID %in% c(3732, 3733, 3737, 3739) ~ "Survodutide"),
         DRUG = as.factor(DRUG)) %>%
  mutate(Lean_pct = 100*((Lean -first(Lean))/first(Lean)),
         Fat_pct = 100*((Fat-first(Fat))/first(Fat)),
         BW_pct = 100*((Weight-first(Weight))/first(Weight)),
         delta_AI= adiposity_index - first(adiposity_index), 
         #n_measurement=as.factor(n_measurement),
         n_measurement = factor(n_measurement, levels = c("Pre-treatment", "After 21 days")))

#Graph with x axis as n_measurement and y axis as Lean
Lean_pct_plot <-ggplot(Body_compt_BWloss, aes(x=n_measurement, y=Lean_pct, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
  #facet_wrap(~DRUG) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot +
  labs(x="Measurement",
       y= "Food intake (g/day)",
       title= "Daily FI (g) during BW loss",
       color="Treatment", fill="Treatment")
Lean_pct_plot

Fat_pct_plot <-ggplot(Body_compt_BWloss, aes(x=n_measurement, y=Fat_pct, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5)
Fat_pct_plot

AI_plot <-ggplot(Body_compt_BWloss, aes(x=n_measurement, y=adiposity_index, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5)
AI_plot

#Boxplots
ggplot(Body_compt_BWloss, aes(x=n_measurement, y=Lean))+ #, group=DRUG, fill=DRUG, color=DRUG) +
geom_boxplot(aes(color=DRUG)) + geom_point(aes(color=DRUG))


