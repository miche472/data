# GLP-1 experiment 2A and 2B

#EchoMRI data after BW loss for NZO mice
#Exp 2A: dose escalation 10nmol TZP -> 20nmol TZP, each for one week in NZO mice
#Exp 2B: dose escalation 10nmol TZP -> 15nmol TZP, each for one week in NZO mice
  #Exp 2B was an exploration and did not have a control group

#Started with obese female NZO mice
#BW loss occurred over 14 days

#Lean mass data requested in 7-1-26 lab meeting by Laurie and Cathy

#Exp 2A: two groups -> TZP and vehicle
  #1st echo was when they reached adulthood
  #2nd echo was pre sable acclimation
  #3rd echo was immediately before start of BW loss (Peak obesity)
  #4th echo was day before BW loss sable recording started for TZP mice. Midway through 
      #saline injections for vehicle mice ocurred on the same date as 4th echo for TZP mice though.
  #5th echo was BW regain for TZP mice. Vehicle mice had been moved on to exp 2B
  #so I don't have a BW regain for vehicle mice in exp 2A

#Exp 2B: one group -> TZP mice. Use mice that had vehicle during Exp 2A
    # Baseline could be considered the 4th echo
    #End BW loss could be considered the 5th echo (July 6th or 8th)
    #End regain in last week of July

#LEFT OFF --> "Measurement" isn't working as an x variable in graphs...n_measurement does work

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
custom_colors_GLP_Exp2 <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#read in metadata
metadata <- read_csv("../data/META.csv") %>% 
select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# echo MRI
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

# compare adiposity index = fat / lean
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
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) 


# Exp. 2A ####
echoMRI_exp2a <- echoMRI_data_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3744,3745,3746,3748,3749,3752) ~ "Vehicle",
                          ID %in% c(3742,3743,3747,3750,3751,3753) ~ "Tirzepatide"),
         DRUG = as.factor(DRUG),
         Measurement = case_when(n_measurement == "1" ~ "Baseline I", 
                                   n_measurement == "2" ~ "Baseline II", 
                                   n_measurement == "3" ~ "Peak obesity", 
                                   n_measurement == "4" ~ "BW loss",
                                   n_measurement == "5" & 
                                   ID %in% c(3742,3743,3747,3750,3751,3753)~ "BW regain")) %>%
  drop_na(Measurement) %>%
  mutate(Measurement = factor(n_measurement, 
                           levels = c("Baseline I", 
                                      "Baseline II", 
                                      "Peak obesity", 
                                      "BW loss",      
                                      "BW regain"))) 



# Calculate change in BW, lean, fat, & adiposity index ####
echmoMRI_exp2A_delta <- echoMRI_exp2a %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(BW = Fat + Lean) %>%
  mutate(delta_lean = Lean - first (Lean), # negative number means lean mass was lost
         delta_fat = Fat - first(Fat),
         delta_BW = BW - first(BW),
         
         Lean_BW = delta_lean/delta_BW, 
         delta_lean_vs_delta_BW = 100*(Lean_BW), # % of BW loss which is lean mass
         
         BW_pct_change = 100*((BW - first(BW)) / first(BW)),         # % change in BW
         Lean_pct_change = 100*((Lean - first(Lean)) / first(Lean)), # % change in lean mass
         Fat_pct_change = 100*((Fat - first(Fat)) / first(Fat)),     # % change in fat mass

         Lean_pct_BW = 100*(Lean/BW), # % of BW comprised of lean mass
         Fat_pct_BW = 100*(Fat/BW))   # % of BW comprised of fat mass


  
####Lean mass (g) ####
#Graph with x axis as n_measurement and y axis as Lean
Lean_plot_exp2A <-ggplot(echmoMRI_exp2A_delta, aes(x=n_measurement, y=Lean, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP_Exp2) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15, angle=25, vjust=0.5),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot +
  labs(#x="Day",
       y= "Lean mass (g)",
       title= "Lean mass (g)",
       color="Treatment", fill="Treatment")
Lean_plot_exp2A


####Fat mass (g) ####
#Graph with x axis as n_measurement and y axis as Fat
Fat_plot_exp2a <-ggplot(echmoMRI_exp2A_delta, aes(x=Measurement, y=Fat, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP_Exp2) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15, angle=25, vjust=0.5),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot +
  labs(y= "Fat mass (g)",
       title= "Fat mass (g)",
       color="Treatment", fill="Treatment")
Fat_plot_exp2a