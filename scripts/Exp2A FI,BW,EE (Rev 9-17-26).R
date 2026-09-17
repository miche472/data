#Experiment 2A = Tirzepatide only ####
#Tirzepatide only. Administered for 1 week at 10nmol/kg and for a 2nd week at 20nmol/kg

#Data for Vijay

#Script started: 9-15-26
#Script revised: 9-16-26

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
library(grid)

#Set working directory
setwd("/Users/laurenmichels/Documents/GitHub/data/data")

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
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")

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
custom_colors1_GLP_Exp2 <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")
custom_colors2_GLP_Exp2 <- c("Tirzepatide" = "#1E90FF", "Vehicle" = "#8B8989")
custom_colors3_GLP_Exp2 <- c("Tirzepatide" = "#104E8B", "Vehicle" = "#403d3c")

# Update BW.csv & FI_LM.csv ####
#Organize data into df ####
# bodyweight and food intake

# libs 
pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand
)

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



## Exp 2A BW & FI ####
Exp2A_FI_BW <- BW_FI_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3748, 3749, 3750, 3751, 3752, 3753)) %>% #exclude 3747
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0)) %>%
  mutate( 
    STAGE = case_when(ID %in% c(3742, 3743, 3744, 3745, 3748, 3749, 3750, 3751) 
                          & DATE< "2026-06-03" ~ "BW gain",
                      ID %in% c(3742, 3743, 3744, 3745, 3748, 3749, 3750, 3751) 
                          & DATE> "2026-06-02" & DATE<"2026-06-17" ~ "BW loss",
                      ID %in% c(3742, 3743, 3750, 3751) 
                          & DATE>"2026-06-16" ~ "BW regain",
                      ID %in% c(3744, 3745, 3748, 3749)
                          & DATE>"2026-06-16" ~ "Exp_2B",
 
                      ID %in% c(3746, 3752, 3753) 
                          & DATE< "2026-06-08" ~ "BW gain",
                      ID %in% c(3746, 3752, 3753) 
                          & DATE>= "2026-06-08" & DATE<"2026-06-22" ~ "BW loss",
                      ID == 3753 
                          & DATE> "2026-06-21" ~ "BW regain",
                      ID %in% c(3746, 3752)
                          & DATE > "2026-06-21" ~ "Exp_2B"),

    DRUG = case_when(ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
                     ID %in% c(3742, 3743, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
  
    # All days since treatment started
   mutate(Treatment_day =case_when(COMMENTS=="INJECT_DAY_1_DOSE_ONE"~0,
                             COMMENTS=="SABLE_DAY_4_INJECT_DAY_1_DOSE_ONE"~0, #For 3752,3753
                             
                             COMMENTS=="INJECT_DAY_2_DOSE_ONE"~1,
                             COMMENTS=="INJECT_DAY_3_DOSE_ONE"~2,
                             COMMENTS=="INJECT_DAY_4_DOSE_ONE"~3,
                             COMMENTS=="INJECT_DAY_5_DOSE_ONE"~4,
                             COMMENTS=="INJECT_DAY_6_DOSE_ONE"~5,
                             COMMENTS=="INJECT_DAY_7_DOSE_ONE"~6,
                             
                             COMMENTS=="INJECT_DAY_1_DOSE_TWO"~7,
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_1_DOSE_TWO"~7, #For 3746
                             
                             COMMENTS=="INJECT_DAY_2_DOSE_TWO"~8,
                             
                             
                             COMMENTS=="INJECT_DAY_3_DOSE_TWO"~9,
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_3_DOSE_TWO"~9, #For 3742,3743,3750,3751,3753
                             
                             
                             COMMENTS=="INJECT_DAY_4_DOSE_TWO"~10,
                             COMMENTS=="SABLE_DAY_2_INJECT_DAY_4_DOSE_TWO"~10, #For 3742,3743,3750,3751,3753
                             
                             
                             COMMENTS=="INJECT_DAY_5_DOSE_TWO"~11,
                             COMMENTS=="SABLE_DAY_3_INJECT_DAY_5_DOSE_TWO"~11, #For 3742,3743,3750,3751,3753
                             
                             COMMENTS=="INJECT_DAY_6_DOSE_TWO"~12,
                             COMMENTS=="SABLE_DAY_4_INJECT_DAY_6_DOSE_TWO"~12, #for 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_1_INJECT_DAY_6_DOSE_TWO"~12, #For 3746,3752
                             
                             COMMENTS=="INJECT_DAY_7_DOSE_TWO"~13,
                             COMMENTS=="SABLE_DAY_5_INJECT_DAY_7_DOSE_TWO" ~13, #for 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_2_INJECT_DAY_7_DOSE_TWO"~13, #For 3746,3752
                             
                             COMMENTS=="SABLE_DAY_6_REGAIN_DAY_1"~14,
                             COMMENTS=="SABLE_DAY_1" & ID %in% c(3744, 3745, 3748) & DATE>"2026-06-02" & DATE < "2026-06-24"~14, #For 3742,3743,3750,3751,3753
                             COMMENTS=="SABLE_DAY_3" & ID %in% c(3746,3752) & DATE>"2026-06-02" & DATE < "2026-06-24"~14,
                             ID==3749 & DATE == "2026-06-17"~14,
                             
                             COMMENTS=="SABLE_DAY_2" & ID %in% c(3744, 3745, 3748) &  DATE>"2026-06-02" & DATE <"2026-06-24"~15,
                             ID == 3749 & DATE == "2026-06-18"~15,
                             COMMENTS=="SABLE_DAY_4" & ID %in% c(3746,3752) & DATE>"2026-06-02" & DATE <"2026-06-24"~15,
                             
                             COMMENTS=="SABLE_DAY_3" & ID %in% c(3744, 3745, 3748) &  DATE>"2026-06-02" & DATE <"2026-06-24"~16,
                             ID == 3749 & DATE == "2026-06-19"~16,
                             COMMENTS=="REGAIN_DAY_3" & ID %in% c(3742,3743,3750,3751,3753)~16, 
                             
                             COMMENTS=="SABLE_DAY_4" & ID %in% c(3744, 3745, 3748) & DATE>"2026-06-02" & DATE<"2026-06-24"~17,
                             ID==3749 & DATE == "2026-06-20"~17,
                             COMMENTS=="REGAIN_DAY_4" & ID %in% c(3742,3743,3750,3751,3753)~17,
                             
                             #ID==3749 & DATE == "2026-06-21"~18,
                             COMMENTS=="REGAIN_DAY_5" & ID %in% c(3742,3743,3750,3751,3753)~18,
                             
                             #ID == 3749 & DATE == "2026-06-22"~19,
                             COMMENTS=="REGAIN_DAY_6" & ID %in% c(3742,3743,3750,3751,3753) ~19,
                             COMMENTS=="REGAIN_DAY_7" & ID %in% c(3742,3743,3750,3751,3753) ~20,
                             COMMENTS=="REGAIN_DAY_8" & ID %in% c(3742,3743,3750,3751,3753) ~21,
                             COMMENTS=="REGAIN_DAY_9" & ID %in% c(3742,3743,3750,3751,3753) ~22,
                             COMMENTS=="REGAIN_DAY_10" & ID %in% c(3742,3743,3750,3751,3753) ~23,
                             COMMENTS=="REGAIN_DAY_11" & ID %in% c(3742,3743,3750,3751,3753) ~24,
                             COMMENTS=="REGAIN_DAY_12" & ID %in% c(3742,3743,3750,3751,3753) ~25,
                             COMMENTS=="REGAIN_DAY_13" & ID %in% c(3742,3743,3750,3751,3753) ~26), 
    Dose_step = case_when(
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-06-02" & DATE < "2026-6-10" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-6-09" & DATE < "2026-06-17" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      
      ID ==3753 & STAGE == "BW loss" & 
        DATE > "2026-06-05" & DATE < "2026-6-14" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID == 3753 & STAGE == "BW loss" & 
        DATE > "2026-6-13" & DATE < "2026-06-21" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      
      ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) 
      & STAGE == "BW loss" ~ "Vehicle")) %>%
  mutate(DRUG = as.factor(DRUG))
  #drop_na(STAGE)

##Df for BW and FI plots ####
Exp2A_BW_loss_plots <- Exp2A_FI_BW %>%
  filter(Treatment_day %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)), #Calculate percent change in BW since BW loss started
         BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR),
         BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal)) 



#------------------------------------------------------------#.
#------------------------------------------------------------#.
#Food intake (FI) ####
# x axis --> treatment day and y axis -->cumulative intake. One line for vehicle and one line for tirzepatide

## Cumulative FI (grams) ####
Exp2A_plot_cumFI_G <-ggplot(Exp2A_BW_loss_plots, aes(x=Treatment_day, y=BWloss_cum_INTAKE_GR, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size= 13, angle=0, vjust=0.5, hjust=0.7)) +
  #geom_vline(xintercept=6, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  labs(x="Day",
       y= "Cumulative intake (g)",
       title= "Cumulative food intake (g)",
       color="Treatment", fill="Treatment")
Exp2A_plot_cumFI_G

#Export plot to Incretin_dose folder
ggsave(Exp2A_plot_cumFI_G,
       filename="Exp2A_TZP_FI_cumG_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")

## Cumulative FI (kcal) ####
Exp2A_plot_cumFI_kcal <-ggplot(Exp2A_BW_loss_plots, aes(x=Treatment_day, y=BWloss_cum_INTAKE_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size= 13, angle=0, vjust=0.5, hjust=0.7)) +
  #geom_vline(xintercept=6, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  labs(x="Day",
       y= "Cumulative intake (kcal)",
       title= "Cumulative food intake (kcal)",
       color="Treatment", fill="Treatment")
Exp2A_plot_cumFI_kcal

#Export plot to Incretin_dose folder
ggsave(Exp2A_plot_cumFI_kcal,
       filename="Exp2A_TZP_FI_cum_kcal_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")



#------------------------------------------------------------#.
#------------------------------------------------------------#.
#Body weight (BW) ####
## BW (grams) ####
Exp2A_plot_BW_G <-ggplot(Exp2A_BW_loss_plots, aes(x=Treatment_day, y=BW, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size= 13, angle=0, vjust=0.5, hjust=0.7)) +
  #geom_vline(xintercept=6, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  labs(x="Day",
       y= "Body weight (g)",
       title= "Body weight (g)",
       color="Treatment", fill="Treatment")
Exp2A_plot_BW_G

#Export plot to Incretin_dose folder
ggsave(Exp2A_plot_BW_G,
       filename="Exp2A_TZP_BW_G_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")

## BW (% change) ####
Exp2A_plot_BW_pct_change <-ggplot(Exp2A_BW_loss_plots, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size= 13, angle=0, vjust=0.5, hjust=0.7)) +
  #geom_vline(xintercept=6, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  geom_hline(yintercept=0, linetype="solid", linewidth=1, color="black")+ 
  labs(x="Day",
       y= "Body weight (% change)",
       title= "Change in body weight (%)",
       color="Treatment", fill="Treatment")
Exp2A_plot_BW_pct_change

#Export plot to Incretin_dose folder
ggsave(Exp2A_plot_BW_pct_change,
       filename="Exp2A_TZP_BW_pct_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")


#------------------------------------------------------------#.
#------------------------------------------------------------#.
#Body composition ####

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


## Exp. 2A Body comp
echoMRI_exp2a <- echoMRI_data_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3744,3745,3746,3748,3749,3752) ~ "Vehicle",
                          ID %in% c(3742,3743,3750,3751,3753) ~ "Tirzepatide"),
         DRUG = as.factor(DRUG),
         Measurement = case_when( n_measurement == "1" ~ "Baseline", 
                                  n_measurement == "2" ~ "Pre acclimation", 
                                   n_measurement == "3" ~ "Peak obesity", 
                                   n_measurement == "4" ~ "During BW loss",
                                   n_measurement == "5" & 
                                   ID %in% c(3742,3743,3749,3750,3751,3753)~ "BW regain")) %>%
  drop_na(Measurement) %>%
  mutate(Measurement = factor(Measurement, 
                           levels = c("Baseline", 
                                      "Pre acclimation", 
                                      "Peak obesity", 
                                      "During BW loss", # 9 days into BW loss for all mice except 3753 (4 days) 
                                      "BW regain"))) 

# Calculate change in BW, lean, fat, & adiposity index
echmoMRI_exp2a_delta <- echoMRI_exp2a %>%
  filter(Measurement %in% c("Peak obesity", "During BW loss", "BW regain"))%>%
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

#df with only obesity and during BW loss
plot_echmoMRI_exp2a_delta <- echmoMRI_exp2a_delta %>%
  filter(Measurement %in% c("Peak obesity", "During BW loss"))
  #filter(!(ID==3753)) #only 4 days into BW loss

#---
## Graph ΔBW (%) ####
Echo_BW_pct_change_Exp2A <-ggplot(plot_echmoMRI_exp2a_delta, aes(x=Measurement, y=BW_pct_change, 
                                                  group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  labs(y= "Δ Body weight (%)", 
       title= "Change in body weight (%)", 
       color="Treatment", fill="Treatment")
Echo_BW_pct_change_Exp2A

#Export plot to Incretin_dose folder on desktop computer 
ggsave(BW_pct_change_exp2a,
       filename="Exp2A_echo_BW_pct_change.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")

## Graph ΔLean (%)  ####
Echo_Lean_pct_change_Exp2a <-ggplot(plot_echmoMRI_exp2a_delta, aes(x=Measurement, y=Lean_pct_change, 
                                                  group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  labs(y= "ΔLean mass (%)", 
       title= "Change in lean mass (%)", 
       color="Treatment", fill="Treatment")
Echo_Lean_pct_change_Exp2a

#Export plot to Incretin_dose folder on desktop computer 
ggsave(Echo_Lean_pct_change_Exp2a,
       filename="Exp2A_echi_lean_pct_change.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")


#------------------------------------------------------------#.
#------------------------------------------------------------#.
#Energy expenditure ####

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

# Code to calculate Energy expenditure data (df called Minute_sum_EE_hr_alldays) ####
#Read in Sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#Functions
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
  filter(!(ID==3747)) %>% #Remove since I had to stop TZP injections early
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
  mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
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
  filter(!(ID==3747)) %>% #Remove since I had to stop TZP injections early
  filter(date<"2026-06-29") %>% #Remove recordings from BW loss for 3744, 45, 46, 48, 52...tried 10->15nmol/kg TZP
    mutate(Time = as_hms(format(DateTime, "%H:%M:%S")),
    lights = if_else(Time >= as_hms("06:30:00") & Time < as_hms("18:30:00"),"on", "off"),
    SABLE= case_when(
      sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                       "SABLE_DAY_4") ~ "Peak obesity",
      ID %in% c(3742,3743,3750,3751,3753) & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7",
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
    ID %in% c(3742, 3743, 3750, 3751, 3753) ~ "Tirzepatide")) %>%
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

#Filter for just complete days 1 and 2
EE_day_1_2 <- new_filter_loc_TEE2 %>%
  filter(complete_days %in% c(1,2))

Minute_sum_EE_hr <- EE_day_1_2 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  # Movement each hr (during how many minutes did mouse move)
  group_by(SABLE, ID, complete_days) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr, lights, complete_days) %>% # group by lights to retain photoperiod in df for later use
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

#Calculate daily sum (add all hours together)
Daily_EE <- Minute_sum_EE_hr %>%
  ungroup() %>%
  group_by(ID, SABLE, complete_days) %>%
  summarise(TEE_kcal_day = sum(TEE_kcal),
            NEAT_kcal_day = sum(NEAT_kcal),
            RMR_kcal_day = sum(RMR_kcal),
            diff= abs(TEE_kcal_day - sum(NEAT_kcal_day + RMR_kcal_day))) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide"))

#Average between complete days for each ID at each Sable time point
Avg_complete_days <- Daily_EE %>%
  ungroup() %>%
  group_by(ID, SABLE) %>%
  summarise(Avg_TEE_kcal_day = mean(TEE_kcal_day),
            Avg_NEAT_kcal_day = mean(NEAT_kcal_day),
            Avg_RMR_kcal_day = mean(RMR_kcal_day)) %>%
  mutate(DRUG = case_when(
    ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
    ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide"))


##EE graphs ####
 # x - axis sable time point, y-axis = TEE, colors = Treatment group

Exp2A_plot_EE_kcal_day <- ggplot(
  Avg_complete_days,
  aes(x = SABLE, y = Avg_RMR_kcal_day, color = DRUG)) +
  geom_col(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7,
    linewidth = 1.5, fill = "white") +
  geom_jitter(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 1) +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size = 13,angle = 0,vjust = 0.5,hjust = 0.7),
        axis.title.x = element_blank()) +
  labs(x = "Stage", 
       y = "TEE (kcal/day)", 
       title = "Total energy expenditure", color = "Treatment")
Exp2A_plot_EE_kcal_day

#Export plot to Incretin_dose folder
ggsave(Exp2A_plot_EE_kcal_day,
       filename="Exp2A_TZP_EE_kcal_day.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures/Incretin_dose")

#-------------------------------------------------------------------------------------#.
#-------------------------------------------------------------------------------------#.
#-------------------------------------------------------------------------------------#.
#-------------------------------------------------------------------------------------#.
#  DIFFERENT STUDY: kcal restriction induced BW loss. ####
#---

format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.minor = element_blank(), # remove background grid lines only
  panel.grid.major = element_blank(),
  axis.line = element_line(color = "black")) # keep axis lines
custom_colors_345 <- c("Control" = "#FAAC41","Weight cycled" = "#3498DB")
custom_colors_345_2 <- c("Control" = "#E67E22","Weight cycled" = "#1d5e8a")

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



echoMRI_data_345 <- echoMRI_data %>%
  filter(COHORT %in% c(3,4,5)) %>%
  mutate(ID= as.factor(ID)) %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) 

echoMRI_data_BWloss_345 <- echoMRI_data_345 %>%
  filter(!(ID %in% c(3712, 3715))) %>% # died during study
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #tech issue with sable cages
  ungroup() %>%
  group_by(ID) %>%
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3711, 3713, 3716, 3719, 3726) ~ "Control",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3727, 3728, 3729) ~ "Weight cycled"),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3711, 3713, 3714, 3720, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3719, 3721, 3722, 3726, 3729) ~ "RTIOXA_47"),
    GROUP = as.factor(GROUP)) %>%
mutate(SABLE = case_when(
           n_measurement == 1 ~ "Baseline",
           Date %in% as.Date(c("2025-01-27", "2025-01-29", "2025-02-07")) ~ "Peak obesity",
           ID %in% c(3711, 3727) & Date == as.Date("2025-02-20") ~ "Peak obesity",
           Date == as.Date("2025-03-28") ~ "BW loss",
           Date == as.Date("2025-05-27") ~ "BW maintenance",
           Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                               "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
           TRUE ~ NA_character_)) %>%
mutate(
    SABLE = factor(SABLE, 
                           levels = c("Baseline", 
                                      "Peak obesity", 
                                      "BW loss", 
                                      "BW maintenance",
                                      "BW regain"))) 
  
  
  
## Calculations 
#change in lean, fat, total body mass, adiposity index
echmoMRI_345_delta <- echoMRI_data_BWloss_345 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(BW = Fat + Lean) %>%
  mutate(Lean_pct_BW = 100*(Lean/BW),
         Fat_pct_BW = 100*(Fat/BW))

# Plot BW, lean, and fat
plot_echmoMRI_345_delta <- echmoMRI_345_delta %>%
  filter(SABLE %in% c("Baseline", "Peak obesity", "BW loss", "BW regain"))
  
#---
#---
##Graph BW (g) ####
BW_plot_G_345 <-ggplot(plot_echmoMRI_345_delta, aes(x=SABLE, y=BW, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=GROUP), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=GROUP), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
  #geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  theme(axis.text.x = element_text(size = 13,angle = 20,vjust = 0.5,hjust = 0.7),
        axis.title.x = element_blank(),
        title = element_blank()) +
  labs(y= "Body weight (g)",
       #title= "Body weight (g)",
       color="Treatment", fill="Treatment")
BW_plot_G_345

##Graph lean mass (g) ####
Lean_plot_G_345 <-ggplot(plot_echmoMRI_345_delta, aes(x=SABLE, y=Lean, 
                                                       group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", fun = "mean", aes(color=GROUP), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=GROUP), linewidth=1) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  #geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  labs(y= "Lean mass (g)",
       title= "Lean mass (g)",
       color="Treatment", fill="Treatment")
Lean_plot_G_345

##Graph fat mass (g) ####
Fat_plot_G_345 <-ggplot(plot_echmoMRI_345_delta, aes(x=SABLE, y=Fat, 
                                                       group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", fun = "mean", aes(color=GROUP), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=GROUP), linewidth=1.5) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.1) +
  #geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM2 +
  labs(y= "Fat mass (g)",
       title= "Fat mass (g)",
       color="Treatment", fill="Treatment")
Fat_plot_G_345


## Percent change in BW/Lean/Fat ####
echmoMRI_345_pct_change <- echoMRI_data_BWloss_345 %>%
  filter(SABLE %in% c("Peak obesity", "BW loss", "BW regain")) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(BW = Fat + Lean) %>%
  mutate(delta_lean = Lean - first (Lean), # negative number means lean mass was lost
         delta_fat = Fat - first(Fat),
         delta_BW = BW - first(BW),
         delta_AI= adiposity_index - first(adiposity_index),
         delta_lean_vs_delta_BW = 100*(delta_lean/delta_BW), # I think this is the percent of BW loss which is lean mass
         BW_pct_change = 100*((BW - first(BW)) / first(BW)),
         Lean_pct_change = 100*((Lean - first(Lean)) / first(Lean)),
         Fat_pct_change = 100*((Fat - first(Fat)) / first(Fat)),
         Lean_pct_BW = 100*(Lean/BW),
         Fat_pct_BW = 100*(Fat/BW))

### Graph ΔBW (%) ####
BW_plot_pct_change_345 <-ggplot(echmoMRI_345_pct_change, aes(x=SABLE, y=BW_pct_change, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=GROUP), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=GROUP), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(size = 13,angle = 20,vjust = 0.5,hjust = 0.7),
        axis.title.x = element_blank(),
        title = element_blank()) +
  labs(y= "Δ Body weight (%)",
       title= "Change in Body weight (%)",
       color="Treatment", fill="Treatment")
BW_plot_pct_change_345

#Export plot to figures folder
ggsave(BW_plot_pct_change_345,
       filename="Cal_restrict_BW_pct_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

### Graph ΔLean (%) ####
plot_Lean_pct_change_345 <-ggplot(echmoMRI_345_pct_change, aes(x=SABLE, y=Lean_pct_change, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=GROUP), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=GROUP), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(size = 13,angle = 20,vjust = 0.5,hjust = 0.7),
        axis.title.x = element_blank()) +
  labs(y= "Δ Lean (%)",
       title= "Change in lean mass (%)",
       color="Treatment", fill="Treatment")
plot_Lean_pct_change_345

#Export plot to figures folder
ggsave(plot_Lean_pct_change_345,
       filename="Cal_restrict_Lean_pct_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")

### Graph ΔFat (%) ####
plot_Fat_pct_change_345 <-ggplot(echmoMRI_345_pct_change, aes(x=SABLE, y=Fat_pct_change, group=GROUP, fill=GROUP, color=GROUP)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=GROUP), size=4, position = position_dodge(width = 0)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=GROUP), linewidth=1, position = position_dodge(width = 0)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_color_manual(values = custom_colors_345) +
  theme_bw(base_size = 14) +
  format.plot_LM +
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(size = 13,angle = 20,vjust = 0.5,hjust = 0.7),
        axis.title.x = element_blank()) +
  labs(y= "Δ Fat (%)",
       title= "Change in fat mass (%)",
       color="Treatment", fill="Treatment")
plot_Fat_pct_change_345

#Export plot to figures folder
ggsave(plot_Fat_pct_change_345,
       filename="Cal_restrict_Fat_pct_plot.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/figures")
