#TZP --> incretin agonist Experiment 2 --> identify protocol for TZP induced weight loss ####
#Objective: monitor change in BW & FI during BW loss

#started script: 6-2-26
#revised script: 



# libraries
install.packages("mmand")
install.packages("pacman")
install.packages("this.path")

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

# Update BW.csv & FI_LM.csv ####
#Create df ####
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
#I changed the META.csv file on my local computer to include COHORT 19...
#NOTE: META.csv --> cohort 19 will be wiped away the next time I pull from origin ####
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

# Make df for analysis
Exp2_tracker <- BW_FI_19 %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0)) %>%
  mutate( 
    # NOTE: dates for the two rounds of sable recording will be different
    STAGE = case_when(ID %in% c(3742, 3743, 3744, 3745, 3748, 3749, 3750, 3751) 
                      & DATE< "2026-06-03" ~ "BW gain",
                      ID %in% c(3742, 3743, 3744, 3745, 3748, 3749, 3750, 3751) 
                      & DATE> "2026-06-02" ~ "BW loss",
                      
                      ID %in% c(3746, 3747, 3752, 3753) 
                      & DATE< "2026-06-07" ~ "BW gain",
                      ID %in% c(3746, 3747, 3752, 3753) 
                      & DATE> "2026-06-06" ~ "BW loss"),
    
    DRUG = case_when(ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) ~ "Vehicle",
                     ID %in% c(3742, 3743, 3747, 3750, 3751, 3753) ~ "Tirzepatide"),
    # All days since treatment started (including days w/o injection)
    Treatment_day =case_when(COMMENTS=="INJECT_DAY_1_DOSE_ONE"~0, 
                             COMMENTS=="INJECT_DAY_2_DOSE_ONE"~1,
                             COMMENTS=="INJECT_DAY_3_DOSE_ONE"~2,
                             COMMENTS=="INJECT_DAY_4_DOSE_ONE"~3,
                             COMMENTS=="INJECT_DAY_5_DOSE_ONE"~4,
                             COMMENTS=="INJECT_DAY_6_DOSE_ONE"~5,
                             COMMENTS=="INJECT_DAY_7_DOSE_ONE"~6,
                             COMMENTS=="INJECT_DAY_1_DOSE_TWO"~7,
                             COMMENTS=="INJECT_DAY_2_DOSE_TWO"~8,
                             COMMENTS=="INJECT_DAY_3_DOSE_TWO"~9,
                             COMMENTS=="INJECT_DAY_4_DOSE_TWO"~10,
                             COMMENTS=="INJECT_DAY_5_DOSE_TWO"~11,
                             COMMENTS=="INJECT_DAY_6_DOSE_TWO"~12,
                             COMMENTS=="INJECT_DAY_7_DOSE_TWO"~13,
                             COMMENTS=="INJECT_DAY_1_DOSE_THREE"~14,
                             COMMENTS=="INJECT_DAY_2_DOSE_THREE"~15,
                             COMMENTS=="INJECT_DAY_3_DOSE_THREE"~16,
                             COMMENTS=="INJECT_DAY_4_DOSE_THREE"~17,
                             COMMENTS=="INJECT_DAY_5_DOSE_THREE"~18,
                             COMMENTS=="INJECT_DAY_6_DOSE_THREE"~19,
                             COMMENTS=="INJECT_DAY_7_DOSE_THREE"~20),
    Dose_step = case_when(
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-06-02" & DATE < "2026-6-10" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE > "2026-6-09" & DATE < "2026-06-17" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      ID %in% c(3742, 3743, 3750, 3751) & STAGE == "BW loss" & 
        DATE >"2026-6-16" & DATE < "2026-6-23" & DRUG == "Tirzepatide" ~ "30nmol/kg",
      
      ID %in% c(3747, 3753) & STAGE == "BW loss" & 
        DATE > "2026-06-06" & DATE < "2026-6-14" & DRUG=="Tirzepatide" ~ "10nmol/kg",
      ID %in% c(3747, 3753) & STAGE == "BW loss" & 
        DATE > "2026-6-13" & DATE < "2026-06-21" & DRUG == "Tirzepatide" ~ "20nmol/kg",
      ID %in% c(3747, 3753) & STAGE == "BW loss" & 
        DATE >"2026-6-20" & DATE < "2026-6-28" & DRUG == "Tirzepatide" ~ "30nmol/kg",
      
      ID %in% c(3744, 3745, 3746, 3748, 3749, 3752) 
      & STAGE == "BW loss" & DRUG == "Vehicle" ~ "Vehicle")) %>%
  mutate(DRUG = as.factor(DRUG))


#Calculate % change in BW and cumulative intake ####
Exp2_BWloss_plot <- Exp2_tracker %>%
  filter(ID %in% c(3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753)) %>% #only IDs undergoing injections
  filter(Treatment_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)), #Calculate percent change in BW since BW loss started
         BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR),
         BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal)) 

# BW - GRAPHS ####

# Body weight graphs by drug group
## BW (g) ####
ggplot(Exp2_BWloss_plot, aes(x = Treatment_day, y = BW, fill = DRUG, color=DRUG)) +
  #geom_line(stat = "summary", 
  # fun = "mean") +
  geom_line(aes(y = BW, group = ID)) +
  geom_point() +
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(size=18, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  #facet_wrap(~ID)+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  labs(x="Day",
       y= "Body weight (g)",
       title= "Body weight during weight loss",
       color="Treatment", fill="Treatment") 

#---
##USE - BW (g) by DRUG during BW loss ####
Exp2_plot_BW_BWloss<-ggplot(Exp2_BWloss_plot, aes(x=Treatment_day, y=BW, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  #stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1), color=DRUG), vjust = -1, size = 4, 
  #fill = "white", fontface="bold", label.padding = unit(0.15, "lines"), show.legend=FALSE) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  format.plot+
  #facet_wrap(~ID)+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Body weight (g)",
       title= "Body weight during weight loss",
       color="Treatment", fill="Treatment")
Exp2_plot_BW_BWloss

#Export plot to Dose response folder
ggsave(Exp2_plot_BW_BWloss,
       filename="Exp2_TZP_BW_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

##USE - BW (% change) by DRUG during BW loss ####
Exp2_plot_pct_BW_BWloss<-ggplot(Exp2_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1), color=DRUG), vjust = -1, size = 4, 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines"), show.legend=FALSE) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  format.plot+
  geom_hline(yintercept=0)+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "ΔBody weight (%) ",
       title= "Percent ΔBW during BW loss",
       color="Treatment", fill="Treatment")
Exp2_plot_pct_BW_BWloss

#Export plot to Dose response folder
ggsave(Exp2_plot_pct_BW_BWloss,
       filename="Exp2_Incretin_BW_pct_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

#-----------------#
#-----------------#
# FI - GRAPHS ####
##Daily FI (g) ####
ggplot(Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  #stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1), color=DRUG), vjust = -1, size = 4, 
  #fill = "white", fontface="bold", label.padding = unit(0.15, "lines"), show.legend=FALSE) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
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
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (g/day)",
       title= "Daily food intake during BW loss",
       color="Treatment", fill="Treatment")

##USE - Daily FI (kcal) ####
Exp2_FI_daily_kcal_plot<-ggplot(Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  #stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1), color=DRUG), vjust = -1, size = 4, 
  #fill = "white", fontface="bold", label.padding = unit(0.15, "lines"), show.legend=FALSE) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary",fun.data = mean_se, aes(width=0.08), width=0.35) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
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
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Daily food intake during BW loss",
       color="Treatment", fill="Treatment")
Exp2_FI_daily_kcal_plot

#Export plot to Dose response folder
ggsave(Exp2_FI_daily_kcal_plot,
       filename="Exp2_FI_daily_kcal_plot.png", 
       width = 14, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

##USE - Cumulative FI (kcal) - by group ####
Exp2_FI_cum_kcal_plot<-ggplot(Exp2_BWloss_plot, aes(x = Treatment_day, y = BWloss_cum_INTAKE_kcal, color = DRUG, fill = DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=3) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID)) +
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
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Cumulative food intake (kcal)",
       title= "Cumulative food intake during BW loss",
       color="Treatment",
       fill="Treatment")
Exp2_FI_cum_kcal_plot

#Export plot to Dose response folder
ggsave(Exp2_FI_cum_kcal_plot,
       filename="Exp2_FI_cum_kcal_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

#---
#Manually check for correlation between drop in FI and plateau in BW loss ####

##CREATE df ####
Exp2_manual_check <- Exp2_BWloss_plot %>%
  select(DATE, ID, corrected_intake_gr, BW_pct_change, BW,  DRUG, Treatment_day) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE)

##GRAPH manual check ####
ggplot(Exp2_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (g)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=14),
        axis.text.x = element_text(),
        axis.title.y =element_text(face="bold", size= 14),
        axis.title.x =element_text(face="bold", size= 14),
        legend.title = element_text(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Exp2: Trajectory of FI & BW loss (NZO, n=12)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

###Graph just TZP ####
Exp2_TZP_manual_check <- Exp2_manual_check %>%
  filter(DRUG=="Tirzepatide")

ggplot(Exp2_TZP_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (grams)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(),
        axis.title.x = element_text(face="bold", size= 14),
        axis.title.y = element_text(face="bold", size= 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Exp 2. Tirzepatide: FI & BW loss (NZO, n=6)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

###Graph just Vehicle####
Exp2_vehicle_manual_check <- Exp2_manual_check %>%
  filter(DRUG=="Vehicle")

ggplot(Exp2_vehicle_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color= "darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=5, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (grams)"), size=5, color="darkgreen") +
  geom_hline(yintercept = 0) +
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(),
        axis.title.x = element_text(face="bold", size= 14),
        axis.title.y = element_text(face="bold", size= 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Exp 2. Vehicle: FI & BW loss (NZO, n=6)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

#------#
# Single drug graphs ####

## Tirzepatide ####
TZP_Exp2_BWloss_plot <- Exp2_BWloss_plot %>%
  filter(DRUG=="Tirzepatide")

### BW TZP ####
Exp2_TZP_plot_pct_BW_BWloss<-ggplot(TZP_Exp2_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        legend.text = element_blank(),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  geom_hline(yintercept=0) + 
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "Exp2. Tirzepatide: ΔBW (%)",
       color="Treatment", fill="Treatment")
Exp2_TZP_plot_pct_BW_BWloss

### FI TZP ####
####TZP FI --> grams ####
Exp2_TZP_FI_daily_gr_plot<-ggplot(TZP_Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (grams/day)",
       title= "Exp 2. Tirzepatide: Daily food intake (g)",
       color="Treatment", fill="Treatment")
Exp2_TZP_FI_daily_gr_plot

#### TZP FI --> kcal ####
Exp2_TZP_FI_daily_kcal_plot<-ggplot(TZP_Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Exp 2. Tirzepatide: Daily food intake (kcal)",
       color="Treatment", fill="Treatment")
Exp2_TZP_FI_daily_kcal_plot

#---------------------#
## Vehicle ####
Veh_Exp2_BWloss_plot <- Exp2_BWloss_plot %>%
  filter(DRUG=="Vehicle")

### BW Vehicle ####
Exp2_Veh_plot_pct_BW_BWloss<-ggplot(Veh_Exp2_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_blank(),
        #legend.text = element_blank(),
        # legend.title=element_text(size=15, face="bold"),
        # legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  geom_hline(yintercept=0) + 
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "Exp2. vehicle: ΔBW (%)",
       color="Treatment", fill="Treatment")
Exp2_Veh_plot_pct_BW_BWloss

###FI Vehicle ####
#### Vehicle FI --> grams ####
Exp2_Veh_FI_daily_gr_plot<-ggplot(Veh_Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (grams/day)",
       title= "Exp 2. Vehicle: Daily food intake (g)",
       color="Treatment", fill="Treatment")
Exp2_Veh_FI_daily_gr_plot

#### Vehicle FI --> kcal ####
Exp2_Veh_FI_daily_kcal_plot<-ggplot(Veh_Exp2_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  stat_summary(fun = mean,geom = "label",aes(label = round(after_stat(y), 1)), vjust = -3.3, size = 4, color = "black", 
               fill = "white", fontface="bold", label.padding = unit(0.15, "lines")) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values = custom_colors_GLP) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size=17, hjust = 0.5, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.text.x = element_text(size= 15),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", size= 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank()) +
  format.plot+
  #geom_vline(xintercept=7, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 10->20
  #geom_vline(xintercept=14, linetype="dotted", linewidth=1, color="#62748E")+ #increased dose 20->30
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Exp 2. Vehicle: Daily food intake (kcal)",
       color="Treatment", fill="Treatment")
Exp2_Veh_FI_daily_kcal_plot
