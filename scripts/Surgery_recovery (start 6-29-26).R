# Track BW and FI after surgery

#Start:
#Rev:

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
#custom_colors_GLP <- c("Tirzepatide" = "#1e6deb", "Vehicle" = "#403d3c")

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

#Taken from postsurgery_care_DREADD_surgery_uncertainty.R
post_surgery_tracker <- BW_FI_19 %>% 
  arrange(DATE) %>% 
  filter(ID %in% c(3731, 3732, 3733, 3735, 3737, 3738, 3739, 3740, 3741)) %>%
  group_by(ID) %>% 
  mutate(
    surgery_date = min(
      DATE[COMMENTS %in% c(
        "SURGERY")],na.rm = TRUE)) %>% 
  filter(DATE >= surgery_date) %>% 
  mutate(day_rel_surg = as.integer(as.Date(DATE) - as.Date(first(DATE)))) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0)) %>%
  mutate(
    STAGE = case_when(day_rel_surg ==0 ~ "Surgery",
                      day_rel_surg %in% c(1,2,3,4,5,6,7,8,9,10) ~ "Surgery recovery",
                      
                      #DATE=="2026-06-29" ~ "Day 1 sable acclimation" ~
                        
                      DATE=="2026-07-04" ~ "First mock",
                      DATE=="2026-07-07" ~ "Second mock",
                      DATE=="2026-07-10" ~ "Third mock"))

#Calculate % change in BW and cumulative intake ####
post_surgery_tracker_19 <- post_surgery_tracker %>%
  filter(day_rel_surg >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)), #Calculate percent change in BW since BW loss started
         BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR),
         BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal)) 

#BW graph ####

# Body weight change --> facet wrap graphs
## BW (g) ####
ggplot(post_surgery_tracker_19, aes(x = day_rel_surg, y = BW_pct_change)) +
  #geom_line(stat = "summary", 
  # fun = "mean") +
  geom_line(aes(y = BW_pct_change, group = ID)) +
  geom_point() +
 # scale_color_manual(values = custom_colors_GLP) +
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
  geom_vline(xintercept=3, linetype="dotted", linewidth=1, color="#62748E")+ #first day without melox
  geom_vline(xintercept=10, linetype="dotted", linewidth=1, color="#62748E")+ #end of 10 day recovery
  geom_hline(yintercept=0)+
  facet_wrap(~ID)+
  labs(x="Day post surgery",
       y= "ΔBody weight (%)",
       title= "ΔBody weight post surgery (%)") 

# Body weight graphs --> all IDs
## BW (g) ####
ggplot(post_surgery_tracker_19, aes(x = day_rel_surg, y = BW_pct_change)) +
  #geom_line(stat = "summary", 
  # fun = "mean") +
  geom_line(aes(y = BW_pct_change, group = ID)) +
  geom_point() +
  # scale_color_manual(values = custom_colors_GLP) +
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
  geom_hline(yintercept=0)+
  #facet_wrap(~ID)+
  geom_vline(xintercept=3, linetype="dotted", linewidth=1, color="#62748E")+ #first day without melox
  geom_vline(xintercept=10, linetype="dotted", linewidth=1, color="#62748E")+ #end of 10 day recovery
  labs(x="Day post surgery",
       y= "ΔBody weight (%)",
       title= "ΔBody weight post surgery (%)") 

#Manually check for correlation between drop in FI and plateau in BW loss ####

##CREATE df ####
post_surgery_tracker_19_manual_check <- post_surgery_tracker_19 %>%
  select(DATE, ID, corrected_intake_gr, BW_pct_change, BW,  day_rel_surg) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE)

##GRAPH manual check ####
ggplot(post_surgery_tracker_19_manual_check, aes(x = day_rel_surg)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (g)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=14),
        axis.text.x = element_text(),
        axis.title.y =element_text(face="bold", size= 14),
        axis.title.x =element_text(face="bold", size= 14),
        legend.title = element_text(),
        legend.text = element_text(size = 14)) +
  labs(x="Day post surgery",
       title= "Post surgery FI & BW (NZO, n=9)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)
