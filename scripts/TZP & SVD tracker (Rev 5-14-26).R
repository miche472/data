#Track BW and FI of NZO mice (n=11) during SVD and TZP treatment for BW loss

#Started:4-20-26
#Revised: 5-14-26

#This script started as Incretin agonist tracker (Rev 5-1-26). On 5-4-26 I made
#major changes to allow for multiple doses of TZP and SVD, so I changed the name

#Objective: monitor change in BW and FI during BW loss for incretin receptor agonist study

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
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb","Survodutide" = "#c93618", "Vehicle" = "#403d3c")
custom_colors_GLPB <- c("Tirzepatide" = "lightblue","Survodutide" = "lightgreen", "Vehicle" = "darkgray")


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

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

# Update BW.csv & FI_LM.csv ####
#Create df ####
# bodyweight and food intake
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
#I changed the META.csv file on my local computer to include COHORT 19...this 
#will be wiped away the next time I pull from origin
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

##Read in BW and FI ####
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

##Create GLP1_tracker df ####
#Modified this to allow for multiple doses of SVD and TZP and for days without injections within the treatment period
GLP1_tracker <- BW_FI_19 %>%
    filter_out(ID == "3734") %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0)) %>%
mutate(
    STAGE = case_when(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741) 
                        & DATE< "2026-04-17" ~ "BW gain",
                      ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741) 
                      & DATE> "2026-04-16" ~ "BW loss"),
    DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                        ID %in% c(3730, 3731, 3738, 3741) ~ "Tirzepatide",
                        ID %in% c(3732, 3733, 3737, 3739) ~ "Survodutide"),
    Treatment_day =case_when(COMMENTS=="INJECT_DAY_ONE"~0, # All days since treatment started (including days w/o injection)
                          COMMENTS=="INJECT_DAY_TWO"~1,
                          COMMENTS=="INJECT_DAY_THREE"~2,
                          COMMENTS=="INJECT_DAY_FOUR"~3,
                          COMMENTS=="INJECT_DAY_FIVE"~4,
                          COMMENTS=="INJECT_DAY_SIX"~5,
                          COMMENTS=="INJECT_DAY_SEVEN"~6,
                          COMMENTS=="INJECT_DAY_EIGHT"~7,
                          COMMENTS=="INJECT_DAY_NINE"~8,
                          COMMENTS=="INJECT_DAY_TEN"~9,
                          COMMENTS=="INJECT_DAY_ELEVEN"~10,
                          COMMENTS=="INJECT_DAY_TWELVE"~11,
                          COMMENTS=="INJECT_DAY_THIRTEEN"~12,
                          COMMENTS=="INJECT_DAY_FOURTEEN"~13,
                          COMMENTS=="INJECT_DAY_FIFTEEN"~14,
                          COMMENTS=="INJECT_DAY_SIXTEEN"~15,          #Sat 5/2
                          COMMENTS=="INJECT_DAY_ONE_DOSE_TWO"~16,     #Sun 5/3
                          COMMENTS=="INJECT_DAY_TWO_DOSE_TWO"~17,     #Mon 5/4
                          COMMENTS=="INJECT_DAY_THREE_DOSE_TWO"~18,   #Tues 5/5
                          COMMENTS=="INJECT_DAY_FOUR_DOSE_TWO"~19,    #Wed 5/6
                          COMMENTS=="NO_INJECT_DAY_ONE"~20,           #Thurs 5/7
                          COMMENTS=="INJECT_DAY_FIVE_DOSE_TWO"~21,    #Fri 5/8
                          COMMENTS=="NO_INJECT_DAY_TWO"~22,           #Sat 5/9
                          COMMENTS=="INJECT_DAY_SIX_DOSE_TWO"~23,     #Sun 5/10
                          COMMENTS=="INJECT_DAY_SEVEN_DOSE_TWO"~24,   #Mon 5/11
                          COMMENTS=="INJECT_DAY_EIGHT_DOSE_TWO"~25,   #Tues 5/12
                          COMMENTS=="INJECT_DAY_NINE_DOSE_TWO"~26,    #Wed 5/13
                          COMMENTS=="INJECT_DAY_TEN_DOSE_TWO"~27,     #Thurs 5/14
                          COMMENTS=="INJECT_DAY_ELEVEN_DOSE_TWO"~28,
                          COMMENTS=="INJECT_DAY_TWELVE_DOSE_TWO"~29),    
    
    Inject_day= case_when(COMMENTS=="INJECT_DAY_ONE"~0, # Only days during treatment when an injection was given
                          COMMENTS=="INJECT_DAY_TWO"~1,
                          COMMENTS=="INJECT_DAY_THREE"~2,
                          COMMENTS=="INJECT_DAY_FOUR"~3,
                          COMMENTS=="INJECT_DAY_FIVE"~4,
                          COMMENTS=="INJECT_DAY_SIX"~5,
                          COMMENTS=="INJECT_DAY_SEVEN"~6,
                          COMMENTS=="INJECT_DAY_EIGHT"~7,
                          COMMENTS=="INJECT_DAY_NINE"~8,
                          COMMENTS=="INJECT_DAY_TEN"~9,
                          COMMENTS=="INJECT_DAY_ELEVEN"~10,
                          COMMENTS=="INJECT_DAY_TWELVE"~11,
                          COMMENTS=="INJECT_DAY_THIRTEEN"~12,
                          COMMENTS=="INJECT_DAY_FOURTEEN"~13,
                          COMMENTS=="INJECT_DAY_FIFTEEN"~14,
                          COMMENTS=="INJECT_DAY_SIXTEEN"~15,
                          COMMENTS=="INJECT_DAY_ONE_DOSE_TWO"~16,
                          COMMENTS=="INJECT_DAY_TWO_DOSE_TWO"~17,
                          COMMENTS=="INJECT_DAY_THREE_DOSE_TWO"~18,
                          COMMENTS=="INJECT_DAY_FOUR_DOSE_TWO"~19,
                          COMMENTS=="INJECT_DAY_FIVE_DOSE_TWO"~21,
                          COMMENTS=="INJECT_DAY_SIX_DOSE_TWO"~23,
                          COMMENTS=="INJECT_DAY_SEVEN_DOSE_TWO"~24,
                          COMMENTS=="INJECT_DAY_EIGHT_DOSE_TWO"~25,
                          COMMENTS=="INJECT_DAY_NINE_DOSE_TWO"~26,
                          COMMENTS=="INJECT_DAY_TEN_DOSE_TWO"~27,     
                          COMMENTS=="INJECT_DAY_ELEVEN_DOSE_TWO"~28,
                          COMMENTS=="INJECT_DAY_TWELVE_DOSE_TWO"~29),
    Dose_step = case_when(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)
                           & STAGE == "BW loss" & DATE< "2026-5-03"& DRUG=="Survodutide" ~ "30nmol/kg",
                           
                           ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741) 
                           & STAGE == "BW loss" & DATE> "2026-5-02" & DRUG == "Survodutide" ~ "40nmol/kg",
                           
                           ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)
                           & STAGE == "BW loss" & DATE< "2026-5-03"& DRUG=="Tirzepatide" ~ "10nmol/kg",
                           
                           ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741) 
                           & STAGE == "BW loss" & DATE> "2026-5-02" & DRUG == "Tirzepatide" ~ "20nmol/kg",
                          
                          ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741) 
                          & STAGE == "BW loss" & DRUG == "Vehicle" ~ "Vehicle")) %>%
    mutate(DRUG = as.factor(DRUG)) 


#Prepare to graph BW during weight loss
#Calculate percent change in BW since BW loss started
GLP1_BWloss_plot <- GLP1_tracker %>%
  filter(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>% #only IDs undergoing injections
  filter(Treatment_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)))
  #filter_out(ID == "3732") #Try this since 3732 has a much lower BW and looks like an outlier


#BW during BW loss####
#Graph BW for each mouse
ggplot(GLP1_BWloss_plot,
       aes(x = Treatment_day, y = BW,
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Body weight (grams)",
       title= "BW (grams) during BW loss")

# Body weight graphs by drug group
##Graph BW (grams) ####
ggplot(GLP1_BWloss_plot, aes(x = Treatment_day, y = BW, fill = DRUG, color=DRUG)) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "BW (g)",
       title= "BW (g) during BW loss",
       color="Treatment", fill="Treatment") 

##USE - Graph BW (g) by group during BW loss ####
plot_BW_BWloss<-ggplot(GLP1_BWloss_plot, aes(x=Treatment_day, y=BW, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  #facet_wrap(~ID)+
  labs(x="Day",
       y= "Body weight (g)",
       title= "Body weight during weight loss",
       color="Treatment", fill="Treatment")
plot_BW_BWloss

#Export plot to Dose response folder
ggsave(plot_BW_BWloss,
       filename="Incretin_BW_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

##USE - Graph BW (% change) by group during BW loss ####
plot_pct_BW_BWloss<-ggplot(GLP1_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
geom_errorbar(stat = "summary", 
              fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "ΔBody weight (%) ",
       title= "Percent ΔBW during BW loss",
       color="Treatment", fill="Treatment")
plot_pct_BW_BWloss

#Export plot to Dose response folder
ggsave(plot_pct_BW_BWloss,
       filename="Incretin_BW_pct_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

#-----------------#
#-----------------#
#FI during BW loss####
GLP1_BWloss_plot2 <- GLP1_tracker %>%
  #filter_out(ID == "3732") %>% #Try this since 3732 has a much lower BW and looks like an outlier
  filter(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>% #only IDs undergoing injections
  filter(Treatment_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW))) %>%
  mutate(BWloss_cum_INTAKE_GR = cum_INTAKE_GR - first(cum_INTAKE_GR)) %>%
  mutate(BWloss_cum_INTAKE_kcal = cum_INTAKE_kcal - first(cum_INTAKE_kcal))


# FI graphs by drug group
##Daily FI (g) ####
ggplot(GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (g/day)",
       title= "Daily food intake during BW loss",
       color="Treatment", fill="Treatment")

##USE - Daily FI (kcal) ####
Incretin_FI_daily_kcal_plot<-ggplot(GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  #geom_point(shape=1) + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Daily food intake during BW loss",
       color="Treatment", fill="Treatment")
Incretin_FI_daily_kcal_plot

#Export plot to Dose response folder
ggsave(Incretin_FI_daily_kcal_plot,
       filename="Incretin_FI_daily_kcal_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

##Cumulative FI (kcal) ####

#By ID
#Cumulative FI during BW loss -- trajectory for each ID 
ggplot(GLP1_BWloss_plot2, aes(x = Treatment_day, y = BWloss_cum_INTAKE_kcal, color = DRUG, fill = DRUG)) +
  #geom_line(linewidth = 1.2) +
  geom_line(aes(y = BWloss_cum_INTAKE_kcal, group = ID)) +
  geom_point(aes(y = BWloss_cum_INTAKE_kcal, group = ID)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_color_manual(values = custom_colors_GLP) +
  #facet_wrap(~ID) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal)",
       title= "Cumulative food intake during BW loss")

#---
##USE - Cumulative FI (kcal) - by group ####
Incretin_FI_cum_kcal_plot<-ggplot(GLP1_BWloss_plot2, aes(x = Treatment_day, y = BWloss_cum_INTAKE_kcal, color = DRUG, fill = DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal)",
       title= "Cumulative food intake during BW loss",
       color="Treatment",
       fill="Treatment")
Incretin_FI_cum_kcal_plot

#Export plot to Dose response folder
ggsave(Incretin_FI_cum_kcal_plot,
       filename="Incretin_FI_cum_kcal_plot.png", 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300,
       path = "/Users/laurenmichels/Desktop/Dose response/Incretin dose/Figures")

#---
##Cumulative FI (g) - by group ####
ggplot(GLP1_BWloss_plot2, aes(x = Treatment_day, y = BWloss_cum_INTAKE_GR, color = DRUG, fill = DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
  geom_point(shape=1) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Cumulative FI (g)",
       title= "Cumulative FI (g) during BW loss")

#---
#Graph FI during BW loss using INTAKE_GR rather than correct_intake_gr --> should be the same for all 
#injection days except the first injection day --> this holds true
ggplot(GLP1_BWloss_plot2, aes(x = Treatment_day, y = corrected_intake_gr, color = DRUG, fill = DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.5) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.25) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Treatment day",
       y= "Intake_GR",
       title= "INTAKE_GR during BW loss (n=11)")

#Statistics ####
##Summary statistics ####
#
GLP1_summary_stats <- GLP1_BWloss_plot %>%
  mutate(Treatment_day = factor(Treatment_day, levels = c("0","1","2","3","4","5","6", 
                    "7","8","9","10","11","12","13","14","15","16","17","18",
                    "19","20","21","22","23","24","25","26","27","28","29","30"))) %>%
  ungroup() %>%
  group_by(DRUG, Treatment_day) %>%
  summarise(n=n(),
            Avg_BW= mean(BW),
            sd_BW=sd(BW), 
            se_BW= sd_BW / sqrt(n),
            
            Avg_BW_pct = mean(BW_pct_change),
            sd_BW_pct = sd(BW_pct_change),
            se_BW_pct = sd_BW_pct / sqrt(n),
            
            Avg_cum_INTAKE_gr=mean(cum_INTAKE_GR),
            Avg_cum_INTAKE_kcal = mean(cum_INTAKE_kcal),
            Avg_daily_FI_gr = mean(corrected_intake_gr),
            
            Avg_daily_FI_kcal = mean(corrected_intake_kcal),
            sd_daily_FI_kcal = sd(corrected_intake_kcal), 
            se_daily_FI_kcal = sd_daily_FI_kcal / sqrt(n))

##LMM: ΔBW (%) ####
#Make df that has Treatment_day as a factor not a numerical variable
GLP1_for_model <- GLP1_BWloss_plot %>%
  mutate(Treatment_day = factor(Treatment_day, levels = c("0", "1","2","3","4","5","6", 
                    "7","8","9","10","11","12","13","14","15","16","17","18",
                    "19","20","21","22","23","24","25","26","27","28","29","30")))
#Make linear mixed model for ΔBW (%)
model_BW_pct <- lmer(BW_pct_change ~ Treatment_day*DRUG + (1 | ID), data = GLP1_for_model)
summary(model_BW_pct)

# Calculate estimated marginal means 
emm_BW_pct <- emmeans(model_BW_pct, ~ Treatment_day*DRUG, cov.reduce = mean)
emm_BW_pct_df <- as.data.frame(emm_BW_pct)

# Pairwise contrasts within each DRUG
contrasts_by_group_BW_pct <- contrast(emm_BW_pct, method = "pairwise", by = "DRUG")
contrasts_by_group_BW_pct_df <- as.data.frame(contrasts_by_group_BW_pct) %>%
  filter(p.value<0.05)

# Pairwise contrasts within each Injection day (time point)
contrasts_by_Treatment_day_BW_pct <- contrast(emm_BW_pct, method = "pairwise", by = "Treatment_day")
contrasts_by_Treatment_day_BW_pct_df <- as.data.frame(contrasts_by_Treatment_day_BW_pct) %>%
  filter(p.value<0.05)

##LMM: Daily FI (kcal) ####
#Make df that has Treatment_day as a factor not a numerical variable
GLP1_for_model <- GLP1_BWloss_plot %>%
  mutate(Treatment_day = factor(Treatment_day, levels = c("0", "1","2","3","4","5","6", 
                                                    "7","8","9","10","11","12","13","14","15","16","17","18",
                                                    "19","20","21","22","23","24","25","26","27","28","29","30")))
#Make linear mixed model for Daily FI (corrected_intake_kcal)
model_FI_kcal <- lmer(corrected_intake_kcal ~ Treatment_day*DRUG + (1 | ID), data = GLP1_for_model)
summary(model_FI_kcal)

# Calculate estimated marginal means 
emm_FI_kcal <- emmeans(model_FI_kcal, ~ Treatment_day*DRUG, cov.reduce = mean)
emm_FI_kcal_df <- as.data.frame(emm_FI_kcal)

# Pairwise contrasts within each DRUG
contrasts_by_group_FI_kcal <- contrast(emm_FI_kcal, method = "pairwise", by = "DRUG")
contrasts_by_group_FI_kcal_df <- as.data.frame(contrasts_by_group_FI_kcal) %>%
  filter(p.value<0.05)

# Pairwise contrasts within each Injection day (time point)
contrasts_by_Treatment_day_FI_kcal <- contrast(emm_FI_kcal, method = "pairwise", by = "Treatment_day")
contrasts_by_Treatment_day_FI_kcal_df <- as.data.frame(contrasts_by_Treatment_day_FI_kcal) %>%
  filter(p.value<0.05)

#Plateau date would be when two consecutive days were not significantly different 

#---------#
#---------#

#Manually check for correlation between drop in FI and plateau in BW loss ####
#select only ID, DATE, Treatment_day, DRUG, BW_pct_change, BW, INTAKE_GR

##CREATE data frame ####
manual_check <- GLP1_BWloss_plot2 %>%
  select(DATE, ID, corrected_intake_gr, BW_pct_change, BW,  DRUG, Treatment_day) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE)

##GRAPH manual check ####
ggplot(manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
    geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                   color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (g)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (g)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=14),
        axis.text.x = element_text(),
        axis.title.y =element_text(face="bold", size= 14),
        axis.title.x =element_text(face="bold", size= 14),
        legend.title = element_text(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Trajectory of FI and BW loss (NZO, n=11)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

####Graph just SVD ####
SVD_manual_check <- manual_check %>%
  filter(DRUG=="Survodutide")

ggplot(SVD_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (grams)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(),
        axis.title.x =element_text(face="bold", size= 14),
        axis.title.y =element_text(face="bold", size= 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Survodutide: FI and BW loss (NZO, n=11)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

###Graph just TZP ####
TZP_manual_check <- manual_check %>%
  filter(DRUG=="Tirzepatide")

ggplot(TZP_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=3, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (grams)"), size=3, color="darkgreen") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(),
        axis.title.x = element_text(face="bold", size= 14),
        axis.title.y = element_text(face="bold", size= 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Tirzepatide: FI and BW loss (NZO, n=11)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)


###Graph just Control####
Vehicle_manual_check <- manual_check %>%
  filter(DRUG=="Vehicle")

ggplot(Vehicle_manual_check, aes(x = Treatment_day)) +
  geom_line(aes(y = BW_pct_change, color = "BW (% change)"), color= "darkviolet") +
  geom_point(aes(y = BW_pct_change, color = "BW (% change)"), color="darkviolet") +
  geom_label(aes(y = BW_pct_change, label = round(BW_pct_change,1), nudge_y = -4, 
                 color="BW (% change)"), size=5, color="darkviolet") +
  geom_line(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_point(aes(y = corrected_intake_gr, color = "Daily FI (grams)"), color="darkgreen") +
  geom_label(aes(y = corrected_intake_gr, label = round(corrected_intake_gr,1), nudge_y = 5, 
                 color="Daily FI (grams)"), size=5, color="darkgreen") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(),
        axis.title.x = element_text(face="bold", size= 14),
        axis.title.y = element_text(face="bold", size= 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  labs(x="Day",
       title= "Vehicle: FI and BW loss (NZO, n=11)",
       y="ΔBW (%) & Daily FI (g)") +
  facet_wrap(~ID)

#------#
# Single drug graphs ####

## Survodutide ####
SVD_GLP1_BWloss_plot <- GLP1_BWloss_plot %>%
  filter(DRUG=="Survodutide")

### BW SVD ####
SVD_plot_pct_BW_BWloss<-ggplot(SVD_GLP1_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  #format.plot+
  geom_hline(yintercept=0) + 
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "Survodutide: ΔBW (%)",
       color="Treatment", fill="Treatment")
SVD_plot_pct_BW_BWloss

###FI SVD ####
#### FI SVD --> gr ####
SVD_Incretin_FI_daily_gr_plot<-ggplot(SVD_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  #format.plot+
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (grams/day)",
       title= "Survodutide: Food intake (grams/day)",
       color="Treatment", fill="Treatment")
SVD_Incretin_FI_daily_gr_plot

#### FI SVD --> kcal ####
SVD_Incretin_FI_daily_kcal_plot<-ggplot(SVD_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Survodutide: Food intake (kcal/day)",
       color="Treatment", fill="Treatment")
SVD_Incretin_FI_daily_kcal_plot

#---
## Tirzepatide ####
TZP_GLP1_BWloss_plot <- GLP1_BWloss_plot %>%
  filter(DRUG=="Tirzepatide")

### BW TZP ####
TZP_plot_pct_BW_BWloss<-ggplot(TZP_GLP1_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "Tirzepatide: ΔBW (%)",
       color="Treatment", fill="Treatment")
TZP_plot_pct_BW_BWloss

###FI TZP ####
####TZP FI --> grams ####
TZP_Incretin_FI_daily_gr_plot<-ggplot(TZP_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  #format.plot+
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (grams/day)",
       title= "Tirzepatide: Food intake (grams/day)",
       color="Treatment", fill="Treatment")
TZP_Incretin_FI_daily_gr_plot

#### TZP FI --> kcal ####
TZP_Incretin_FI_daily_kcal_plot<-ggplot(TZP_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Tirzepatide: Food intake (kcal/day)",
       color="Treatment", fill="Treatment")
TZP_Incretin_FI_daily_kcal_plot


#---
## Vehicle ####
Veh_GLP1_BWloss_plot <- GLP1_BWloss_plot %>%
  filter(DRUG=="Vehicle")

### BW Vehicle ####
Veh_plot_pct_BW_BWloss<-ggplot(Veh_GLP1_BWloss_plot, aes(x=Treatment_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point(aes(group=ID), shape=1) + #this adds individual values from each ID to the graph
  geom_line(aes(group=ID), alpha=0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se, aes(width=0.08), width=0.45, linewidth=1, color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "ΔBW (%) ",
       title= "Vehicle: ΔBW (%)",
       color="Treatment", fill="Treatment")
Veh_plot_pct_BW_BWloss

###FI Vehicle ####
#### Vehicle FI --> grams ####
Veh_Incretin_FI_daily_gr_plot<-ggplot(Veh_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (grams/day)",
       title= "Vehicle: Food intake (grams/day)",
       color="Treatment", fill="Treatment")
Veh_Incretin_FI_daily_gr_plot

#### Vehicle FI --> kcal ####
Veh_Incretin_FI_daily_kcal_plot<-ggplot(Veh_GLP1_BWloss_plot, aes(x=Treatment_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", fun = "mean", aes(color=DRUG), size=4) +
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
  geom_vline(xintercept=16, linetype="dotted", linewidth=1, color="#62748E")+ #increased doses
  geom_vline(xintercept=20, linetype="dashed", linewidth=1, color="#62748E")+ #started doing injections every other day 
  geom_vline(xintercept=23, linetype="dashed", linewidth= 1, color="#62748E")+ #final day of every other day injections 
  labs(x="Day",
       y= "Food intake (kcal/day)",
       title= "Vehicle: Food intake (kcal/day)",
       color="Treatment", fill="Treatment")
Veh_Incretin_FI_daily_kcal_plot




#-----------------------------#
#-----------------------------#
# Double check that I did the sorting correctly ####
#ANOVAs for AI sorted 
#One way ANOVA for AI sorted: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkBW <- aov(Weight ~ DRUG, data = echoMRI_dataC)
summary(ANOVA_AI_sorted_checkBW)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkBW, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkBW))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.4669 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(Weight ~ DRUG, xlab='DRUG', ylab='Weight', data=echoMRI_dataC)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(Weight ~ DRUG, data=echoMRI_dataC)
#p value is 0.8058 --> samples have equal variance


#-----------------------------------------------#
#-----------------------------------------------#
#Which method for calculating corrected FI is better? ####
#Used this to determine if the original method or my new method for calculated corrected FI was better
df <- cohort_open_files %>%
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
  mutate(delta_measurement = DATE - lag(DATE)) %>% #delete after comparing to delta_alt
  #drop_na(delta_alt) %>%
  mutate(corrected_intake_gr_alt = INTAKE_GR / as.numeric(delta_alt)) %>%
  mutate(corrected_intake_gr_original = INTAKE_GR / as.numeric(delta_measurement)) %>%
  mutate(diff=corrected_intake_gr_original - corrected_intake_gr_alt)
left_join(., food_desc, by = "DIET_FORMULA") %>%
  mutate(corrected_intake_kcal = corrected_intake_gr * KCAL_G) %>%
  left_join(., metadata, by = "ID")

#########-----#


#Original FI code  
# output food-intake file
#In this version drop_na(delta_measurement) only removes the first observation
FI <- cohort_open_files %>% 
  select(ID, DIET_FORMULA, INTAKE_GR, DATE,COMMENTS) %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(
    delta_measurement = DATE - lag(DATE)
  ) %>% 
  drop_na(delta_measurement) %>% 
  mutate(
    corrected_intake_gr = INTAKE_GR / as.numeric(delta_measurement) #DAILY INTAKE
  ) %>% 
  left_join(., food_desc, by = "DIET_FORMULA") %>% 
  mutate(
    corrected_intake_kcal = corrected_intake_gr * KCAL_G
  ) %>% 
  left_join(., metadata, by = "ID")

###########-----#
