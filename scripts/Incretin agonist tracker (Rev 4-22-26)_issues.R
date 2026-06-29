#Track BW and FI of mice once they start injections

#Started:4-20-26
#Revised: 4-21-26

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

# Write BW.csv & FI.csv --> up to date versions ####
# libraries
install.packages("mmand")
install.packages("pacman")
install.packages("this.path")

library(mmand)
library(pacman)
library(this.path)

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

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
         day_rel = as.numeric(day_rel)) %>%
  replace_na(list(INTAKE_GR=0, delta_alt=0, corrected_intake_gr=0,
                  corrected_intake_kcal=0, KCAL_G=3.82))

GLP1_tracker <- BW_FI_19 %>%
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
                      & DATE> "2026-04-15" ~ "BW loss"),
    DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                        ID %in% c(3730, 3731, 3738, 3741) ~ "Survodutide",
                        ID %in% c(3732, 3733, 3737, 3739) ~ "Tirzepatide"),
    Inject_day= case_when(COMMENTS=="INJECT_DAY_ONE"~0, 
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
                          COMMENTS=="INJECT_DAY_TWELVE"~11), 
    Group_number= case_when(ID %in% c(3730, 3731, 3732, 3733,3735, 3736, 3737, 
                                         3738, 3739, 3740, 3741) ~"First",
                            ID %in% c(3747, 3748, 3749, 3750, 3751, 3752, 3753)~"Second",
                            ID %in% c(3742, 3743, 3744, 3745, 3746)~"Third")) %>%
    mutate(DRUG = as.factor(DRUG))


#Prepare to graph BW during weight loss
#Calculate percent change in BW since BW loss started
GLP1_BWloss_plot <- GLP1_tracker %>%
  filter(Inject_day >=0) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)))


#BW during weight loss####
#Graph BW for each mouse
ggplot(GLP1_BWloss_plot,
       aes(x = Inject_day, y = BW,
           color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
 # facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x="# of days with injections",
       y= "Body weight (grams)",
       title= "NZO BW during injections")

# Body weight graphs by drug group
##Graph BW (grams) ####
ggplot(GLP1_BWloss_plot, aes(x = Inject_day, y = BW, fill = DRUG, color=DRUG)) +
  geom_line(stat = "summary", 
           fun = "mean") +
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
  facet_wrap(~ID)+
  labs(x="Day",
       y= "BW (g)",
       title= "BW (g) during BW loss",
       color="Treatment", fill="Treatment") 

##Graph BW (% change) ####
ggplot(GLP1_BWloss_plot, aes(x=Inject_day, y=BW_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1.8) +
  geom_point(shape=1) + #this adds individual values from each ID to the graph
geom_errorbar(stat = "summary", 
              fun.data = mean_se, aes(width=0.08)) +
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
  geom_hline(yintercept=0, linetype="dashed", color="darkgray")+
  labs(x="Day",
       y= "Body weight (%) ",
       title= "Change in BW during BW loss (%)",
       color="Treatment", fill="Treatment")


#FI ####
# FI graphs by drug group
##Graph FI (g) ####
ggplot(GLP1_BWloss_plot, aes(x=Inject_day, y=corrected_intake_gr, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point() + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08)) +
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
  labs(x="Day",
       y= "Daily food intake (g)",
       title= "Daily FI (g) during BW loss",
       color="Treatment", fill="Treatment")


##Graph FI (kcal) ####
ggplot(GLP1_BWloss_plot, aes(x=Inject_day, y=corrected_intake_kcal, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=2) +
  geom_point() + #this adds individual values from each ID to the graph
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08)) +
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
  labs(x="Day",
       y= "Daily food intake (kcal)",
       title= "Daily FI (kcal) during BW loss",
       color="Treatment", fill="Treatment")

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

#########


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

###########
