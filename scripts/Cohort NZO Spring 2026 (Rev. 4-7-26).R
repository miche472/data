#Process BW and FI data for COHORT_19

#Started:3-31-26
#Last revised: 4-7-26 
#Last run on: 4-14-26

#Objectives ####
#Purpose of this script is to: 
  #1. View daily changes in BW/Fi for incretin agonist study
      #Involves: Manually enter data into COHORT_19.csv, located in GitHub, then 
      #run script to write updated BW.csv and FI.csv
      #IMPORTANT NOTE: when I pull from origin in GitHub, the modified version of 
      #META.csv and updated COHORT_19.csv will be overwritten
  #2. Assign mice to drug treatment groups based on adiposity index or BW
    #Assign and check whether the assignment successfully created treatment 
    #groups with the same mean BW and adiposity index
  #3. Health checks: quickly identify if a mouse is reducing FI or losing BW (prior to start of BW loss)
  #4. Track BW loss during injections and quickly spot if a mouse is approaching a 20% reduction in BW

# Write BW.csv & FI.csv --> up to date versions ####
# libraries
install.packages("mmand")
install.packages("pacman")
install.packages("this.path")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

# change the directory to source file location
#setwd(this.path::here())

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
FI <- cohort_open_files %>% 
  select(ID, DIET_FORMULA, INTAKE_GR, DATE,COMMENTS) %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(
    delta_measurement = DATE - lag(DATE)) %>% 
  drop_na(delta_measurement) %>% 
  mutate(
    corrected_intake_gr = INTAKE_GR / as.numeric(delta_measurement) #DAILY INTAKE
  ) %>% 
  left_join(., food_desc, by = "DIET_FORMULA") %>% 
  mutate(
    corrected_intake_kcal = corrected_intake_gr * KCAL_G) %>% 
  left_join(., metadata, by = "ID")

# output bodyweight file
BW <- cohort_open_files %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, BW, DATE,COMMENTS) %>% 
  drop_na(BW) %>% 
  left_join(., metadata, by = "ID")

write_csv(x = FI, "../data/FI.csv")
write_csv(x = BW, "../data/BW.csv")

#Read in BW and FI ####
#Read in BW and filter for cohort 19 (Spring 2026 NZO mice)
BW_COHORT19 <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>%
  filter(COHORT == 19)

#Read in FI and filter for cohort 19 (Spring 2026 NZO mice)
FI_COHORT19 <- read_csv("~/Documents/GitHub/data/data/FI.csv") %>%
  filter(COHORT == 19)

#Join FI and BW
BW_FI_19 <- BW_COHORT19 %>%
  left_join(
    FI_COHORT19 %>% 
      select(ID, INTAKE_GR, DATE, delta_measurement, corrected_intake_gr, corrected_intake_kcal, KCAL_G),
    by = c("ID", "DATE")) %>%
  mutate(ID = as.factor(ID)) %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(day_rel = DATE - first(DATE),
         day_rel = as.numeric(day_rel)) %>%
  replace_na(list(INTAKE_GR=0, delta_measurement=0, corrected_intake_gr=0,
                  corrected_intake_kcal=0, KCAL_G=3.82))
  

# Create BW & FI graphs ####

##BW ####
#Graph BW
ggplot(BW_FI_19,
       aes(x = day_rel, y = BW,
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
  labs(x="Days with ad libitum LFD",
       y= "Body weight (grams)",
       title= "NZO bodwy weight (Cohort 19)")

#---
##Cumulative FI ####

#Calculate cumulative FI
FI_19 <- BW_FI_19 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0))

#Graph cumulative FI
ggplot(FI_19, aes(x = day_rel, y = cum_INTAKE_kcal, color = ID, fill = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ID) +
  format.plot+
  theme_bw(base_size = 14) +
  theme( legend.position = "right",
         plot.title = element_text(hjust = 0.5, face="bold", size=15),
         axis.text.x = element_text(angle = 45, hjust = 1),
         legend.title = element_text(hjust=0.5, size = 13, face = "bold"),
         legend.text = element_text(size = 11)) +
  labs(x="Days with ad libitum LFD",
       y= "Cumulative energy intake (kcal)",
       title= "NZO Cumulative energy intake (kcal)")

#---
##Daily FI (aka corrected) ####

#Graph corrected intake --> get a sense of when this stabilizes
ggplot(FI_19, aes(x = day_rel, y = corrected_intake_kcal, color = ID, fill = ID)) +
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
  labs(x="Days with ad libitum LFD",
       y= "Energy intake (kcal/day)",
       title= "NZO Daily energy intake (kcal/day)")


#Assign to treatment groups ####
#Try using multiple methods. (either adiposity index or BW) ##
#Drug groups: TZP (10nmol/kg), SVD (30 nmol/kg), vehicle (1X PBS, 5 ml/kg)

#After assigning groups, interrogate group assignation methods 
    #1. If I use BW for assignation, do the resulting groups have different average 
        #adiposity indices?
    #2. If I use adiposity index for assignation, do the resulting groups have 
        #different avg BW?
        
#---
##Method 1: Use adiposity index to assign mice to groups ####
#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_dataA<- echoMRI_data %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3734, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  #filter(Date == "2026-04-16") %>%
  select(ID, adiposity_index, Weight, Date) %>%
  ungroup()

#Create pattern for group assignment
pattern1 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

echoMRI_dataB <- echoMRI_dataA %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern1, length.out = n())) %>%
  ungroup()

#Have R randomly assign TZP, SVD, and Vehicle to Groups 1-3
#[Assume that I use the adiposity index group assignment (i.e. method 1)]
# Create a random mapping between treatment groups and drugs
Group_to_drug_map <- setNames(
  sample(c("TZP", "SVD", "Vehicle")),
  c("Group1", "Group2", "Group3"))

# Apply the mapping to my df called echoMRI_dataB
echoMRI_dataC <- echoMRI_dataB %>%
  mutate(Drug = Group_to_drug_map[treatment])

#One way ANOVA: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW <- aov(Weight ~ treatment, data = echoMRI_dataC)
summary(ANOVA_AI_sorted_checkBW)
         
#One way ANOVA: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkAI <- aov(adiposity_index ~ treatment, data = echoMRI_dataC)
summary(ANOVA_AI_sorted_checkAI)

#---
##Method 2: Use BW to assign mice to groups ####
#Use echoMRI_data which was read in for method 1 --> 
          #for method 2 I will use body "Weight" included in this df
Weight_dataA <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3734, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  #filter(Date == "2026-04-16") %>%
  select(ID, adiposity_index, Weight, Date) %>%
  ungroup()

#Create pattern for group assignment
#Have R randomly assign TZP, SVD, and Vehicle to Groups 1-3
pattern1 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

Weight_dataB <- Weight_dataA %>%
  arrange(desc(Weight)) %>%
  mutate(treatment = rep(pattern1, length.out = n()))

#One way ANOVA: same mean BW for all three groups?
ANOVA_BW_sorted_checkBW <- aov(Weight ~ treatment, data = Weight_dataB)
summary(ANOVA_BW_sorted_checkBW)

#One way ANOVA: same mean adiposity index for all three groups?
ANOVA_BW_sorted_checkAI <- aov(adiposity_index ~ treatment, data = Weight_dataB)
summary(ANOVA_BW_sorted_checkAI)
#---






#----------#
#----------#
#Used this section to figure out how to write the code for assignation (didn't have 
#Echo data at this point)
#Based the code for assignation off of a script called 
            #"Treatment gorup assignment NZO (May 2025).R"

#Use most recent BW
BW_dataA <- FI_19 %>%
  filter(DATE == "2026-03-30") %>%
  select(ID, BW, DATE) %>%
  ungroup()

#Create pattern for group assignment
#Have R randomly assign TZP, SVD, and Vehicle to Groups 1-3
pattern2 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

BW_dataB <- BW_dataA %>%
  arrange(desc(BW)) %>%
  mutate(treatment = rep(pattern2, length.out = n()))

#Now, have R randomly assign TZP, SVD, and Vehicle to Groups 1-3
# Create a random mapping between treatment groups and drugs
Group_to_drug_map <- setNames(
  sample(c("TZP", "SVD", "Vehicle")),
  c("Group1", "Group2", "Group3"))

# Apply the mapping to my df
BW_dataC <- BW_dataB %>%
  mutate(Drug = Group_to_drug_map[treatment])

#One way ANOVA: same mean BW for all three groups?
ANOVA_BW_sorted_checkBW <- aov(BW ~ Drug, data = BW_dataC)
summary(ANOVA_BW_sorted_checkBW)

#One way ANOVA: same mean adiposity index for all three groups?
#Can't do this until I actually have the echoMRI data for the mice


#---#
#Written on 4/7/26

#BW Tracker: LFD -> percent change since start of LFD ####
df_LFD_19 <- BW_FI_19 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    LFD_date = if_else(any(COMMENTS == "FIRST_DAY_LFD"), # LFD start date
                       DATE[COMMENTS == "FIRST_DAY_LFD"][1], 
                       as.Date(NA)),
    day_rel_LFD = as.numeric(DATE - LFD_date), # Days relative to LFD
    LFD_day0 = BW[which.min(abs(day_rel_LFD))], #BW on first day of LFD
    BW_pct_change = 100*((BW - LFD_day0) / LFD_day0),
    BW_pct_change = as.numeric(BW_pct_change)) %>% #% change from start LFD
    ungroup()

##Health checks - LFD ####
#BW --> More than 10% decrease since start of LFD?
BW_Health_check <- df_LFD_19 %>%
  mutate(Dropped_10_percent = if_else(BW_pct_change < -10, "Yes", "No")) %>%
  filter(Reached_10 == "Yes")

#---#
#FI --> Identify if a mouse is eating less than it has been eating previously?
# in other words, is corrected_intake_gr less than the rolling average intake
FI_Health_check <- BW_FI_19 %>%
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, DATE, corrected_intake_gr) %>%
  drop_na(corrected_intake_gr) %>% 
  mutate(moving_avg = rollmean(corrected_intake_gr, k=3, fill=NA, align='right'))
  
#Took the rolling avg calculation method from "Food restriction (2-19-25) LM, LL (LM version).R"
#Confirm the rollmean using spaghetti plots (verify no errors in data entry)
  #Need multiple aes because you want multiple lines. 
  #Black line shows raw corrected intake values. Red line shows moving average
sp1 <- FI_Health_check %>% 
  mutate(date = lubridate::ymd(DATE)) %>% 
  ggplot(aes(date, corrected_intake_gr)) +
  geom_line(aes(group=ID)) +
  geom_line(aes(date, moving_avg, group=ID), color="red") +
  facet_wrap(~ID, scale="free_y")
sp1

FI_Health_check_2g <- FI_Health_check %>%
  ungroup() %>%
  group_by(ID) %>%
  filter(DATE > '2026-04-3') %>%
  filter(moving_avg > corrected_intake_gr) 

#---------------------------------------------#
#BW tracker: Injections -> percent change since start of BW loss treatment ####
df_inject_19 <- BW_FI_19 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    inject_date = if_else(any(COMMENTS == "FIRST_DAY_INJECT"), #Treatment start date
                     DATE[COMMENTS == "FIRST_DAY_INJECT"][1], 
                     as.Date(NA)),
    day_rel_inject = as.numeric(DATE - inject_date), # Days relative to start of injections
    Inject_day0 = BW[which.min(abs(day_rel_inject))], #BW on first day of injections
    BW_pct_change = 100*((BW - LFD_day0) / LFD_day0)) %>% #% change from start of injections 
  mutate(Reached_20 = if_else(BW_pct_change < -20), "Yes", "No")

#BW loss progress ####
#BW --> More than 15% decrease since start of injections?
BW_loss_goal <- df_LFD_19 %>%
  mutate(Lost_20 = if_else(BW_pct_change < -15, "Close", "No")) %>%
  filter(Lost_20 == "Close")

#---#
#Maybe do this later 
#Change in FI
#Rolling avg FI for 5 measurements prior to start of injections? Percent change?


