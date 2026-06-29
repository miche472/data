#Process BW and FI data for COHORT_19

#Started:3-31-26
#Last revised: 4-16-26 
#Last run on: 4-16-26

#I decided to use adiposity index to sort mice into treatment groups
#Data frame with mice assigned to groups is called 

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
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>%
  filter(COHORT == 19)

#Join FI and BW
BW_FI_19 <- FI_LM_COHORT19 %>%
  left_join(
    BW_COHORT19 %>% 
      select(ID, BW, DATE),
    by = c("ID", "DATE")) %>%
  mutate(ID = as.factor(ID)) %>%
  filter(!ID =='3734') %>% #died during study
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(day_rel = DATE - first(DATE),
         day_rel = as.numeric(day_rel)) #%>%
 # replace_na(list(INTAKE_GR=0, delta_alt=0, corrected_intake_gr=0,
                  #corrected_intake_kcal=0, KCAL_G=3.82))
  

# Create BW & FI graphs ####
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

##BW ####
#Graph BW
BW_FI_19_plot <- BW_FI_19 %>%
  filter(BW>0) #remove rows where there isn't a BW measurement
  
ggplot(BW_FI_19_plot,
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
#It may be best to use BW_FI_19 here and drop rows (dates) for which corrected_intake_kcal is NA
FI_plot_BW_FI_19 <- BW_FI_19 %>%
  drop_na(corrected_intake_gr)

ggplot(FI_plot_BW_FI_19, aes(x = day_rel, y = corrected_intake_gr, color = ID, fill = ID)) +
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
       y= "Energy intake (grams/day)",
       title= "NZO Daily energy intake (grams/day)")


#Assign to treatment groups ####
##1st group of NZO mice (n=11) ####
#Try using multiple methods. (either adiposity index or BW) ##
#Drug groups: TZP (10nmol/kg), SVD (30 nmol/kg), vehicle (1X PBS, 5 ml/kg)

#After assigning groups, interrogate group assignation methods 
    #1. If I use BW for assignation, do the resulting groups have different average 
        #adiposity indices?
    #2. If I use adiposity index for assignation, do the resulting groups have 
        #different avg BW?
        
#---
###Method 1: Use adiposity index to assign mice to groups ####
#Need to save the echoMRI data files with the correct formatting and title
#to the echoMRI folder in the GitHub data folder. Then I need to re-run
#data_proc so that an updated echoMRI.csv is written.
#The code below will pull data from echoMRI.csv

#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_dataA<- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  filter(Date == "2026-04-16") %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index) %>%
  ungroup()

#Create pattern for group assignment
pattern1 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

echoMRI_dataB <- echoMRI_dataA %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern1, length.out = n())) %>%
  ungroup()

#I will use adiposity index (method 1) for group assignment (decided on 4/16/26)
#Have R randomly assign TZP and SVD. Make the group with n=3 be vehicle. 

#Assign the group of 3 mice to be vehicle
group_sizes <- echoMRI_dataB %>% # Count how many mice per group
  count(treatment)
group_3 <- group_sizes %>% # Find the group with 3 mice
  filter(n == 3) %>%
  pull(treatment)
groups_4 <- group_sizes %>% # Find the groups with 4 mice
  filter(n == 4) %>%
  pull(treatment)

#Create a random mapping between treatment groups and drugs
# - 3 mice → Vehicle
# - 4 mice → randomly TZP or SVD
Group_to_drug_map <- c(
  setNames("Vehicle", group_3),
  setNames(sample(c("Tirzepatide", "Survodutide")), groups_4))

# Apply the mapping to my df called echoMRI_dataB
echoMRI_dataC <- echoMRI_dataB %>%
  mutate(DRUG = Group_to_drug_map[treatment])

#--
#One way ANOVA for AI sorted: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW <- aov(Weight ~ treatment, data = echoMRI_dataC)
summary(ANOVA_AI_sorted_checkBW)

  #Perform Tukey's Test for multiple comparisons
  TukeyHSD(ANOVA_AI_sorted_checkBW, conf.level=.95) 
  
  #Shapiro Wilk test to check normality of residuals
  shapiro.test(residuals(ANOVA_AI_sorted_checkBW))
  #Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
  #the p value is 0.3557 suggesting that they do come from a normal distribution
  
  #Check for equal variance
  #Create box plots that show distribution of adiposity index for each group
  boxplot(Weight ~ treatment, xlab='Treatment group', ylab='Body weight', data=echoMRI_dataC)
  #Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
  #Check for equal variances
  bartlett.test(Weight ~ treatment, data=echoMRI_dataC)
  #p value is 0.8456 --> samples have equal variance
         
#--
####ANOVAs for AI sorted ####
#One way ANOVA for AI sorted: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkAI <- aov(adiposity_index ~ treatment, data = echoMRI_dataC)
summary(ANOVA_AI_sorted_checkAI)

  #Perform Tukey's Test for multiple comparisons
  TukeyHSD(ANOVA_AI_sorted_checkAI, conf.level=.95) 
  
  #Shapiro Wilk test to check normality of residuals
  shapiro.test(residuals(ANOVA_AI_sorted_checkAI))
  #Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
  #the p value is 0.4669 suggesting that they do come from a normal distribution
  
  #Check for equal variance
  #Create box plots that show distribution of adiposity index for each group
  boxplot(adiposity_index ~ treatment, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC)
  #Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
  #Check for equal variances
  bartlett.test(adiposity_index ~ treatment, data=echoMRI_dataC)
  #p value is 0.8058 --> samples have equal variance

#---
###Method 2: Use BW to assign mice to groups ####
#Use echoMRI_data which was read in for method 1 --> 
          #for method 2 I will use body "Weight" included in this df
Weight_dataA <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  #filter(Date == "2026-04-16") %>%
  select(ID, adiposity_index, Weight, Date) %>%
  ungroup()

#Create pattern for group assignment
#Have R randomly assign TZP, SVD, and Vehicle to Groups 1-3
pattern1 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

Weight_dataB <- Weight_dataA %>%
  arrange(desc(Weight)) %>%
  mutate(treatment = rep(pattern1, length.out = n()))

#--
#### ANOVAs for BW sorted ####
#One way ANOVA: same mean BW for all three groups?
ANOVA_BW_sorted_checkBW <- aov(Weight ~ treatment, data = Weight_dataB)
summary(ANOVA_BW_sorted_checkBW)

  #Perform Tukey's Test for multiple comparisons
  TukeyHSD(ANOVA_BW_sorted_checkBW, conf.level=.95) 
  
  #Shapiro Wilk test to check normality of residuals
  shapiro.test(residuals(ANOVA_BW_sorted_checkBW))
  #Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
  #the p value is 0.8036 suggesting that they do come from a normal distribution
  
  #Check for equal variance
  #Create box plots that show distribution of adiposity index for each group
  boxplot(Weight ~ treatment, xlab='Treatment group', ylab='Body weight', data=Weight_dataB)
  #Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
  #Check for equal variances
  bartlett.test(Weight ~ treatment, data=Weight_dataB)
  #p value is 0.7189 --> samples have equal variance

#--
#One way ANOVA: same mean adiposity index for all three groups?
ANOVA_BW_sorted_checkAI <- aov(adiposity_index ~ treatment, data = Weight_dataB)
summary(ANOVA_BW_sorted_checkAI)

  #perform Tukey's Test for multiple comparisons
  TukeyHSD(ANOVA_BW_sorted_checkAI, conf.level=.95) 

  #Shapiro Wilk test to check normality of residuals
  shapiro.test(residuals(ANOVA_BW_sorted_checkAI))
  #Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
  #the p value is 0.6483 suggesting that they do come from a normal distribution

  #Check for equal variance
  #Create box plots that show distribution of adiposity index for each group
  boxplot(adiposity_index ~ treatment, xlab='Treatment group', ylab='Adiposity index', data=Weight_dataB)
  #Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
  #equal variances
  bartlett.test(adiposity_index ~ treatment, data=Weight_dataB)
  #p value is 0.3856 --> samples have equal variance

#---
#Calculate and graph mean BW, AI, and rolling average FI for all mice
#when using the adiposity index sorting method --> does it look like the groups are fairly even?
###Graph treatment groups BW and AI ####
  
#Sorted into groups by adiposity index. Graph of mean adiposity index by treatment group
ggplot(echoMRI_dataC, aes(x = treatment, y = adiposity_index)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        panel.border = element_blank()) +
  labs(title="Sorted by adiposity index", y = "Adiposity index)")

#Sorted into groups by adiposity index. Graph of mean BW by treatment group
ggplot(echoMRI_dataC, aes(x = treatment, y = Weight)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.text.x = element_text(color="black"),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(),
        axis.title.x = element_text(face="bold"),
        panel.border = element_blank()) +
  labs(title="Weight (Sorted by adiposity index)", 
       x="Treatment group ", y = "Body Weight (g)")

#Sorted into groups by Weight. Graph of mean BW by treatment group
ggplot(Weight_dataB, aes(x = treatment, y = Weight)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
        axis.title.y = element_text(size= 20, face="bold"),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  labs(title="Sorted by Weight", y = "Weight (g)")

#Sorted into groups by adiposity index. Graph of mean adiposity index by treatment group
ggplot(Weight_dataB, aes(x = treatment, y = adiposity_index)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8), width=0.73) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.65, color="#454441") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.text.x = element_text(size= 20, color="black", angle=20, hjust=1, vjust=1),
        axis.title.y = element_text(size= 20, face="bold"),
        axis.text.y = element_text(size = 20),
        panel.border = element_blank()) +
  labs(title="Sorted by Weight", y = "Adiposity index (fat/lean)")

### FINAL df for first n= 11 NZO mice ####
df_incretin_groups <- echoMRI_dataC %>%
  select(-treatment) %>%
  ungroup()
  
#-----------------------------------------------#
#-----------------------------------------------#
#-----------------------------------------------#
##2nd group of NZO mice (n=7) ####

###Method 1: Use adiposity index to assign mice to groups ####

#Need to save the echoMRI data files with the correct formatting and title
#to the echoMRI folder in the GitHub data folder. Then I need to re-run
#data_proc so that an updated echoMRI.csv is written.
#The code below will pull data from echoMRI.csv

#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_dataA2<- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  filter(ID %in% c(3747, 3748,3749, 3750, 3751, 3752, 3753)) %>%
  filter(Date == "2026-04-29") %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index) %>%
  ungroup()

#Create pattern for group assignment
pattern12 <- c("Group1", "Group2", "Group3", "Group3", "Group2", "Group1")

echoMRI_dataB2 <- echoMRI_dataA2 %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern12, length.out = n())) %>%
  ungroup()

#Use adiposity index (method 1) for group assignment 
#R randomly assign: Tirzepatide (n=3), Survodutide (n=2), and Vehicle (n=2). 

#Assign the group of 3 mice to be Tirzepatide
group_sizes <- echoMRI_dataB2 %>% # Count how many mice per group
  count(treatment)
groups_3 <- group_sizes %>% # Find the group with 3 mice
  filter(n == 3) %>%
  pull(treatment)
groups_2 <- group_sizes %>% # Find the groups with 2 mice
  filter(n == 2) %>%
  pull(treatment)

#Create a random mapping between treatment groups and drugs
# - 3 mice → Tirzepatide
# - 2 mice → randomly Survodutide or Vehicle
Group_to_drug_map <- c(
  setNames("Tirzepatide", groups_3),
  setNames(sample(c("Survodutide", "Vehicle")), groups_2))

# Apply the mapping to my df called echoMRI_dataB2
echoMRI_dataC2 <- echoMRI_dataB2 %>%
  mutate(DRUG = Group_to_drug_map[treatment])

#--
#One way ANOVA for AI sorted: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW_2 <- aov(Weight ~ treatment, data = echoMRI_dataC2)
summary(ANOVA_AI_sorted_checkBW_2)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkBW_2, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkBW_2))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.9842 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(Weight ~ treatment, xlab='Treatment group', ylab='Body weight', data=echoMRI_dataC2)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(Weight ~ treatment, data=echoMRI_dataC2)
#p value is 0.7083 --> samples have equal variance

#--
####ANOVAs for AI sorted ####
#One way ANOVA for AI sorted: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkAI_2 <- aov(adiposity_index ~ treatment, data = echoMRI_dataC2)
summary(ANOVA_AI_sorted_checkAI_2)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkAI_2, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkAI_2))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.2491 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(adiposity_index ~ treatment, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC2)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(adiposity_index ~ treatment, data=echoMRI_dataC2)
#p value is 0.729 --> samples have equal variance


###Method: Try doing a different assignment method ####
pattern_new <- c("Group1", "Group2", "Group3", "Group3", "Group3", "Group2", "Group1")

echoMRI_data_Bnew <- echoMRI_dataA2 %>%
  arrange(desc(adiposity_index)) %>%
  mutate(treatment = rep(pattern_new, length.out = n())) %>%
  ungroup()

#Use adiposity index (method 1) for group assignment 
#R randomly assign: Tirzepatide (n=3), Survodutide (n=2), and Vehicle (n=2). 

#Assign the group of 3 mice to be Tirzepatide
group_sizes <- echoMRI_data_Bnew %>% # Count how many mice per group
  count(treatment)
group_3 <- group_sizes %>% # Find the group with 3 mice
  filter(n == 3) %>%
  pull(treatment)
groups_2 <- group_sizes %>% # Find the groups with 2 mice
  filter(n == 2) %>%
  pull(treatment)

#Create a random mapping between treatment groups and drugs
# - 3 mice → Tirzepatide
# - 2 mice → randomly Vehicle or Survodutide
Group_to_drug_map <- c(
  setNames("Tirzepatide", group_3),
  setNames(sample(c("Vehicle", "Survodutide")), groups_2))

# Apply the mapping to my df called echoMRI_dataB
echoMRI_data_Cnew <- echoMRI_data_Bnew %>%
  mutate(DRUG = Group_to_drug_map[treatment])

#--
#One way ANOVA for AI sorted: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW_new <- aov(Weight ~ treatment, data = echoMRI_data_Cnew)
summary(ANOVA_AI_sorted_checkBW_new)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkBW_new, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkBW_new))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.8829 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(Weight ~ treatment, xlab='Treatment group', ylab='Body weight', data=echoMRI_data_Cnew)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(Weight ~ treatment, data=echoMRI_data_Cnew)
#p value is 0.3285 --> samples have equal variance

#Combine 1st and 2nd ####
#Which approach (2 or new) gives more even groups when combined with the n=11 mice?
#Use this to make the decision

##Approach "2" (echo_MRI_C2) ####
#Compile dfs echoMRI_dataC and echoMRI_dataC2 --> now n=18
echoMRI_dataC_joined_2 <- bind_rows(echoMRI_dataC, echoMRI_dataC2)
check_echoMRI_dataC_joined_2 <- bind_rows(echoMRI_dataC, echoMRI_dataC2)
boxplot(Weight ~ DRUG, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC_joined_2)

####ANOVAs for AI sorted ####
#One way ANOVA for AI sorted: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkAI_joined2 <- aov(adiposity_index ~ DRUG, data = echoMRI_dataC_joined_2)
summary(ANOVA_AI_sorted_checkAI_joined2)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkAI_joined2, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkAI_joined2))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.2859 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(adiposity_index ~ DRUG, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC_joined_2)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(adiposity_index ~ DRUG, data=echoMRI_dataC_joined_2)
#p value is 0.9162 --> samples have equal variance

#---
#One way ANOVA for AI sorted: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW_joined2 <- aov(Weight ~ DRUG, data = echoMRI_dataC_joined_2)
summary(ANOVA_AI_sorted_checkBW_joined2)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkBW_joined2, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkBW_joined2))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.5384 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(Weight ~ DRUG, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC_joined_2)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(Weight ~ DRUG, data=echoMRI_dataC_joined_2)
#p value is 0.7211 --> samples have equal variance


##Approach "new" (echoMRI_data_Cnew) ####
#Compile dfs echoMRI_dataC and echoMRI_data_Cnew --> now n=18
echoMRI_dataC_joined_new <- bind_rows(echoMRI_dataC, echoMRI_data_Cnew)

####ANOVAs for AI sorted ####
#One way ANOVA for AI sorted: same mean adiposity index for all three groups?
ANOVA_AI_sorted_checkAI_joinednew <- aov(adiposity_index ~ DRUG, data = echoMRI_dataC_joined_new)
summary(ANOVA_AI_sorted_checkAI_joinednew)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkAI_joinednew, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkAI_joinednew))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.4618 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(adiposity_index ~ DRUG, xlab='Treatment group', ylab='Adiposity index (fat/lean)', data=echoMRI_dataC_joined_new)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(adiposity_index ~ DRUG, data=echoMRI_dataC_joined_new)
#p value is 0.5074 --> samples have equal variance

#--
#One way ANOVA for AI sorted: same mean BW for all three groups?
ANOVA_AI_sorted_checkBW_new <- aov(Weight ~ DRUG, data = echoMRI_data_Cnew)
summary(ANOVA_AI_sorted_checkBW_new)

#Perform Tukey's Test for multiple comparisons
TukeyHSD(ANOVA_AI_sorted_checkBW_new, conf.level=.95) 

#Shapiro Wilk test to check normality of residuals
shapiro.test(residuals(ANOVA_AI_sorted_checkBW_new))
#Shapiro Wilk test tests the null hypothesis that the samples come from a normal distribution
#the p value is 0.8829 suggesting that they do come from a normal distribution

#Check for equal variance
#Create box plots that show distribution of adiposity index for each group
boxplot(Weight ~ DRUG, xlab='Treatment group', ylab='Body mass (g)', data=echoMRI_data_Cnew)
#Conduct Bartlett's test for variance -> tests the null hypothesis that the samples have 
#Check for equal variances
bartlett.test(Weight ~ DRUG, data=echoMRI_data_Cnew)
#p value is 0.3285 --> samples have equal variance

#--


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
  sample(c("Tirzepatide", "Survodutide", "Vehicle")),
  c("Group1", "Group2", "Group3"))

# Apply the mapping to my df
BW_dataC <- BW_dataB %>%
  mutate(DRUG = Group_to_drug_map[treatment])

#One way ANOVA: same mean BW for all three groups?
ANOVA_BW_sorted_checkBW <- aov(BW ~ DRUG, data = BW_dataC)
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
  filter(Dropped_10_percent == "Yes")

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
  mutate(avg_vs_current = moving_avg - corrected_intake_gr)%>%
  #filter(DATE > '2026-04-3') %>% #to see from beginning of LFD
  filter(DATE > '2026-04-08') %>% #to see most recent
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
  mutate(Lost_20_percent = if_else(BW_pct_change < -15, "Close", "No")) %>%
  filter(Lost_20_percent == "Close")

#---#
#Maybe do this later 
#Change in FI
#Rolling avg FI for 5 measurements prior to start of injections? Percent change?
