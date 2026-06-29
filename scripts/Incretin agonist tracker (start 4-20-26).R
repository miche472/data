#Track BW and FI of mice once they start injections


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

#Calculate cumulative FI
FI_19 <- BW_FI_19 %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(INTAKE_GR = if_else(INTAKE_GR >= 0, INTAKE_GR, 0, missing=0),
         cum_INTAKE_GR= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR), 0, missing=0),
         cum_INTAKE_kcal= if_else(INTAKE_GR >=0, cumsum(INTAKE_GR*3.82), 0, missing=0))

#Make BW gain, BW loss, and BW regain time periods
df_19 <- BW_FI_19 %>%
mutate(Inject_day= case_when(COMMENTS=="INJECT_DAY_ONE"~1,
                               COMMENTS=="INJECT_DAY_TWO"~2,
                               COMMENTS=="INJECT_DAY_THREE"~3,
                               COMMENTS=="INJECT_DAY_FOUR"~4)) %>%
mutate(stage= case_when(ID %in% c(3730, 3731, 3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740) 
& DATE< "2026-04-17" ~ "BW gain")) 
  
df_incretin_groups <- echoMRI_dataC %>%
  select(-treatment) %>%
  ungroup()

#Join DRUG groups from df_incretin_groups to BW_FI_19
GLP1_BWloss <-BW_FI_19  %>%
  left_join(echoMRI_dataC %>% 
      select(ID, DRUG),by = c("ID")) %>%
  mutate(DRUG = as.factor(DRUG)) %>%
  mutate(Inject_day= case_when(COMMENTS=="INJECT_DAY_ONE"~1,
                               COMMENTS=="INJECT_DAY_TWO"~2,
                               COMMENTS=="INJECT_DAY_THREE"~3,
                               COMMENTS=="INJECT_DAY_FOUR"~4)) %>%
  filter(Inject_day >0)
  


##BW ####
#Graph BW
ggplot(GLP1_BWloss,
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


#BW by group
ggplot(GLP1_BWloss, aes(x = Inject_day, y = BW, fill = DRUG, color=DRUG)) +
  geom_line(stat = "summary", 
           fun = "mean") +
  geom_point()

#BW percent change by drug group
#Make data frame with percent change
GLP1_BW_percent <- GLP1_BWloss %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(BW_pct_change = 100*((BW - first(BW)) / first(BW)))

#Make graph of percent change
ggplot(GLP1_BW_percent, aes(x = Inject_day, y = BW_pct_change, fill = DRUG)) +
    geom_point(stat = "summary", 
               fun = "mean") +
  geom_line()+
    labs(x="# of days with injections",
         y= "% change in BW",
         title= "Percent change in BW during injections")

library(dplyr)
library(ggplot2)

# Calculate averages
df_mean <- GLP1_BW_percent %>%
  group_by(DRUG, Inject_day) %>%
  summarise(mean_val = mean(BW_pct_change)) 

# Plot
ggplot(df_mean, aes(x = Inject_day, y = mean_val), fill=="DRUG") +
  geom_point(stat = "summary", 
             fun = "mean") +
  labs(x="# of days with injections",
       y= "% change in BW",
       title= "Percent change in BW during injections")


#### Left off here ####
ggplot(GLP1_BW_percent, aes(x=Inject_day, y=BW_pct_change, group=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG)) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG)) +
geom_errorbar(stat = "summary", 
              fun.data = mean_se, aes(width=0.1))







