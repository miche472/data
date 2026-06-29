#Calculate 24hr FI 

#For each mouse calculate 24hr FI during the Sable measurement points
  #Baseline, peak obesity, BW loss, BW maintenance, BW regain
#This script calculates 24hr FI for mice during Baselin, peak obesity & BW regain
#Pending: 
  #Remove calculated FI for GROUP=restricted during SABLE=BW loss, BW maintenance 
  #Add daily FI for GROUP=restricted during SABLE=BW loss, BW maintenance

library(dplyr)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(ggpubr)
install.packages("write_csv")

# --- Function ---
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr - 20, hr + 4))
}

# --- Load data ---
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds")
#META <- read.csv(file = "../data/META.csv")
#FI <-read.csv(file = "../data/FI.csv")

# --- Data prep for FI measured by Sable food hopper ---
sable_FI_LM_data <- sable_dwn %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  mutate(SABLE = case_when(
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                             "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                             "SABLE_DAY_7") ~ "baseline",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3","SABLE_DAY_4") ~ "peak obesity",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7","SABLE_DAY_8") ~ "BW loss",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12") ~ "BW maintenance", 
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain"
  )           
  ) %>%          
  filter(!is.na(SABLE)) %>%
  
  filter(parameter=="FoodA") %>%
  mutate(food_mass= value) %>%
  group_by(ID, SABLE,STRAIN,SEX,DIET_FORMULA) %>%
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>%
  ungroup() %>%
  group_by(ID, complete_days, STRAIN, SEX, DIET_FORMULA) %>%
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>%
  ungroup() %>%
  
  #Calculate the minute by minute changes in food intake. Sum all of these minute by minute changes to get FI for the day
  group_by(ID, complete_days, SABLE,STRAIN, SEX,DIET_FORMULA) %>%
  mutate(FI_minute = lag(value) - value) %>% 
  #Mass at one minute previous minus mass at the current minute (ex. mass at 11:00a - mass at 11:01am)
  filter(FI_minute >= 0) %>% #this corrects for instances when fresh food was added to the hopper in Sable
  summarise(FI_total = sum(FI_minute), .groups = "drop") %>%
  
#### Can't figure out how to do this (issue 1) ####
  #For GROUP=restricted mice during SABLE=BW loss and SABLE=BW maintenance I need to remove the 
  #FI recorded by SABLE (recall that restricted mice did not eat out of the Sable food hopper)
  #The daily mass of food consumed by these mice is located in the FI.csv
  
  # exclude problematic IDs (same as TEE)
  filter(!ID %in% c(3712, 3715)) %>% # NZO mice which died during study
  
#### Can't figure out how to do this (issue 2) ####
  #When SABLE= BW loss and BW maintenance exclude the FI measured by Sable for the GROUP=restricted mice 
  #filter(!ID %in% c(3708, 3710, 3714, 3720, 3721, 3722, 3723, 3724, 3725, 3727, 3728, 3729, 7861,
  #7863,7865,7878,7872,7874,7877,,7866)) %>% 
  
  filter(complete_days %in% c(1, 2)) %>%
  
  # average across 2 complete days
  group_by(ID, SABLE,STRAIN,SEX,DIET_FORMULA) %>%
  summarise(Day_totalFI = mean(FI_total), .groups = "drop") %>%
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728, 7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729, 7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )
  ) %>%
  mutate(
    SABLE = factor(SABLE,
                   levels = c("baseline", "peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))
  ) %>% 
  filter(!is.na(SABLE)) %>% 
  filter(SEX=="F") 

write.csv(sable_FI_LM_data, "../data/sable_FI_LM_data.csv") # Save as CSV

#----Mixed model for 24hr FI (Strain, Sable, Group, Drug)
FI_model <- lmer(Day_totalFI ~ STRAIN* SABLE * GROUP * DRUG  + (1|ID), data = sable_FI_LM_data)
summary(FI_model)

# --- Estimated marginal means for 24hr FI ---
FI_emm <- emmeans(FI_model, pairwise ~ STRAIN* SABLE * GROUP * DRUG, adjust = "tukey")

# --- Convert to dataframes for results for 24hr FI ---
df_emm_FI <- as.data.frame(FI_emm$emmeans)
df_pairs_FI <- as.data.frame(FI_emm$contrasts)

# --- Keep significant results only for 24hr FI ---
df_sig_FI <- df_pairs_FI %>%
  filter(p.value <= 0.05)
print(df_sig_FI, n = Inf)
View(df_sig_FI)

#Other linear mixed models to create
  #lmer(Day_totalFI ~ STRAIN* SABLE * GROUP  + (1|ID), data = sable_FI_LM_data)
  #lmer(Day_totalFI ~ STRAIN* SABLE + (1|ID), data = sable_FI_LM_data)
