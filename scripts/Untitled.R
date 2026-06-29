#This script is a work in progress/attempt to do the following: 

#Calculate TEE (avg for a group) for each hour of a 24hr sable recording
#Distinguish between light and dark period
#Do this for each time point
  #baseline to peak obesity,
  #peak obesity to BW loss
  #BW los to BW maintenance
  #BW maintenance to BW regain
#Graph in multiple formats
  #histogram 
  #scatterplot with connecting line

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#TEE####
# build the summarized dataset
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  #summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep only one complete day (ie one complete 24hr period)
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE) %>% 
  #summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"))

sable_TEE_data <- sable_TEE_data %>%
  mutate(SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))

#echoMRI
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      Date == as.Date("2025-02-20") ~ "peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join lean mass info into TEE dataset
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean),
    by = c("ID", "SABLE")
  )

# Sort by SABLE, GROUP, hr (not by ID)
#Calculate mean kcal_hr for all mice for each hour I need to get the mean value for ad lib and the mean value for restricted
#then graph each hour
sable_TEE_adj2 <-sable_TEE_adj %>%
  ungroup() %>%
group_by(SABLE, GROUP, hr) %>%
mutate(hourly_mean=mean(value), .groups="drop") #maybe the .groups="drop" will eliminate IDs?

ggplot(sable_TEE_adj2, aes(x = SABLE, y = hourly_mean, color = GROUP, group = GROUP)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_text(aes(label = ID), size = 2.5, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", aes(group = GROUP), size = 1.2) +
  facet_wrap(~lights) +
  labs(y = "TEE (kcal/day)", color = "Group") +
  theme_minimal()

ggplot() +
  geom_jitter(data = df2, 
              aes(x = SABLE, y =value, color = hr),
              width = 0.2, alpha = 0.4, size = 2) +
  geom_point(data = emm_TEE_df,
             aes(x = SABLE, y = emmean, color = hr),
             position = position_dodge(0.2), size = 3) +
  geom_line(data = emm_TEE_df,
            aes(x = SABLE, y = emmean, color = hr, group = GROUP),
            position = position_dodge(0.2), linewidth = 1) +
  geom_errorbar(data = emm_TEE_df,
                aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                width = 0.1, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(y = "TEE (adjusted for Lean mass)", x = "SABLE phase",
       color = "Group",
       title = "TEE across SABLE phases (adjusted for Lean mass)") +
  theme(legend.position = "top")

df2 <- sable_TEE_adj %>%
  filter(SABLE == "bw regain")
 

plot <- df2 %>%
  ggplot(aes(x = hr, y = value, fill = GROUP)) +
  
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7, alpha = 0.7) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  
  #geom_point(aes(color = ),
            # position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
            # alpha = 0.6, size = 2) +
  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  theme_minimal() +
  labs(y = "Body Weight (g)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  format.plot
plot

####need a regression line here instead of a connecting line...might need to do a simple linear regression calculation####
plot1 <- sable_TEE_adj %>%
  filter(SABLE == "BW regain")
  
  ggplot() +
  geom_jitter(data = plot1, 
              aes(x = Lean, y = tee, color = GROUP),
              width = 0.2, alpha = 0.4, size = 2) +
  geom_point(data = plot1,
             aes(x = Lean, y = tee, color = GROUP),
             position = position_dodge(0.2), size = 3) +
  geom_line(data = plot1,
            aes(x = Lean, y = tee, color = GROUP, group = GROUP),
            position = position_dodge(0.2), linewidth = 1) +
  #geom_errorbar(data = emm_TEE_df,
               # aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                #width = 0.1, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(y = "TEE (adjusted for Lean mass)", x = "SABLE phase",
       color = "Group",
       title = "TEE across SABLE phases (adjusted for Lean mass)") +
  theme(legend.position = "top")
