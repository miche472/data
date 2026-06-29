####Aim is to modify figure5_paper to show light and dark periods
#1. At each time point
#2. divided into diet groups
#3. What statistical analyses do I want to do?
#percent of total movement occurring during the dark and light period. Maybe calculate for each mouse 
  #individually at each time point and then average at each time point within diet groups


#Libraries
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels
library(lme4)
library(stringr)
library(lme4)
library(car)  # car for Anova(), vif()
library(effsize) # for Cohen's d; install if needed

zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 


##mean of day 1 and 2 light collapsed ----

sable_loc_data <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Dark", "Light"),
    SABLE = case_when(
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                               "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                               "SABLE_DAY_7") ~ "Baseline",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  group_by(ID, DRUG,GROUP,DIET_FORMULA,SABLE) %>% #LM removed SEX and Strain
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,STRAIN,DIET_FORMULA,SABLE) %>% #LM removed SEX and STRAIN
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter( is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain")))

sable_loc_data_minutes <- sable_loc_data %>%
  #filter(SABLE =="BW regain") %>% 
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,GROUP,DIET_FORMULA,SABLE) %>% # group per animal, drug, and day #LM deleted SEX and STRAIN
  mutate(
    locomotion = value - lag(value),       # change in meters
    locomotion = if_else(locomotion < 0, 0, locomotion),  # remove negative jumps
    moving_min = if_else(locomotion > 0, 1, 0)             # 1 min per movement
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    # total minutes moved per day
    total_distance = sum(locomotion),      # total meters per day
    .groups = "drop"
  ) %>%
  group_by(ID,DRUG,GROUP,DIET_FORMULA,SABLE) %>% #LM deleted SEX and STRAIN
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  ) 

sable_loc_data_minutes%>%
  group_by(SABLE, GROUP) %>% #LM deleted strain
  summarise(n_ID = n_distinct(ID)) #this is good

sable_loc_data_minutes<- sable_loc_data_minutes %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("Vehicle", "RTIOXA_47")),
    #SEX = factor(SEX), #LM commented out SEX
    GROUP = factor(GROUP),
    DIET_FORMULA = factor(DIET_FORMULA),
    #STRAIN = factor(STRAIN), #LM commented out STRAIN
    SABLE = factor(SABLE))


  
#Format plot  
  scaleFill <- scale_fill_manual(values = c("#FAAC41", "#3498DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines  
  
#### Plot locomotion in meters for Control and Weight cycled at time points####
ggplot(sable_loc_data_minutes, aes(x = SABLE, y =  avg_distance , fill = GROUP)) +
  stat_summary(fun = mean, 
               geom = "col", 
               position = position_dodge(width = 0.8), 
               color = "black", 
               width = 0.7, 
               alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 0.8),
               width = 0.3) +
  geom_point(  
    aes(color = GROUP), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, 
    size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  format.plot+
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))+
  labs(title = "Locomotion (meters) in 24hrs", 
      y = "Locomotion (meters in 24h)", 
      x = "Time point",
      color= "Treatment group",
      fill= "Treatment group")


#### Evaluate if Control is different to Weight cycled in terms of meters/24h ####
#Data analysis...need to do. This is for 24hrs rather than distinguishing between light and dark


#### Time spent moving: Plot time spent moving (units?) ####
ggplot(sable_loc_data_minutes, aes(x = SABLE, y =  avg_moving_hr , fill = GROUP)) +
  #  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, 
               geom = "col", 
               position = position_dodge(width = 0.8),
               color = "black", 
               width = 0.7, 
               alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  geom_point(aes(group = GROUP),
             position = position_dodge(width = 0.8), 
             size = 2) +
  scaleFill +
  scaleColor +
  theme_minimal() +
  format.plot+
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Time spent moving (minutes)",
       y = "Time spent moving (hr in 24h)", 
       x = "Time point", 
       color= "Treatment group", 
       fill="Treatment group") +
  #geom_text_repel(aes(label = ID), size = 3, alpha = 0.7)

#data analysis to evaluate if ID 3708 is outlier

subset_data <- sable_loc_data_minutes %>%
  filter(GROUP == "Weight cycled")


ggplot(subset_data, aes(y = avg_moving_hr, x = "")) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
  geom_point(aes(label = ID), position = position_jitter(width = 0.1)) +
  geom_text_repel(aes(label = ID)) +
  labs(y = "Avg Moving Hours", x = "") +
  theme_minimal() #IT SEEMS LIKE ID 3708 IS AN OUTLIER

Q1 <- quantile(subset_data$avg_moving_hr, 0.25)
Q3 <- quantile(subset_data$avg_moving_hr, 0.75)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

subset_data %>%
  mutate(outlier_IQR = avg_moving_hr < lower_bound | avg_moving_hr > upper_bound) 

#conclusion 3708 IS an outlier so lets run the analysis without those ID

sable_loc_data_minutes<- sable_loc_data_minutes %>% 
  filter(!ID == 3708) 

# Plot time spent moving in min without ID 3708 ----
ggplot(sable_loc_data_minutes, aes(x = SABLE, y =  avg_moving_hr , fill = GROUP)) +
  stat_summary(fun = mean, 
               geom = "col", 
               position = position_dodge(width = 0.8),
               color = "black", 
               width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 0.8),
               width = 0.3) +
  geom_point(aes(color = GROUP), alpha = 0.7, size = 2, 
             position = position_dodge(width = 0.8)) +
  scaleFill +
  scaleColor +
  theme_minimal() +
  format.plot+
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))+
  labs(y = "Time spent moving (hr in 24h)", 
       x = "Time point",
       title = "Time spent moving in food restricted animals",
       fill="Treatment group",
       color="Treatment group") +
  #geom_text_repel(aes(label = ID), size = 3, alpha = 0.7)

  #### Evaluate if Control is different to Weight cycled in terms of time spent moving in 24hrs ####
#Data analysis...need to do. This is for 24hrs rather than distinguishing between light and dark



#### Divided into light and dark period ####
# mean of day 1 and 2 separated by lights----

sable_loc_data_lights <- sable_dwn %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "Off", "On"),
         SABLE = case_when(
           STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                                    "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                                    "SABLE_DAY_7") ~ "Baseline",
           STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
           STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
           STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )) %>%
  filter(!ID %in% c(3715,3712,3720,3721)) %>% #just testing what happens if we eliminate 20 and 21
  filter(grepl("AllMeters_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,lights,GROUP,DIET_FORMULA) %>% #LM deleted SEX and STRAIN
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,DRUG,lights,GROUP,DIET_FORMULA) %>% #LM deleted SEX and STRAIN
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("Baseline", "Peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))) 

sable_loc_data_minutes_lights <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,lights,GROUP,DIET_FORMULA,SABLE) %>%  #group per animal, drug, and day #LM deleted SEX, STRAIN; added SABLE
  mutate(
    locomotion = value - lag(value),       # change in meters
    locomotion = if_else(locomotion < 0, 0, locomotion),  # remove negative jumps
    moving_min = if_else(locomotion > 0, 1, 0)             # 1 min per movement
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    # total minutes moved per day
    total_distance = sum(locomotion),      # total meters per day
    .groups = "drop"
  ) %>%
  group_by(ID, DRUG,lights,GROUP,DIET_FORMULA, SABLE) %>% #LM deleted SEX, STRAIN and added SABLE
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  )

sable_loc_data_minutes_lights%>%
  group_by(lights) %>% #LM deleted SEX and STRAIN
  summarise(n_ID = n_distinct(ID)) #this is good. verify number of mice is the same during dark and light

sable_loc_data_minutes_lights<- sable_loc_data_minutes_lights %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("Vehicle", "RTIOXA_47")),
    #SEX = factor(SEX),
    GROUP = factor(GROUP),
    DIET_FORMULA = factor(DIET_FORMULA)
    #STRAIN = factor(STRAIN)
  )

# Plot total locomotion in meters LIGHTS ON AND OFF ----
ggplot(sable_loc_data_minutes_lights, aes(x = SABLE, y =  avg_distance , fill = GROUP)) +
  stat_summary(fun = mean, 
               geom = "col", 
               position = position_dodge(width = 0.8),
               color = "black", 
               width = 0.7, 
               alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 0.8),
               width = 0.3) +
  geom_point(aes(color = GROUP), 
             alpha = 0.7, 
             size = 2,   
             position = position_dodge(width = 0.8)) +
  scaleFill +
  scaleColor +
  theme_minimal() +
  theme(legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title ="title",
       y = "Total locomotion (meters in 24h)",
       x = "Time point",
       color="Treatment group",
       fill="Treatment group") +
  facet_wrap(~lights)


## Plot time spent moving in min LIGHTS ON AND OFF ----
ggplot(sable_loc_data_minutes_lights, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  #  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(lights~ GROUP) + #LM deleted SEX and STRAIN
  scale_fill_manual(values = c(
    "Vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "Time spent moving (hr in 24h)", x = "")+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)

####separate based on lights and sable phase
ggplot(sable_loc_data_minutes_lights, aes(x = SABLE, y =  avg_moving_hr , fill = SABLE)) +
  #  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(lights~ GROUP) + #LM deleted SEX and STRAIN
  scale_fill_manual(values = c(
    "Vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "Time spent moving (hr in 24h)", x = "")+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)

####just wrap by light or dark period
ggplot(sable_loc_data_minutes_lights, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  #  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(lights~ GROUP) + #LM deleted SEX and STRAIN
  scale_fill_manual(values = c(
    "Vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "Time spent moving (hr in 24h)", x = "")+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)

#### SABLE phase analysis ####
#NOte that 3712, 3715, and 3708 have been removed from sable_loc_data_minutes because 3708 was an outlier
sable_loc_data_minutes%>%
  group_by(SABLE,GROUP) %>%
  summarise(n_ID = n_distinct(ID)) #this is good



## Plot locomotion in meters ----
#Format plot
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#3498DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  #panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines


ggplot(sable_loc_data_minutes, aes(x = SABLE, y = avg_distance , fill = GROUP)) + 
  stat_summary( # mean bars
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  stat_summary( # error bars (mean ± SE)
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  geom_point(  # individual data points
    aes(color = GROUP), 
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + 
  scaleColor +
  theme_minimal() +
  labs(title = "NZO locomotion", 
       y = "locomotion (meters in 24hrs)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))
