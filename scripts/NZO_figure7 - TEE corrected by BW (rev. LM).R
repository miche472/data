#we aim to evaluate changes in 24h TEE in middle age NZO after different stages of feeding:
#This script is TEE adjusted for BW (grams)
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

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
                     "SABLE_DAY_7") ~ "Baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "Peak obesity",
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
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE 
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "Control",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "Weight cycled"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "Vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("Baseline", 
                              "Peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )

#echoMRI
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
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
  dplyr::select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "Baseline",
      Date == as.Date("2025-02-20") ~ "Peak obesity",
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
                         levels = c("Baseline", "Peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join lean mass info into TEE dataset
sable_TEE_adj_BW <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% 
      dplyr::select(ID, SABLE, Lean, Weight),
    by = c("ID", "SABLE")
  )

#### Build linear mixed model for TEE (adj for Body weight) #### 
model_TEE_BW1 <- lmer(tee ~ SABLE * GROUP + Weight + (1 | ID), data = sable_TEE_adj_BW)
summary(model_TEE_BW1)

#Confirm the number of mice in the data frame is correct
n_distinct(sable_TEE_adj_BW$ID) #good we have 16 animals

#### Calculate estimated marginal means (emmeans) ####
emm_TEE_BW1 <- emmeans(model_TEE_BW1, ~ SABLE * GROUP, cov.reduce = mean)
emm_TEE_BW1_df <- as.data.frame(emm_TEE_BW1)

# Pairwise contrasts within each GROUP
      contrasts_by_group_TEE_BW1 <- contrast(emm_TEE_BW1, method = "pairwise", by = "GROUP")
      # Convert to a data frame
      contrasts_TEE_BW1_df <- as.data.frame(contrasts_by_group_TEE_BW1)
      
      # Filter for restricted group and baseline vs __
      restricted_contrasts_TEE_BW1 <- contrasts_TEE_BW1_df %>%
        filter(GROUP == "Weight cycled") %>%
               filter(contrast %in% c("Baseline - BW loss", "Baseline - BW maintenance", "Baseline - BW regain"))
      restricted_contrasts_TEE_BW1
      
      #Filter for significant contrasts
      Sig_contrasts_TEE_BW1 <- contrasts_TEE_BW1_df %>%
        filter(p.value <=0.05)
      Sig_contrasts_TEE_BW1

# Pairwise contrasts within each SABLE (time point)
      contrasts_by_SABLE_TEE_BW1 <- contrast(emm_TEE_BW1, method = "pairwise", by = "SABLE")
      # Convert to a data frame
      contrasts_SABLE_TEE_BW1_df <- as.data.frame(contrasts_by_SABLE_TEE_BW1)

#### Plots ####
#Format plot (optional)
  format.plot <- theme(
    strip.background = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),          
    panel.spacing.y = unit(1.5, "lines"),  
    axis.text = element_text(family = "Helvetica", size = 13),
    axis.title = element_text(family = "Helvetica", size = 14),
    #panel.grid.major = element_blank(), # remove background grid lines only
    panel.grid.minor = element_blank(), # remove background grid lines only
    axis.line = element_line(color = "black")) # keep axis lines

####Scatter plot - Graph predicted TEE adjusted for BW (NZO) ####
#Commented out the measured values for TEE

ggplot() +
  #geom_jitter(data = sable_TEE_adj_BW, 
             #aes(x = SABLE, y = tee, color = GROUP),
              #width = 0.2, alpha = 0.4, size = 2) +
  geom_point(data = emm_TEE_BW1_df,
             aes(x = SABLE, y = emmean, color = GROUP),
             position = position_dodge(0.2), size = 4) +
  geom_line(data = emm_TEE_BW1_df,
            aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
            position = position_dodge(0.2), linewidth = 1.5) + 
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  geom_errorbar(data = emm_TEE_BW1_df,
                aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                width = 0.15, position = position_dodge(0.2)) +
  theme_minimal(base_size = 14) +
  labs(y = "Adjusted TEE (kcal/day)", x = "Time point",
       color = "Restriction group",
       title = "Total energy expenditure (TEE) adj. for bodyweight") +
    format.plot +
  theme(legend.position = "top", 
        plot.title = element_text(hjust=0.5), 
        axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#### Bar plot - Graph predicted TEE adjusted for BW ####
#Commented out the measured values for TEE

# Define custom colors
custom_colors <- c(
  "Control" = "#FAAC41",              
  "Weight cycled" = "#3498DB")

barplot_emm_TEE_BW1 <- emm_TEE_BW1_df %>%
  ggplot(aes(x = SABLE, y = emmean, fill = GROUP)) +
  # mean bars
  geom_col(position = position_dodge(width = 0.8),
           color = "black", width = 0.7, alpha = 0.7) +
  # error bars using SE
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    position = position_dodge(width = 0.8),
    width = 0.3) +
  #scale_color_manual(values=c('#FAAC41','#5392DB'))+
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Total energy expenditure adjusted for BW",
    y = "Adjusted TEE (kcal/day)",
    x = "Time point",
    fill = "Treatment group") +
  format.plot +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
barplot_emm_TEE_BW1
