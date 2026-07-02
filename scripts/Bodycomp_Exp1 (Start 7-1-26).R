# EchoMRI data after BW loss for NZO mice
#requested in 7-1-26 lab meeting by Laurie and Cathy
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
custom_colors_GLP <- c("Tirzepatide" = "#1e6deb","Survodutide" = "#c93618", "Vehicle" = "#403d3c", "Naive" ="darkgray")
custom_colors_GLPB <- c("Tirzepatide" = "lightblue","Survodutide" = "lightgreen", "Vehicle" = "darkgray")

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo,
  robustlmm,
  mmand)

#read in metadata
metadata <- read_csv("../data/META.csv") %>% 
select(ID, SEX, COHORT, STRAIN, AIM, DIET_FORMULA) %>% 
  mutate(ID=as.factor(ID))

# echo MRI
echomri_csv_files <- tibble(
  filepath = list.files("../data/echoMRI", full.names = TRUE)) %>% 
  filter(grepl("*.xlsx", filepath)) 
echomri_csv_files

echomri_open_files <- echomri_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map(., function(X){
    readxl::read_xlsx(X$filepath) %>% 
      select(Label, Fat, Lean, Weight, TimeDateDura) %>% 
      rename(ID = Label) %>% 
      separate_wider_delim(TimeDateDura, delim = ";", names = c("Date", "A", "B")) %>% 
      select(-A, -B) %>% 
      separate_wider_delim(Date, delim = " ", names = c("hms", "month", "day", "year")) %>% 
      mutate(day = gsub(",", "", day),
             Date = paste(year, month, day, sep = "-"),
             Date = lubridate::ymd(Date),
             ID =  as.factor(ID)) %>% 
      select(-hms, -month, -day, -year)}) %>% 
  bind_rows() %>% 
  left_join(., metadata, by = "ID")
echomri_open_files

# compare adiposity index = fat / lean
echomri_data <- echomri_open_files %>% 
  mutate(adiposity_index = Fat / Lean) %>% 
  group_by(ID) %>% 
  mutate(n_measurement = as.numeric(as.factor(Date)))
echomri_data

write_csv(x = echomri_data, "../data/echomri.csv")

##The code below will pull data from echoMRI.csv
#Read in echo_mri_data
echoMRI_data <- read_csv("../data/echomri.csv")

echoMRI_data_19 <- echoMRI_data %>%
  filter(COHORT == 19) %>%
  mutate(ID= as.factor(ID)) %>%
  select(ID, Fat, Lean, Weight, Date, adiposity_index, n_measurement) 

echoMRI_data_BWloss <- echoMRI_data_19 %>%
  filter(ID %in% c(3730, 3731,3732, 3733, 3735, 3736, 3737, 3738, 3739, 3740, 3741)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(DRUG = case_when(ID %in% c(3735, 3736, 3740) ~ "Vehicle",
                          ID %in% c(3730, 3731, 3738, 3741) ~ "Tirzepatide",
                          ID %in% c(3732, 3733, 3737, 3739) ~ "Survodutide"),
         DRUG = as.factor(DRUG),
         n_measurement = case_when(n_measurement == "1" ~ "Pre BW loss", # Day 0
                                   n_measurement == "2" ~ "During BW loss", # Day 21
                                   n_measurement == "3" ~ "End BW loss", # Day 34
                                   n_measurement == "4" ~ "Post BW regain")) %>%
  mutate(
    n_measurement = factor(n_measurement, 
                           levels = c("Pre BW loss", #pre-treatment
                                      "During BW loss", #During BW loss
                                      "End BW loss", #End BW loss
                                      "Post BW regain"))) #BW regain

## CORRECT I think - Calculate change in lean, fat, Weight, adiposity index ####
echmoMRI_Exp1_delta <- echoMRI_data_BWloss %>%
  ungroup() %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(FatplusLean = Fat + Lean) %>%
  mutate(delta_lean = Lean - first (Lean), # negative number means lean mass was lost
         delta_fat = Fat - first(Fat),
         delta_FatplusLean = FatplusLean - first(FatplusLean),
         delta_AI= adiposity_index - first(adiposity_index),
         Lean_FatplusLean = delta_lean/delta_FatplusLean,
         BW_pct_change = 100*((FatplusLean - first(FatplusLean)) / first(FatplusLean)),
         Lean_pct_change = 100*((Lean - first(Lean)) / first(Lean)),
         Lean_pct_BW = 100*(Lean/FatplusLean))

plot_echo_delta <-echmoMRI_Exp1_delta %>%
  filter_out(n_measurement=="Day 0") %>%
  ungroup() %>%
  group_by(DRUG, n_measurement) %>%
  summarise(Avg_delta_Lean_G = mean(delta_lean),
            SD_delta_Lean_G = sd(delta_lean),
            Avg_delta_Fat_G = mean(delta_fat),
            SD_delta_Fat_G = sd(delta_fat),
            Avg_delta_FatplusLean = mean(delta_FatplusLean),
            SD_delta_FatplusLean = sd(delta_FatplusLean),
            Avg_Lean_FatplusLean = mean(Lean_FatplusLean),
            SD_Lean_FatplusLean = sd(Lean_FatplusLean),
            Avg_Lean_pct_change = mean(Lean_pct_change),
            SD_Lean_pct_change = sd(Lean_pct_change),
            Avg_Lean_pct_BW = mean(Lean_pct_BW),
            SD_Lean_pct_BW = sd(Lean_pct_BW)) %>%
  select(n_measurement, Avg_delta_Lean_G, SD_delta_Lean_G, 
         Avg_delta_FatplusLean, SD_delta_FatplusLean, 
         Avg_Lean_FatplusLean, SD_Lean_FatplusLean, Avg_Lean_pct_change, SD_Lean_pct_change, 
         Avg_Lean_pct_BW, SD_Lean_pct_BW)

###Lean mass ####
#Graph with x axis as n_measurement and y axis as Lean
Lean_plot <-ggplot(echoMRI_data_BWloss, aes(x=n_measurement, y=Lean, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1, position = position_dodge(width = 0.3)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.15, position = position_dodge(width = 0.3)) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
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
  labs(x="Day",
       y= "Lean mass (g)",
       title= "Lean mass (g) during BW loss",
       color="Treatment", fill="Treatment")
Lean_plot

### Adiposity index ####
AI_plot <-ggplot(echoMRI_data_BWloss, aes(x=n_measurement, y=adiposity_index, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
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
  labs(x="Day",
       y= "Adiposity index (fat mass/lean mass)",
       title= "Adiposity index during BW loss",
       color="Treatment", fill="Treatment")
AI_plot

#### Percent change in Lean mass ####
Lean_pct_change_plot <-ggplot(echmoMRI_Exp1_delta, 
                              aes(x=n_measurement, y=Lean_pct_change, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
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
  labs(x="Day",
       y= "ΔLean mass (%)",
       title= "ΔLean mass (%) during BW loss & regain",
       color="Treatment", fill="Treatment")
Lean_pct_change_plot

#### Percent of BW comprised of lean mass ####
Lean_pct_of_BW_plot <-ggplot(echmoMRI_Exp1_delta, 
                              aes(x=n_measurement, y=Lean_pct_BW, group=DRUG, fill=DRUG, color=DRUG)) +
  geom_point(stat = "summary", 
             fun = "mean", aes(color=DRUG), size=4) +
  geom_line(stat = "summary", 
            fun = "mean", aes(color=DRUG), linewidth=1) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, aes(width=0.08), width=0.1) +
  #geom_jitter(width = 0.08, alpha = 0.6) +
  #geom_point(shape=1) +
  scale_color_manual(values = custom_colors_GLP) +
  #geom_point(aes(group=ID)) + geom_line(aes(group=ID, alpha=0.6)) +
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
  labs(x="Day",
       y= "Lean mass/Total body mass (%)",
       title= "Percent of BW comprised of lean mass",
       color="Treatment", fill="Treatment")
Lean_pct_of_BW_plot
