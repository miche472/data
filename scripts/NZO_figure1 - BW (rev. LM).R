# This script aims to explore changes in body weight in middle age NZO after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)

#####Organize data for BW analysis####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    ),
    day_rel = DATE - first(DATE)
  ) %>%
  mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance", 
      COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      day_rel == 161 ~ "BW loss",
      DATE == as.Date("2025-02-24") ~ "peak obesity",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(STATUS)) %>% 
  filter(!(ID == 3726 & DATE == as.Date("2025-04-28")))  #repeated


# Make STATUS an ordered factor
BW_data <- BW_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))
BW_data_csv %>%
  select(ID, BW, DATE, )
write.csv()

#format plot

####Plot 1: NZO BW, collapsed drug groups####

#LM:I use a mix of "color" and "fill" in this code chunk since it was originally 
#written to show both diet group and drug group. Now I want the bars to show mean values
#for diet groups and the points to show individual values for diet groups. Therefore I 
#manually "synced" the color and fill using the code directly below this note
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))

format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  # remove background grid lines only
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # keep axis lines
  axis.line = element_line(color = "black"))


#### bar plot of Measured BW values ####
plot1 <- BW_data %>%
  ggplot(aes(x = STATUS, y = BW, fill = GROUP)) + 
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
    aes(color = GROUP), ##LM changed color= DRUG to color = GROUP
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) + 
  scaleFill + scaleColor +
  theme_minimal() +
  labs(title = "NZO body weight (drug group collapsed)", 
       y = "Body Weight (grams)", 
       x= "Time point", 
       color = "Diet group",
       fill = "Diet group") +
  # set y-axis breaks every 10 grams
  scale_y_continuous(breaks = seq(0, max(BW_data$BW, na.rm = TRUE), by = 10)) +
  format.plot +
  theme(
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))
plot1

####Stats Plot 1: NZO BW, collapsed drug group####
# Fit model
model1 <- lmer(BW ~ STATUS * GROUP + (1|ID), data = BW_data)
summary(model1)

# Save emmeans results
emmeans_model1_results <- emmeans(model1, pairwise ~ STATUS * GROUP, adjust = "tukey")

#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm_model1 <- as.data.frame(emmeans_model1_results$emmeans)   # estimated means
df_pairs_model1 <- as.data.frame(emmeans_model1_results$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm_model1, n = Inf)
print(df_pairs_model1, n = Inf)
# Keep only significant contrasts
df_sig_model1 <- df_pairs_model1 %>%
  filter(p.value <= 0.05)

# View the results
print(df_sig_model1, n = Inf)
View(df_sig_model1)


#### Use this to check accuracy of graph --> 
#See labeled points on drug collapsed graph (LM) #
#This lets me see the labeled points (the portion of code for this is from ChatGPT)
library(ggplot2)
library(dplyr)

check_plot1 <- BW_data %>%
  ggplot(aes(x = STATUS, y = BW, fill = GROUP)) + 
  # mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  # error bars (mean ± SE)
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  # individual data points
  geom_point(
    aes(color = GROUP),
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) +
  # labels for each point (y-value)
  geom_text_repel(
    aes(label = round(BW, 1), color = GROUP),
    position = position_dodge(width = 0.8),
    size = 3,
    show.legend = FALSE) +
  scaleFill +
  theme_minimal() +
  labs(y = "Body Weight (grams)", color = "GROUP") +
  format.plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))
check_plot1
#-----------

#### Plot 2: NZO BW (diet & drug groups)####
BW_plotdata <- BW_data %>%
  mutate(
    PlotGroup = case_when(
      STATUS %in% c("baseline", "peak obesity") ~ "all",          # collapse all
      STATUS %in% c("BW loss", "BW maintenance") ~ GROUP,         # separate by GROUP
      STATUS == "BW regain" ~ paste(GROUP, DRUG, sep = "_")       # GROUP × DRUG
    )
  )

# Define custom colors (updated by LM)
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#BB6509",              # orange
  "restricted" = "#246997",          # sky blue
  "ad lib_vehicle" = "#BB6509",      # darker orange
  "ad lib_RTIOXA_47" = "#F39C12",    # lighter orange
  "restricted_vehicle" = "#246997",  # darker blue
  "restricted_RTIOXA_47" = "#8CD3FF" # lighter blue
)

# Define custom colors (original from CS)
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#E67E22",              # orange
  "restricted" = "#3498DB",          # sky blue
  "ad lib_vehicle" = "#E67E22",      # darker orange
  "ad lib_RTIOXA_47" = "#F39C12",    # lighter orange
  "restricted_vehicle" = "#3498DB",  # darker blue
  "restricted_RTIOXA_47" = "#5DADE2" # lighter blue
)

plot2 <- BW_plotdata %>%
  ggplot(aes(x = STATUS, y = BW, fill = PlotGroup)) +
  
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7, alpha = 0.7) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  
  geom_point(aes(color = PlotGroup),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  theme_minimal() +
  labs(title= "NZO body weight (drug & diet)", 
       x = "Time point", 
       y = "Body Weight (g)", 
       fill = "Treatment group", 
       color = "Treatment group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  format.plot

plot2

####Stats Plot2: NZO BW, diet & drug groups####
# Fit model
model <- lmer(BW ~ STATUS * GROUP * DRUG + (1|ID), data = BW_data)
summary(model)

# Save emmeans results
emmeans_results <- emmeans(model, pairwise ~ STATUS * GROUP, adjust = "tukey")

#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm <- as.data.frame(emmeans_results$emmeans)   # estimated means
df_pairs <- as.data.frame(emmeans_results$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm, n = Inf)
print(df_pairs, n = Inf)
# Keep only significant contrasts
df_sig <- df_pairs %>%
  filter(p.value <= 0.05)

# View the results
print(df_sig, n = Inf)
View(df_sig)

#### Plot 3: NZO BW, vehicle ####
  
scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  
  # remove background grid lines only
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  # keep axis lines
  axis.line = element_line(color = "black")
)

#Select only mice that received vehicle (filtered out agonist mice)
plot3_data <- BW_data %>%
  filter(DRUG == "vehicle")
plot3 <- plot3_data %>%
  ggplot(aes(x = STATUS, y = BW, fill = GROUP)) + 
  # mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7) +
  # error bars (mean ± SE)
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3) +
  # individual data points
  geom_point(
    aes(color = GROUP), ##LM changed color= DRUG to color = GROUP
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2) +
  scaleFill + scaleColor +
  theme_minimal() +
  labs(title = "NZO body weight (vehicle only)", 
       y = "Body Weight (grams)", 
       x= "Time point", 
       color = "Restriction group",
       fill = "Restriction group") +
  format.plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))
plot3

####Stats Plot3: NZO BW, vehicle ####
model3 <- lmer(BW ~ STATUS * GROUP + (1|ID), data = plot3_data)
summary(model3)

# Save emmeans results
emmeans_results3 <- emmeans(model3, pairwise ~ STATUS * GROUP, adjust = "tukey")

#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm3 <- as.data.frame(emmeans_results3$emmeans)   # estimated means
df_pairs3 <- as.data.frame(emmeans_results3$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm3, n = Inf)
print(df_pairs3, n = Inf)
# Keep only significant contrasts
df_sig3 <- df_pairs3 %>%
  filter(p.value <= 0.05)

# View the results
print(df_sig3, n = Inf)
View(df_sig3)

#baseline: ad lib and restricted are ns different
#peak obesity: ad lib and restricted are ns different
#BW loss: ad lib and restricted, p=3.227269e-04
#BW maintenance: ad lib and restricted, p= 9.342443e-04
#BW regain: ad lib and restricted, p=4.828269e-02

#BW maintenance and BW regain for restricted: estimate= 4.706667 , p=2.166625e-02

#### Plot 4: NZO BW, RTIOXA-47 ####

scaleFill <- scale_fill_manual(values = c("#FAAC41", "#5392DB"))
scaleColor <- scale_color_manual(values = c("#C77314", "#183873"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  
  # remove background grid lines only
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  # keep axis lines
  axis.line = element_line(color = "black")
)

#Select only mice that received vehicle (no mice with agonist)
plot4_data <- BW_data %>%
  filter(DRUG == "RTIOXA_47")

plot4 <- plot4_data %>%
  ggplot(aes(x = STATUS, y = BW, fill = GROUP)) + 
  
  # mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7
  ) +
  # error bars (mean ± SE)
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3
  ) +
  
  # individual data points
  geom_point(
    aes(color = GROUP), #
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2
  ) + 
  
  scaleFill + scaleColor +
  theme_minimal() +
  labs(title = "NZO body weight (RTIOXA-47 only)", 
       y = "Body Weight (grams)", 
       x= "Time point", 
       color = "Restriction group",
       fill = "Restriction group") +
  format.plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

plot4

####Stats Plot 4: NZO BW, RTIOXA47 ####
model4 <- lmer(BW ~ STATUS * GROUP + (1|ID), data = plot4_data)
summary(model4)

# Save emmeans results
emmeans_results4 <- emmeans(model4, pairwise ~ STATUS * GROUP, adjust = "tukey")

#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm4 <- as.data.frame(emmeans_results4$emmeans)   # estimated means
df_pairs4 <- as.data.frame(emmeans_results4$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm4, n = Inf)
print(df_pairs4, n = Inf)
# Keep only significant contrasts
df_sig4 <- df_pairs4 %>%
  filter(p.value <= 0.05)

# View the results
print(df_sig4, n = Inf)
View(df_sig4)

#baseline: ad lib and restricted, p=
#peak obesity: ad lib and restricted, p=
#BW loss: ad lib and restricted, p=
#BW maintenance: ad lib and restricted, p=
#BW regain: ad lib and restricted, p=

#BW maintenance and BW regain for restricted: 
