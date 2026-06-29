#Graphs and simple linear regression of EE (unadjusted) vs BW at each time point

#1.Five graphs: at each sable time point-->TEE vs BW 
  #simple linear regression for restricted mice and for ad libitum mice
  #y-intercept, slope, R^2, p value for each regression line

#2.Five graphs: at each sable time point-->TEE vs lean mass 
  #simple linear regression for restricted mice and for ad libitum mice
  #y-intercept, slope, R^2, p value for each regression line

#3.Five graphs: at each sable time point-->REE vs BW 
  #simple linear regression for restricted mice and for ad libitum mice
  #y-intercept, slope, R^2, p value for each regression line

#4.Five graphs: at each sable time point-->REE vs lean mass 
  #simple linear regression for restricted mice and for ad libitum mice
  #y-intercept, slope, R^2, p value for each regression line
  
#Libraries####
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
library(broom)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

#Load in sable data
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#Run script called: NZO_Figure7b-RMR_correctedbyLean (LM, 10-17).R
  #generates a df containing RMR, TEE, Lean, Weight for NZO
  #data frame is named: sable_TEE_adj_RMR 
#Use sable_TEE_adj_RMR for regression and graphs

#Write into a csv (not necessary for this script)
write_csv(x = sable_TEE_adj_RMR, "../data/sable_TEE_adj_RMR.csv")

#### 1. TEE vs BW (baseline, peak obesity, BW loss, BW maintenance, BW regain), NZO ####

TEEvsBW <- sable_TEE_adj_RMR %>%
  ungroup()

# Fit separate models for each SABLE × GROUP combination
model_stats_TEEvsBW <- TEEvsBW %>%
  group_by(SABLE, GROUP) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(tee ~ Weight, data = .x)),
    glance_stats = map(model, glance),
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    r_squared = map_dbl(glance_stats, "r.squared"),
    p_value = map_dbl(glance_stats, "p.value")) %>%
  select(SABLE, GROUP, intercept, slope, r_squared, p_value)

# Add label text for each regression
labels <- model_stats_TEEvsBW %>%
  mutate(
    label = sprintf(
      "y = %.2f + %.2f·x\nR² = %.3f\np = %.3f",
      intercept, slope, r_squared, p_value))

# Compute placement coordinates (inside each facet, near the top-right)
label_positions <- TEEvsBW %>%
  group_by(SABLE, GROUP) %>%
  summarise(
    x_pos = max(Weight, na.rm = TRUE) * 0.85,
    y_pos = max(tee, na.rm = TRUE) * 0.9,
    .groups = "drop") %>%
  left_join(labels, by = c("SABLE", "GROUP"))

#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

# Plot all time points with regression lines + formatted text
ggplot(TEEvsBW, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_text(
    data = label_positions,
    aes(x = x_pos, y = y_pos, label = label),
    #color = "black",            # text color
    size = 4.5,                 # text size (larger)
    hjust = 1, vjust = 1,
    lineheight = 1.1,
    inherit.aes = FALSE
  ) + format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction group") + 
    format.plot +
  facet_wrap(~ SABLE, scales="free") 

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#### 2. TEE vs lean (baseline, peak obesity, BW loss, BW maintenance, BW regain), NZO ####

TEEvsLean <- sable_TEE_adj_RMR %>%
  ungroup()

# Fit separate models for each SABLE × GROUP combination
model_stats_TEEvsLean <- TEEvsLean %>%
  group_by(SABLE, GROUP) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(tee ~ Lean, data = .x)),
    glance_stats = map(model, glance),
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    r_squared = map_dbl(glance_stats, "r.squared"),
    p_value = map_dbl(glance_stats, "p.value")) %>%
  select(SABLE, GROUP, intercept, slope, r_squared, p_value)

# Add label text for each regression
labels <- model_stats_TEEvsLean %>%
  mutate(
    label = sprintf(
      "y = %.2f + %.2f·x\nR² = %.3f\np = %.3f",
      intercept, slope, r_squared, p_value))

# Compute placement coordinates (inside each facet, near the top-right)
label_positions <- TEEvsLean %>%
  group_by(SABLE, GROUP) %>%
  summarise(
    x_pos = max(Lean, na.rm = TRUE) * 0.85,
    y_pos = max(tee, na.rm = TRUE) * 0.9,
    .groups = "drop") %>%
  left_join(labels, by = c("SABLE", "GROUP"))

#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

# Plot all time points with regression lines + formatted text
ggplot(TEEvsLean, aes(x = Lean, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_text(
    data = label_positions,
    aes(x = x_pos, y = y_pos, label = label),
    #color = "black",            # text color
    size = 4.5,                 # text size (larger)
    hjust = 1, vjust = 1,
    lineheight = 1.1,
    inherit.aes = FALSE
  ) + format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  facet_wrap(~ SABLE, scales="free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "TEE vs. Lean mass",
    x = "Lean mass (g)",
    y = "TEE (kcal/day)",
    color = "Restriction group") + format.plot


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#### 3. REE vs BW (baseline, peak obesity, BW loss, BW maintenance, BW regain), NZO ####

REEvsBW <- sable_TEE_adj_RMR %>%
  ungroup()

# Fit separate models for each SABLE × GROUP combination
model_stats_REEvsBW <- REEvsBW %>%
  group_by(SABLE, GROUP) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(RMR_kcal_day ~ Weight, data = .x)),
    glance_stats = map(model, glance),
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    r_squared = map_dbl(glance_stats, "r.squared"),
    p_value = map_dbl(glance_stats, "p.value")) %>%
  select(SABLE, GROUP, intercept, slope, r_squared, p_value)

# Add label text for each regression
labels <- model_stats_REEvsBW %>%
  mutate(
    label = sprintf(
      "y = %.2f + %.2f·x\nR² = %.3f\np = %.3f",
      intercept, slope, r_squared, p_value))

# Compute placement coordinates (inside each facet, near the top-right)
label_positions <- REEvsBW %>%
  group_by(SABLE, GROUP) %>%
  summarise(
    x_pos = max(Weight, na.rm = TRUE) * 0.85,
    y_pos = max(RMR_kcal_day, na.rm = TRUE) * 0.9,
    .groups = "drop") %>%
  left_join(labels, by = c("SABLE", "GROUP"))

#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

# Plot all time points with regression lines + formatted text
ggplot(REEvsBW, aes(x = Weight, y = RMR_kcal_day, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_text(
    data = label_positions,
    aes(x = x_pos, y = y_pos, label = label),
    #color = "black",            # text color
    size = 4.5,                 # text size (larger)
    hjust = 1, vjust = 1,
    lineheight = 1.1,
    inherit.aes = FALSE
  ) + format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  facet_wrap(~ SABLE, scales="free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "REE vs. Body weight",
    x = "Body weight (g)",
    y = "REE (kcal/day)",
    color = "Restriction group"
  ) + format.plot

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#### 4. REE vs Lean (baseline, peak obesity, BW loss, BW maintenance, BW regain), NZO ####

REEvsLean <- sable_TEE_adj_RMR %>%
  ungroup()

# Fit separate models for each SABLE × GROUP combination
model_stats_REEvsLean <- REEvsLean %>%
  group_by(SABLE, GROUP) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(RMR_kcal_day ~ Lean, data = .x)),
    glance_stats = map(model, glance),
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    r_squared = map_dbl(glance_stats, "r.squared"),
    p_value = map_dbl(glance_stats, "p.value")) %>%
  select(SABLE, GROUP, intercept, slope, r_squared, p_value)

# Add label text for each regression
labels <- model_stats_REEvsLean %>%
  mutate(
    label = sprintf(
      "y = %.2f + %.2f·x\nR² = %.3f\np = %.3f",
      intercept, slope, r_squared, p_value))

# Compute placement coordinates (inside each facet, near the top-right)
label_positions <- REEvsLean %>%
  group_by(SABLE, GROUP) %>%
  summarise(
    x_pos = max(Lean, na.rm = TRUE) * 0.85,
    y_pos = max(RMR_kcal_day, na.rm = TRUE) * 0.9,
    .groups = "drop") %>%
  left_join(labels, by = c("SABLE", "GROUP"))

#Format plot (optional)
format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  panel.grid.major = element_blank(), # remove background grid lines only
  panel.grid.minor = element_blank(), # remove background grid lines only
  axis.line = element_line(color = "black")) # keep axis lines

# Plot all time points with regression lines + formatted text
ggplot(REEvsLean, aes(x = Lean, y = RMR_kcal_day, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_text(
    data = label_positions,
    aes(x = x_pos, y = y_pos, label = label),
    #color = "black",            # text color
    size = 4.5,                 # text size (larger)
    hjust = 1, vjust = 1,
    lineheight = 1.1,
    inherit.aes = FALSE
  ) + format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  facet_wrap(~ SABLE, scales="free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "REE vs. Lean mass",
    x = "Lean mass (g)",
    y = "REE (kcal/day)",
    color = "Restriction group"
  ) + format.plot

####---------------------------------------------------------------------####
####---------------------------------------------------------------------####
#Reorganization of the figure
library(tidyverse)
library(broom)
library(ggrepel)  # optional but recommended if you use mouse ID labels

# Fit regressions for each SABLE × GROUP
model_stats <- TEEvsBW %>%
  group_by(SABLE, GROUP) %>%
  do({
    fit <- lm(tee ~ Weight, data = .)
    glance_fit <- glance(fit)
    tibble(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      r_squared = glance_fit$r.squared,
      p_value = glance_fit$p.value
    )
  }) %>%
  ungroup()

# Create label text for each regression
labels <- model_stats %>%
  mutate(
    label = sprintf("y = %.2f + %.2f·x\nR² = %.3f\np = %.3f",
                    intercept, slope, r_squared, p_value)
  )

# Compute label positions (left and near top inside each facet)
label_positions <- TEEvsBW %>%
  group_by(SABLE, GROUP) %>%
  summarise(
    x_pos = min(Weight, na.rm = TRUE) * 1.05,  # shift slightly right of min x
    y_pos = max(tee, na.rm = TRUE) * 0.95,     # slightly below max y
    .groups = "drop"
  ) %>%
  left_join(labels, by = c("SABLE", "GROUP"))

# Plot
ggplot(TEEvsBW, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_text(
    data = label_positions,
    aes(x = x_pos, y = y_pos, label = label),
    color = "black",
    size = 4.8,          # large text
    hjust = 0, vjust = 1,
    lineheight = 1.1,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ SABLE, ncol = 1, scales = "free") +  # stack vertically
  scale_color_manual(values = c("ad libitum" = "darkgray",
                                "restricted" = "darkblue")) +
  labs(
    title = "TEE vs. BW",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
