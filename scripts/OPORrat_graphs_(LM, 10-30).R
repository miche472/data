#For Laurie graphs to show cathy in meeting on 10-30-25





reg<-lm(tee ~ Weight,
        data=sable_TEE_adj_RMR)                      

#get intercept and slope value
coeff<-coefficients(reg)          
intercept<-coeff[1]
slope<- coeff[2]

reg_line <- lm(tee ~ Weight, data=sable_TEE_adj_RMR)
summary(reg_line)
reg_line <- reg_line %>%
  mutate()

#### This works
reg_line <- lm(tee ~ Weight, data=sable_TEE_adj_RMR)
summary(reg_line)

my_coef <- coef(reg_line)

ggp <- ggplot(sable_TEE_adj_RMR, aes(x=Weight, y=tee, color=GROUP)) +   
  geom_point() +
  geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  lims(x = c(0, 60), y = c(0, 30)) #this allows you to choose the x and y axis ranges on graph
ggp

#### Graphs of TEE vs BW simple linear regression using geom_smooth() ####

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

####  Baseline: TEE vs BW ####
data_plot1 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Baseline")
ggplot(data_plot1, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Baseline: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction")

####  Peak Obesity: TEE vs BW ####
data_plot2 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Peak obesity")
ggplot(data_plot2, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peak Obesity: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction")

####  BW loss: TEE vs BW ####
data_plot3 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="BW loss")
ggplot(data_plot3, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "BW loss: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction")

####  BW regain: TEE vs BW ####
data_plot4 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="BW regain")
ggplot(data_plot4, aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  format.plot +
  scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "BW regain: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction")

####  Baseline (combined): TEE vs BW ####
data_plot5 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Baseline")
ggplot(data_plot5, aes(x = Weight, y = tee, color=GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2, color="black") +
  format.plot +
  #scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Baseline: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)")

####  Baseline (combined): TEE vs BW ####
n_distinct(data_plot5$ID) #good we have 16 animals
data_plot5 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Baseline")
ggplot(data_plot5, aes(x = Weight, y = tee, color=GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2, color="black") +
  format.plot +
  #scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Baseline: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)")

data5 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Baseline")

data5_lm <-lm(tee~Weight, data=data5)
summary(data5_lm)

####  Peak obesity (combined): TEE vs BW ####
data_plot6 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Peak obesity") %>%
  ungroup()
ggplot(data_plot6, aes(x = Weight, y = tee)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) + 
  #geom_abline(intercept = my_coef[1], slope = my_coef[2], size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2, color="black") +
  format.plot +
  #scale_color_manual(values=c('#FAAC41','#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peak obesity: TEE vs. Body weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)") +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))
n_distinct(data_plot6$ID) #good we have 16 animals

data6 <- sable_TEE_adj_RMR %>%
  filter(SABLE=="Peak obesity")

data6_lm <-lm(tee~Weight, data=data6)
summary(data6_lm)

#------------------------------------------------------------------------
####more efficient way to make graphs and equations ####
    #for BW loss, BW maintenance, and BW regain

library(tidyverse)
library(broom)

# Define the SABLE time points of interest
timepoints <- c("BW loss", "BW maintenance", "BW regain")

# Create a single pipeline:
model_results <- sable_TEE_adj_RMR %>%
  filter(SABLE %in% timepoints) %>%
  group_by(SABLE, GROUP) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(tee ~ Weight, data = .x)),
    tidy = map(model, tidy),
    glance = map(model, glance)
  ) %>%
  unnest(c(tidy, glance), names_sep = "_") %>%
  select(SABLE, GROUP, term = tidy_term, estimate = tidy_estimate, r.squared = glance_r.squared) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`, slope = Weight) %>%
  arrange(SABLE, GROUP)

#This data frame contains slope, intercept, and R² per SABLE × GROUP
model_results

#Make graphs for BW loss, BW maintenance, and BW regain 

sable_TEE_adj_RMR %>%
  filter(SABLE %in% timepoints) %>%
  left_join(model_results, by = c("SABLE", "GROUP")) %>%
  ggplot(aes(x = Weight, y = tee, color = GROUP)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  # 🔹 Plot the regression lines using the stored slopes & intercepts
  geom_abline(
    data = model_results,
    aes(intercept = intercept, slope = slope, color = GROUP),
    size = 1.2
  ) +
  scale_color_manual(values = c('#FAAC41', '#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TEE vs. Body Weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
    color = "Restriction"
  ) +
  lims(x = c(0, 60), y = c(0, 30)) +
  facet_wrap(~ SABLE)

#------------------------------------------------------------------------
####more efficient way to make graphs and equations ####
#for Baseline and peak obesity

# Define the SABLE time points of interest
timepoints2 <- c("Baseline", "Peak obesity")

# Create a single pipeline:
model_results2 <- sable_TEE_adj_RMR %>%
  filter(SABLE %in% timepoints2) %>%
  ungroup() %>%
  group_by(SABLE) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(tee ~ Weight, data = .x)),
    tidy = map(model, tidy),
    glance = map(model, glance)
  ) %>%
  unnest(c(tidy, glance), names_sep = "_") %>%
  select(SABLE, term = tidy_term, estimate = tidy_estimate, r.squared = glance_r.squared, p.value= glance_p.value) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`, slope = Weight) %>%
  arrange(SABLE)

#This data frame contains slope, intercept, and R² per SABLE
model_results2

#Make graphs for Baseline & Peak obesity 
sable_TEE_adj_RMR %>%
  filter(SABLE %in% timepoints2) %>%
  left_join(model_results2, by = c("SABLE")) %>%
  ggplot(aes(x = Weight, y = tee)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(
    aes(label = ID),
    hjust = -0.1, vjust = 0.5,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  # 🔹 Plot the regression lines using the stored slopes & intercepts
  geom_abline(
    data = model_results2,
    aes(intercept = intercept, slope = slope),
    size = 1.2
  ) +
  #scale_color_manual(values = c('#FAAC41', '#5392DB')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TEE vs. Body Weight",
    x = "Body Weight (g)",
    y = "TEE (kcal/day)",
  ) +
  lims(x = c(0, 60), y = c(0, 30)) +
  facet_wrap(~ SABLE)

