#### GRAPH: EE (kcal/hr) vs BW (Hrs 0-2) ####

#Try creating a df with the higher EE and the lower EE mice and look at EE vs doses
#Alternatively, could have created a new column for HIGH or LOW EE and then used facet_wrap

low_EE_group <- Summary_Combined_sum_loc_avg_EE_hr_bins %>%
  filter(ID %in% c(3731, 3732, 3733, 3735))

high_EE_group <- Summary_Combined_sum_loc_avg_EE_hr_bins %>%
  filter(ID %in% c(3737, 3738, 3739, 3740, 3741)) #removed 3740 due to outlier status

## Hours 0-2 post start of recording ####
### EE (kcal/hr) vs Dose (including baseline) ####
ggplot(low_EE_group, 
       aes(x = factor(DOSE), y = Avg_EE_0_24hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = low_EE_group %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "Low EE group: 0 - 24 hrs post")


## Hours 0-2 post start of recording ####
### EE (kcal/hr) vs Dose (including baseline) ####
ggplot(high_EE_group, 
       aes(x = factor(DOSE), y = Avg_EE_0_24hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = high_EE_group %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Energy expenditure (kcal/hr)", title = "High EE group: 0 - 24 hrs post")


## Hours 0-2 post start of recording ####
### Total distance (m) vs Dose (including baseline) ####
ggplot(low_EE_group, 
       aes(x = factor(DOSE), y = sum_loc_0_2hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = low_EE_group %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Total movement (m)", title = "Low EE group: 0 - 2 hrs post")

ggplot(high_EE_group, 
       aes(x = factor(DOSE), y = sum_loc_0_2hr, fill = factor(DOSE))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(DOSE)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = high_EE_group %>% group_by(ID) %>% slice_max(DOSE, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  scale_fill_manual(values = custom_colors_OXA) +
  scale_color_manual(values = custom_colors2_OXA) +
   theme(legend.position = "none") +
  labs(x = "Dose", y = "Total movement (m)", title = "High EE group: 0 - 2 hrs post")


#-------------------------------------------------#.
#-------------------------------------------------#.
#Cumulative sum of EE...does the time at which a mouse expends energy just shift when they are given orexin?
#They could have a compensatory mechanism that recognizes that brings down their EE in later hours to compensate for the 
#use df Combined_sum_EE_hr_bins. For each ID after each dose calculate the cumsum for sum_EE_60min

#delete
cumulative_EE <- Combined_sum_EE_hr_bins %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(recording_bin) %>%
  summarise(cum_sum_EE = cumsum(sum_EE_60min))
#end delete

cumulative_EE <- Combined_sum_EE_hr_bins %>%
  ungroup() %>%
  group_by(ID, DOSE) %>%
  arrange(recording_bin, .by_group = TRUE) %>%
  mutate(
    cum_sum_EE = cumsum(coalesce(sum_EE_60min, 0))
  ) %>%
  ungroup() %>%
  filter(recording_bin <5)
  

# Graphs ####
ggplot(
  cumulative_EE,
  aes(
    x = recording_bin,
    y = cum_sum_EE,
    color = DOSE,
    group = DOSE
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  scale_color_manual(values = custom_colors_OXA) +
  
  facet_wrap(~ID) +
  
  labs(
    x = "Recording Bin (60 min)",
    y = "Cumulative Energy Expenditure (kcal)",
    color = "Dose"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggplot(
  cumulative_EE,
  aes(
    x = recording_bin,
    y = cum_sum_EE,
    color = ID,
    group = ID
  )
) +
  geom_line(linewidth = 1) +
  
  scale_color_viridis_d() +
  
  facet_wrap(~PERIOD) +
  
  labs(
    x = "Recording Bin (60 min)",
    y = "Cumulative Energy Expenditure (kcal)",
    color = "Mouse ID"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

## Hours 0-2 post start of recording (kcal/hr) ####
###Avg EE (kcal/hr) vs Dose (including baseline) ####
plot_EE_0_2hr_avg <- ggplot(Summary_Combined_avg_EE_hr_bins, 
       aes(x = factor(PERIOD), y = Avg_EE_0_4hr, fill = factor(PERIOD))) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.2) +
  # Lines connecting the same mouse across doses
  geom_line(aes(group = ID),color = "gray50",linewidth = 0.7, alpha = 0.6) +
  geom_jitter(aes(color = factor(PERIOD)),width = 0.12,size = 2,alpha = 0.7) +
  #label lines with ID
  geom_text(data = Summary_Combined_avg_EE_hr_bins %>% group_by(ID) %>% slice_max(PERIOD, n = 1), aes(label = ID), hjust = -0.5, size = 3) +
  theme_bw(base_size = 14) +
  format.plot_LM3 +
  #scale_fill_manual(values = custom_colors_OXA) +
  #scale_color_manual(values = custom_colors2_OXA) +
  theme(legend.position = "none") +
  labs(x = "PERIOD", y = "Energy expenditure (kcal/hr)", title = "0 - 4 hrs post")
plot_EE_0_2hr_avg