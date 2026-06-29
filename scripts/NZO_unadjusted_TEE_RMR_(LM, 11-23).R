#### Raw TEE (i.e. TEE not adjsted for BW or FFM) (use script from NZO_Fat_summary but with y=TEE) ####

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

#####----- Barplot - NZO raw TEE (measured values from echoMRI_data/sable_TEE_adj_RMR) ####
TEEraw_barplot_measured <- sable_TEE_adj_RMR %>%
  ggplot(aes(x = SABLE, y = tee, fill = GROUP)) + 
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
  labs(title = "NZO TEE (raw)", 
       y = "Raw TEE (kcal/day)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  #scale_y_continuous( # set y-axis breaks every 10 grams
  #breaks = seq(0, max(sable_TEE_adj_RMR$tee, na.rm = TRUE), by = 10)) +
  format.plot +
  theme(legend.position = "top",
        axis.ticks.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color="black", size=14),
        plot.title = element_text(hjust = 0.5, face = "bold"))
TEEraw_barplot_measured


####-----T-test, repeated measures ANOVA, post hoc for ANOVA -----####
#T-test: pairwise comparison between GROUP (ad lib and restricted) at each STATUS
ttest_results_rawTEE <- sable_TEE_adj_RMR %>%
  group_by(SABLE) %>%
  t_test(tee ~ GROUP, var.equal = TRUE) %>%   # or var.equal = FALSE if not assumed
  adjust_pvalue(method = "bonferroni") %>%   # optional multiple-comparison correction
  add_significance("p.adj")                  # adds stars based on adjusted p-values
ttest_results_rawTEE

#Repeated measures ANOVA: comparison between groups across the 5 time points
# Perform repeated-measures ANOVA separately for each GROUP
anova_results_rawTEE <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  group_by(GROUP) %>%
  anova_test(dv = tee, wid = ID, within = SABLE)
anova_results_rawTEE

#Post-hoc test for ANOVA
pairwise_results_rawTEE <- sable_TEE_adj_RMR %>%
  group_by(GROUP) %>%
  pairwise_t_test(tee ~ SABLE, paired = TRUE, p.adjust.method = "bonferroni") %>%
  filter(p.adj.signif <= 0.05)
pairwise_results_rawTEE

#### Raw RMR (i.e. RMR not adjsted for BW or FFM) (use script from NZO_Fat_summary but with y=RMR_kcal_day) ####

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

#####----- Barplot - NZO raw RMR (measured values from echoMRI_data/sable_TEE_adj_RMR) ####
rawRMR_barplot_measured <- sable_TEE_adj_RMR %>%
  ggplot(aes(x = SABLE, y = RMR_kcal_day, fill = GROUP)) + 
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
  labs(title = "NZO RMR (raw)", 
       y = "raw RMR (kcal/day)", 
       x= "Time point", 
       color = "Treatment group",
       fill = "Treatment group") +
  #scale_y_continuous( # set y-axis breaks every 10 grams
  #breaks = seq(0, max(sable_TEE_adj_RMR$RMR_kcal_day, na.rm = TRUE), by = 10)) +
  format.plot +
  theme(legend.position = "top",
        axis.ticks.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color="black", size=14),
        plot.title = element_text(hjust = 0.5, face = "bold"))
rawRMR_barplot_measured


####-----T-test, repeated measures ANOVA, post hoc for ANOVA -----####
#T-test: pairwise comparison between GROUP (ad lib and restricted) at each STATUS
ttest_results_rawRMR <- sable_TEE_adj_RMR %>%
  group_by(SABLE) %>%
  t_test(RMR_kcal_day ~ GROUP, var.equal = TRUE) %>%   # or var.equal = FALSE if not assumed
  adjust_pvalue(method = "bonferroni") %>%   # optional multiple-comparison correction
  add_significance("p.adj")                  # adds stars based on adjusted p-values
ttest_results_rawRMR

#Repeated measures ANOVA: comparison between groups across the 5 time points
# Perform repeated-measures ANOVA separately for each GROUP
anova_results_rawRMR <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  group_by(GROUP) %>%
  anova_test(dv = RMR_kcal_day, wid = ID, within = SABLE)
anova_results_rawRMR

#Post-hoc test for ANOVA
pairwise_results_rawRMR <- sable_TEE_adj_RMR %>%
  group_by(GROUP) %>%
  pairwise_t_test(RMR_kcal_day ~ SABLE, paired = TRUE, p.adjust.method = "bonferroni") %>%
  filter(p.adj.signif <= 0.05)
pairwise_results_rawRMR