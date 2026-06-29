#Compare two methods for calculating NEAT: ####
ID_TEE3_hr_compare <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr) %>%
  summarise(
    # --- basic quantities ---
    n_obs = n(),
    total_minutes = n(),
    
    # --- total energy ---
    kcal_total = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # --- estimate RMR (from inactive periods) ---
    RMR_rate = mean(Kcal_Hr[move == 0], na.rm = TRUE),
    
    # --- Method 1: your current approach -> NEAT is the TEE when the mouse moves (overestimates NEAT) ---
    kcal_active_raw = sum((Kcal_Hr / 60)[move == 1], na.rm = TRUE),
    
    # --- Method 2: NEAT=TEE-RMR when the mouse is moving---
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE),
    
    # --- RMR total energy ---
    kcal_RMR = RMR_rate * (total_minutes / 60),
    
    # --- behavior ---
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(difference = kcal_active_raw - NEAT_kcal) %>%
  mutate(check = kcal_total - (kcal_RMR + NEAT_kcal))


#version for comparing two methods for calculating NEAT which also has a fallback to avoid NaN ####
ID_TEE3_hr_compare <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr) %>%
  summarise(
    n_obs = n(),
    total_minutes = n(),
    kcal_total = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # --- robust RMR ---
    RMR_rate = ifelse(
      sum(move == 0, na.rm = TRUE) > 0,
      mean(Kcal_Hr[move == 0], na.rm = TRUE),
      mean(Kcal_Hr, na.rm = TRUE)),
    
    # --- method 1 ---
    kcal_active_raw = sum((Kcal_Hr / 60)[move == 1], na.rm = TRUE),
    
    # --- method 2 (fixed) ---
    NEAT_kcal = sum(
      pmax((Kcal_Hr - RMR_rate) / 60, 0)[move == 1],
      na.rm = TRUE),
    kcal_RMR = RMR_rate * (total_minutes / 60),
    
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    difference = kcal_active_raw - NEAT_kcal,
    check = kcal_total - (kcal_RMR + NEAT_kcal)) 

#Version for checking coverage --> this tells me how many hour have less than 80% of the possible 60 observations
ID_TEE3_hr_compare <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  group_by(SABLE, ID, hr) %>%
  summarise(coverage = n() / 60) %>%
filter(!coverage >= 0.8)
#Total of 84 hours that have more than 20% of observations missing. 33 of these hours are from 
    #hour 18 (i.e. when recording was stopped for feeding, injections, and BW)
#Perhaps I should remove hr 18 from analyses? (of 80 total observations at hr 18, only 30 have 100% coverage)

#Version using percentile, 
#Calculate: TEE, RMR, and NEAT for IDs at each SABLE and during each hr of a 24hr day
#This version is directly from ChatGPT; below I have a version that is slightly modified
ID_TEE3_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  
  # Define movement
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  
  # Summarise per hour
  group_by(SABLE, ID, hr) %>%
  summarise(
    # -----------------------------
    # Basic info
    # -----------------------------
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    
    # -----------------------------
    # Total Energy Expenditure (TEE)
    # -----------------------------
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # -----------------------------
    # RMR (10th percentile method)
    # -----------------------------
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    
    # Total RMR energy across observed time
    RMR_kcal = RMR_rate * (n_obs / 60),
    
    # -----------------------------
    # NEAT (RMR-subtracted, clamped ≥ 0)
    # -----------------------------
    NEAT_kcal = sum(pmax((Kcal_Hr - RMR_rate) / 60, 0),
      na.rm = TRUE),
    
    # -----------------------------
    # Optional: raw activity EE (your original method)
    # -----------------------------
    activity_kcal_raw = sum((Kcal_Hr / 60)[move == 1], na.rm = TRUE),
    .groups = "drop") %>%
  
  # -----------------------------
# Diagnostics
# -----------------------------
mutate(
  # Difference between raw activity EE and true NEAT
  diff_activity_vs_NEAT = activity_kcal_raw - NEAT_kcal,
  
  # Check decomposition (should be ≥ 0 due to clamping)
  check = TEE_kcal - (RMR_kcal + NEAT_kcal))


#---------------------------------------------------------------------#
# My version using percentile for RMR --> happy with this...just check if global or light/dark RMR is better ####
ID_TEE3_hr <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  
  # Define movement
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  
  # Summarise per hour
  group_by(SABLE, ID, hr) %>%
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    
    # Total Energy Expenditure (TEE) --> minute by minute summation
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # Resting metabolic rate for each hour calculated using 10th percentile method for RMR:
    #Note: this is the rate of RMR rate within each hr not for the entire day
    #Calculating this will allow for calculation of NEAT
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    
    # Total RMR energy across observed time in the hour
    RMR_kcal = RMR_rate * (n_obs / 60),
    
    # NEAT= TEE-RMR when the mouse is moving
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE), .groups = "drop") %>%
  
    # Verify that TEE = RMR + NEAT --> TEEvsRMR_NEAT should be close to zero
    mutate(TEEvsRMR_NEAT = TEE_kcal - (RMR_kcal + NEAT_kcal))

#Copy of code to check if calculating global RMR or hourly RMR is better (compare the two options)
# First, calculate the global (SABLE × ID) percentile RMR
RMR_global <- filter_locom_energy3 %>%
  group_by(SABLE, ID) %>%
  summarise(
    RMR_rate_global = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

# Now summarise per hour using this global RMR
ID_TEE3_hr_globalRMR <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  
  # Define movement
  group_by(SABLE, ID, GROUP) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  
  # Summarise per hour
  group_by(SABLE, ID, hr) %>%
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # Join in the global RMR
    RMR_rate_global = RMR_global$RMR_rate_global[RMR_global$SABLE == unique(SABLE) &
                                                   RMR_global$ID == unique(ID)],
    
    # Total RMR energy using the global rate
    RMR_kcal_global = RMR_rate_global * (n_obs / 60),
    
    # NEAT using your preferred method
    NEAT_kcal_global = sum(((Kcal_Hr - RMR_rate_global) / 60)[move == 1], na.rm = TRUE),
    
    # Optional: raw activity EE
    activity_kcal_raw = sum((Kcal_Hr / 60)[move == 1], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  # Diagnostics
  mutate(
    diff_activity_vs_NEAT = activity_kcal_raw - NEAT_kcal_global,
    check_global = TEE_kcal - (RMR_kcal_global + NEAT_kcal_global)
  )

#Do the comparison of Global vs hourly 
comparison <- ID_TEE3_hr %>%
  left_join(ID_TEE3_hr_globalRMR %>%
              select()
            by = c("SABLE", "ID", "hr")) %>%
  mutate(
    NEAT_diff = NEAT_kcal - NEAT_kcal_global,
    RMR_diff  = RMR_kcal - RMR_kcal_global
  )

#Create a plot to visualize differences between hourly and global

# Combine hourly and global NEAT for plotting
NEAT_plot_data <- ID_TEE3_hr %>%
  select(SABLE, ID, hr, NEAT_kcal) %>%
  rename(NEAT_hourly = NEAT_kcal) %>%
  left_join(
    ID_TEE3_hr_globalRMR %>%
      select(SABLE, ID, hr, NEAT_kcal_global) %>%
      rename(NEAT_global = NEAT_kcal_global),
    by = c("SABLE", "ID", "hr")
  ) %>%
  pivot_longer(
    cols = c(NEAT_hourly, NEAT_global),
    names_to = "NEAT_type",
    values_to = "NEAT_value"
  )

# Plot
ggplot(NEAT_plot_data, aes(x = hr, y = NEAT_value, color = NEAT_type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ID, scales = "free_y") +
  scale_color_manual(values = c("NEAT_hourly" = "blue", "NEAT_global" = "red"),
                     labels = c("Hourly percentile RMR", "Global percentile RMR")) +
  labs(
    title = "Comparison of NEAT calculated using Hourly vs Global Percentile RMR",
    x = "Hour",
    y = "NEAT (kcal)",
    color = "NEAT Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

#Calculate RMR by light vs dark cycle (rather than hourly or globally)
ID_TEE_lightdark <- filter_locom_energy3 %>%
  ungroup() %>%
  arrange(DateTime) %>%
  
  # Define movement
  group_by(SABLE, ID, GROUP, lights) %>%
  mutate(move = if_else(All_meters > lag(All_meters), 1, 0, missing = 0)) %>%
  ungroup() %>%
  
  # Summarise per hour within light/dark cycle
  group_by(SABLE, ID, lights, hr) %>%
  summarise(
    n_obs = n(),
    minutes_active = sum(move == 1, na.rm = TRUE),
    minutes_rest   = sum(move == 0, na.rm = TRUE),
    
    # Total energy expenditure (minute-by-minute summation)
    TEE_kcal = sum(Kcal_Hr / 60, na.rm = TRUE),
    
    # RMR for light/dark cycle using 10th percentile
    # RMR_rate is the 10th percentile within this light/dark cycle and SABLE × ID
    RMR_rate = quantile(Kcal_Hr, probs = 0.10, na.rm = TRUE, names = FALSE),
    
    # Total RMR energy in this hour
    RMR_kcal = RMR_rate * (n_obs / 60),
    
    # NEAT: subtract RMR only during movement (your preferred method)
    NEAT_kcal = sum(((Kcal_Hr - RMR_rate) / 60)[move == 1], na.rm = TRUE),
    
    # Optional: raw activity EE for comparison
    activity_kcal_raw = sum((Kcal_Hr / 60)[move == 1], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  # Diagnostics
  mutate(
    diff_activity_vs_NEAT = activity_kcal_raw - NEAT_kcal,
    check = TEE_kcal - (RMR_kcal + NEAT_kcal)
  )

#Plot comparing hourly, light/dark, and global RMR calculation method
# Prepare the three NEAT datasets
NEAT_compare <- ID_TEE3_hr %>%
  select(SABLE, ID, hr, NEAT_kcal) %>%
  rename(NEAT_hourly = NEAT_kcal) %>%
  left_join(
    ID_TEE_lightdark %>%
      select(SABLE, ID, hr, lights, NEAT_kcal) %>%
      rename(NEAT_lightdark = NEAT_kcal),
    by = c("SABLE", "ID", "hr")
  ) %>%
  left_join(
    ID_TEE3_hr_globalRMR %>%
      select(SABLE, ID, hr, NEAT_kcal_global) %>%
      rename(NEAT_global = NEAT_kcal_global),
    by = c("SABLE", "ID", "hr")
  )

# Pivot longer for plotting
NEAT_plot_long <- NEAT_compare %>%
  pivot_longer(
    cols = c(NEAT_hourly, NEAT_lightdark, NEAT_global),
    names_to = "NEAT_method",
    values_to = "NEAT_value"
  )
ggplot(NEAT_plot_long, aes(x = hr, y = NEAT_value, color = NEAT_method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ID, scales = "free_y") +
  scale_color_manual(
    values = c(
      "NEAT_hourly" = "blue",
      "NEAT_lightdark" = "green",
      "NEAT_global" = "red"
    ),
    labels = c(
      "Hourly percentile RMR",
      "Light/Dark percentile RMR",
      "Global percentile RMR"
    )
  ) +
  labs(
    title = "Comparison of NEAT calculated using different RMR methods",
    x = "Hour",
    y = "NEAT (kcal)",
    color = "NEAT Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )