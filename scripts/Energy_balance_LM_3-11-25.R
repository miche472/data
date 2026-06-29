
#Import data for 24hr energy expenditure (created in TEE_NZO.R and called 
  # sable_tee_data. Variable for 24 hr TEE is called "tee")

sable_tee_data <- sable_dwn %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% # just to see TEE in kcal first
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
  group_by(ID,complete_days,SABLE,is_complete_day) %>% 
  summarise(tee = sum(value)*(1/60)) %>% 
  filter(!ID %in% c(3715,3727, 3718, 3711, 3723), is_complete_day ==1) %>% #3715 died and the rest of the animals were measured in cage 5 
  ungroup() %>% 
  group_by(SABLE, ID) %>% 
  slice_max(order_by = complete_days, n=1)

#Import data for 24hr FI (created in TEE_NZO.R and called 
# sable_min_data. Variable for 24 hr food intake in kcal is called "kcal")

sable_min_data <- sable_dwn %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
  filter(grepl("FoodA_*", parameter)) %>% # just to see TEE in kcal first
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    value = cumsum(if_else(replace_na(value - lag(value),0)>0, 0, replace_na(value - lag(value),0))),
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  group_by(ID,complete_days,SABLE,is_complete_day) %>% 
  summarise(intake_gr = abs(max(value) - min(value))) %>% 
  filter(!ID %in% c(3715,3727, 3718, 3711, 3723), is_complete_day ==1) %>% #3715 died and the rest of the animals were measured in cage 5 
  ungroup() %>% 
  group_by(SABLE, ID) %>% 
  slice_max(order_by = complete_days, n=1) %>% 
  mutate(
    kcal = if_else(SABLE=="Before", intake_gr*3.1, intake_gr*3.82)
  )
-----------------------------
#Energy balance####
#Calculate energy balance for each mouse both before and after LFD
intake_kcal <- total_tee %>% 
  left_join(., total_intake, by = c("ID", "drug"))

March <- sable_tee_data %>% 
  ungroup() %>%
  left_join(., sable_min_data, by = "ID", "SABLE", "complete_days", "is_complete_day") %>%
  group_by(ID, SABLE) %>%
  mutate(balance = (tee - kcal))

March <- sable_tee_data %>% 
  ungroup() %>%
  left_join(., sable_min_data) %>%
  group_by(ID, SABLE) %>%
  mutate(balance = (tee - kcal))

March <- sable_tee_data %>% 
  ungroup() %>%
  left_join(., sable_min_data) %>%
  group_by(ID, SABLE) %>%
  mutate(balance = (kcal - tee))

#graph energy balance
March_plot  <- March %>% 
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, balance)) +
  
  # Bar for mean value
  stat_summary(
    fun = mean, 
    # geom = "col", 
    fill = "gray", 
    alpha = 0.5  # Transparency for the bar
  ) +
  
  # Individual points with transparency
  geom_point(aes(group = ID), alpha = 0.5) +
  
  # Lines connecting paired observations
  geom_line(aes(group = ID), alpha = 0.5) +
  
  # Mean with SEM as point and error bars
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.7,
    shape = 21,
    color = "black",
    fill = "red"
  ) +
  # Axis labels
  labs(x = NULL, y = "Energy Balance (kcal)") +
  
  # White background
  theme_classic()  
March_plot 

#Linear mixed model for energy balance
mdl_March <- lmer(
  data = March,
  balance ~ SABLE + (1|ID)
)
summary(mdl_March)

emmeans(
  mdl_March,
  pairwise ~ SABLE,
  type = "response"
)

#Is there a correlation between % change in TEE and % change in adiposity index?

#Energy flux####



