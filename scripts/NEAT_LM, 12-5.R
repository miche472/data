#The goal here is to calculate NEAT so that I can add it to the calculated RMR and see how close 
    #it is to the calculated TEE. This is following up from the meeting with Jen Teske

#Trying to directly calculate NEAT by adding up the energy expenditure during minutes 
    #during which the mouse moved

#------- Get only the kcal_hr parameter --> value column has kcal_hr data only
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

filter_EE <-sable_dwn %>%
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
  filter(grepl("kcal_hr_*", parameter))


#------- Get only the AllMeters parameter --> value column has AllMeters data
filter_loc <-sable_dwn %>%
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
  filter(grepl("AllMeters_*", parameter))

#---
#In df filter_loc use mutate to make a column called AllMeters_ using the value column data
#In df filter_EE use mutate to make a column called kcal_hr_ using the value column data

filter_loc_ <- filter_loc %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(All_meters = value) %>%
  rename(parameter_AllMeters = parameter) %>%
  rename(fix_value_AllMeters = fix_value) %>%
  ungroup() %>%
  group_by(DateTime, ID)

filter_EE_ <- filter_EE %>%
  ungroup() %>%
  group_by(DateTime, ID) %>%
  rename(Kcal_Hr = value) %>%
  rename(parameter_kcal_hr = parameter) %>%
  rename(fix_value_kcal_hr = fix_value)
  ungroup() %>%
  group_by(DateTime, ID)
  

#Join filter_loc and filter_EE into a df called filter_loc_EE

# Left join Lean, Fat, and Weight info into TEE dataset
filter_loc_EE <- filter_loc_ %>%
  left_join(
    filter_EE_ %>% 
      select(Kcal_Hr, ID, DateTime),
    by = c("ID", "DateTime"))
####**** this is giving a warning so i need to fix it ####
# Next step is to do if then statement...all meters is cumulative so if the lag (All_meters) < All_meters
  # then a new column called move is assigned a value of 1 and if not then 0. 
  # the idea is to identify if the mouse moved during that minute of measurement
  # if move=1 for a given ID and DateTime then the mouse moved, so the Kcal_Hr for that mouse/minute
  #should be added to a column called NEAT (this is a way of "directly" calculating non RMR energy expenditure)



