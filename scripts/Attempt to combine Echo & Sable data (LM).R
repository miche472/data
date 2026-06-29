
#LM attempt on 3-14-2025

#Obective####
  #Use the 24hr TEE taken Before and After LFD
  #Must include lean mass as a covariate. Need to have lean mass at the time of Before and After
  #Which is going to be different for each animal


#Create TEE data frame #### (this is taken from the script called TEE_NZO.R 
  #The only thing I changed was adding back in 3723, 3711, 3718, 3727). Also I ungrouped the data in the last line of this block
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
  filter(!ID %in% c(3715), is_complete_day ==1) %>% #3715 died 
  ungroup() %>% 
  group_by(SABLE, ID) %>% 
  slice_max(order_by = complete_days, n=1)
  (ungroup)

#Add a column called SABLE to the EchoMRI data so that it can be joined with the Sable data#### 
  #This code doesn't work yet, but I'm not sure why not
library(dplyr)
energy_echo <- echomri_data  %>% 
    ungroup() %>%
    mutate(SABLE = case_when(
    between(Date, as.Date("2024-11-11"), as.Date("2024-12-14"))~ "Before"),
    between(Date, as.Date("2025-01-26"), as.Date("2025-2-07"))~"After") %>%
  drop_na(SABLE) %>%
    left_join(., sable_tee_data)
  group_by(ID, SABLE)
 
#Create a mixed multiple linear regression model for TEE with SABLE, lean, and ID as explanatory variables####
mdl_energy_echo <- lmer(
    data = energy_echo,
    tee ~ SABLE + lean + (1|ID)
  )
  summary(mdl_energy_echo)
  
  emmeans(
    mdl_energy_echo,
    pairwise ~ SABLE,
    type = "response"
  )


#--------------------------------------------------------
#Various attempts
#Attempt 1-- to add SABLE as a column in EchoMRI data and then combine with Sable data
    energy_echo <- echomri_data  %>% 
      mutate(SABLE = if_else(Date >as.Date(2025-01-26) & as.Date(Date<2025-02-08), "After", "NA", if_else(Date>as.Date("2024-11-11") & Date<as.Date("2024-12-15"), "Before", "NA")))

#Attempt 2          
  Energy_echo <- echomri_data  %>% 
               mutate(SABLE =if_else(Date>2025-01-26 & Date<2025-02-08, "After", "NA"))
             
#Attempt 3 (Attempted to use method recommended on stack over flow)
             ifelse (Date >= 2025- | d < SE, "Winter",
                     ifelse (d >= SE & d < SS, "Spring",
                             ifelse (d >= SS & d < FE, "Summer", "Fall")))
             
             #This link might have information that could work
             #https://stackoverflow.com/questions/61199314/if-else-statements-for-dates-in-r
             
             
             getSeason <- function(d) {
               WS <- as.Date("2024-11-11") # Date of first pre-sable "Before" EchoMRI
               SE <- as.Date("2024-12-14") # Date of final pre-sable "Before" EchoMRI
               SS <- as.Date("2025-1-26") # Date of first pre-sable "After" EchoMRI
               FE <- as.Date("2025-2-8") # Date of final pre-sable "After" EchoMRI
               
               ifelse (d >= WS & d < SE, "Before",
                       ifelse (d >= SS & d <= FE, "After"))
             }
             #Change to standard date format
             SampleData$date <- as.Date(SampleData$SightDate, format = '%m/%d/%Y')
             #Make date of the same year i.e 2016
             SampleData$date <- as.Date(format(SampleData$date, "2016-%m-%d"))
             #Get season for each date. 
             SampleData$SightSeason <- getSeason(SampleData$date)

             getSeason <- function(d) {
               WS <- as.Date("2024-11-11") # Date of first pre-sable "Before" EchoMRI
               SE <- as.Date("2024-12-14") # Date of final pre-sable "Before" EchoMRI
               SS <- as.Date("2025-1-26") # Date of first pre-sable "After" EchoMRI
               FE <- as.Date("2025-2-8") # Date of final pre-sable "After" EchoMRI
               
               ifelse (d >= WS & d < SE, "Before",
                       ifelse (d >= SS & d <= FE, "After"))
             }
             #Change to standard date format
             SampleData$date <- as.Date(SampleData$SightDate, format = '%m/%d/%Y')
             #Make date of the same year i.e 2016
             SampleData$date <- as.Date(format(SampleData$date, "2016-%m-%d"))
             #Get season for each date. 
             SampleData$SightSeason <- getSeason(SampleData$date)