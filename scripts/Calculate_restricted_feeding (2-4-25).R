#Average food consumption calculations for food restriction

pacman::p_load(
  tidyverse,
  ggplot2,
  ggpubr,
  furrr,
  zoo,
  TTR
)

#Excerpt of code that I used to sort NZO cohorts 3-5 by adiposity index
#Recall that you made "Sorted_by_adiposity.csv" in python. It is stored in GitHub folder of Documents
open_files_AI <- read_csv("/Users/laurenmichels/Documents/GitHub/data/scripts/Sorted_by_adiposity.csv") %>% 
  select(-`...1`)

# base R way
sort_files_AI <- open_files_AI %>%
  mutate (pythongroup = as.character(group)) 
sort_files_AI$pythongroup[sort_files_AI$pythongroup == '-1'] <- 'Restrict'
sort_files_AI$pythongroup[sort_files_AI$pythongroup == '1'] <- 'Ad_libitum'

# tidyverse example
var <- open_files_AI %>% 
  mutate(
    restriction_group = if_else(group == 1, "Ad_Libitum", "Restrict")
  ) %>% 
  rename(echomri_date = Date)

#Selecting variables isn't necessary here because I am not working with a large data file
Group_assignment_AI <- sort_files_AI %>%
  select(ID, pythongroup, adiposity_index, Fat, Lean, Weight, Date, COHORT) %>%
  rename(Restriction_group = pythongroup) %>%
  rename(EchoMRI_Date = Date)

#Add the Group_assignment_AI data to FI.csv data so that you have FI and the diet group for each mouse (ie restricted or ad lib) 
#Calculate rolling average of daily food intake using ad libitum diet group
FI <-read.csv("/Users/laurenmichels/Documents/GitHub/data/data/FI.csv")
df <-FI %>%
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, DATE, corrected_intake_gr) %>%
  drop_na(corrected_intake_gr) %>% 
  filter(ID > 3705, ID<3730) %>%
  mutate(moving_avg = rollmean(corrected_intake_gr, k=3, fill=NA, align='right')) %>%
  filter(DATE > '2025-01-21') %>%
  left_join(., Group_assignment_AI, by = "ID") %>%
  filter(Restriction_group == "Ad_libitum")

# confirm the rollmean using spaghetti plots (verify no errors in data entry)
#Need multiple aes because you want multiple lines. 
#Black line shows raw corrected intake values. Red line shows moving average
#facet_wrap gives you multiple graphs (one for each ID rather than showing one compiled graph)
sp1 <- df %>% 
  mutate(date = lubridate::ymd(DATE)) %>% 
  ggplot(aes(
    date, corrected_intake_gr
  )) +
  geom_line(aes(group=ID)) +
  geom_line(aes(date, moving_avg, group=ID), color="red") +
  facet_wrap(~ID, scale="free_y")
sp1

#Restricted group will be fed 60% of the average daily food consumed by free fed mice in the past three measurements

# taking the last rolling mean measurement of each animal (df only includes Ad libitum mice)
#Slice_tail with n=1 selects the most recent rolling average value for each mouse. If you did n=2 it would select the two most recent rolling means
df2 <- df %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  slice_tail(., n=1) %>%
  ungroup() %>%
  summarise(
    last_measurement_mean = mean(moving_avg)
  ) %>% 
  mutate(
    restricted_daily_food_gr = last_measurement_mean*0.6,
    restricted_per_meal_gr = restricted_daily_food_gr*0.5
  )
df2
#Running this file will show the spaghetti plot and in the console it will give the mass 
#of food to feed the restricted mice for each meal
