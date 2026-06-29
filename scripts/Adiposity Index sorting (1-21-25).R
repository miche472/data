#Start with csv file generated in python called "Sorted_by_adiposity"
pacman::p_load(
  tidyverse,
  ggplot2,
  ggpubr,
  furrr
)

open_files_AI <- read.csv("/Users/laurenmichels/Documents/GitHub/data/scripts/Sorted_by_adiposity.csv")

sort_files_AI <- open_files_AI %>%
  mutate (pythongroup = as.character(group)) 
  sort_files_AI$pythongroup[sort_files_AI$pythongroup == '-1'] <- 'Restrict'
  sort_files_AI$pythongroup[sort_files_AI$pythongroup == '1'] <- 'Ad_libitum'
  
Group_assignment_AI <- sort_files_AI %>%
  select(ID, pythongroup, adiposity_index, Fat, Lean, Weight, Date, COHORT) %>%
  rename(Restriction_group = pythongroup) %>%
  rename(EchoMRI_Date = Date)

#Summary statistics for adiposity index for the two groups within pythongroup: Restrict and Ad-libitum
tapply(sort_files_AI$adiposity_index, sort_files_AI$pythongroup, summary)

#T-test to see if the mean adiposity index is significantly different for Restrict and Ad libitum groups
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  
#Get summary stats
sort_files_AI %>%
    group_by(pythongroup) %>%
    get_summary_stats(adiposity_index, type = "mean_sd")
  
##Test pre-requisites for T Test and perform T test
#Test for outliers  
sort_files_AI %>%
    group_by(pythongroup) %>%
    identify_outliers(adiposity_index)

#Test for normality
sort_files_AI %>%
  group_by(pythongroup) %>%
  shapiro_test(adiposity_index)

#qqplot by group
ggqqplot(sort_files_AI, x = "adiposity_index", facet.by = "pythongroup")

#Equal variance (p value is 0.623 suggesting that there is not a significant difference in the variance within each group)
sort_files_AI %>% levene_test(adiposity_index ~ pythongroup)

#Student t test (assumes equal variance)
stat.test2 <- sort_files_AI %>%
  t_test(adiposity_index ~ pythongroup, var.equal = TRUE) %>%
  add_significance()
stat.test2

#Graph adiposity index for Restrict and Ad libitum group
ggplot(sort_files_AI, aes(x = pythongroup, y = adiposity_index, fill = pythongroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) +
    theme(legend.position = "right") + labs(x=NULL) +labs(y= "Adiposity Index (fat mass/lean mass)")
  
