####Objective is to graph distribution of glucose measurements in NZO mice####
####date is March 25 2025
# libs ----

pacman::p_load(
  tidyverse,
  googledrive,
  ggplot2,
  ggpubr,
  furrr,
  zoo,
  robustlmm,
  mmand,
  lme4,
  robustlmm,
  coxme,
  ggsurvfit,
  tidycmprsk
)

# change the directory to source file location
setwd(this.path::here())

#format plot
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        #    strip.text = element_blank(),
        plot.margin = unit(rep(0.2,4), "cm"), 
        #legend.position = "none",
        axis.text.x = element_text(family = "Helvetica", size = 14),
        axis.text.y = element_text(family = "Helvetica", size = 14),
        axis.title.y.left =element_text(family = "Helvetica", size = 16),
        axis.title.x.bottom =element_text(family = "Helvetica", size = 16))

# import data ----
## remember to set the path to script location
META <- read_csv("../data/META.csv")
BW_RAW <- read_csv("../data/BW.csv")
FI_RAW <- read_csv("../data/FI.csv")
ECHOMRI_RAW <- read_csv("../data/echomri.csv")
food_desc <- read_csv("../data/food_description.csv")


#### Create data frame with COHORT.csv files and meta data####

cohort_csv_files <- tibble(
  filepath = list.files("../data", full.names = TRUE)
) %>% 
  filter(
    grepl("COHORT_[0-9]+[0-9]*.csv", filepath) #now we can used cohort > 10
  )
cohort_csv_files

cohort_open_files2 <- cohort_csv_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map_dfr(
    ., function(X){
      read_csv(X$filepath) %>% 
        select(ID, FOOD_WEIGHT_START_G, FOOD_WEIGHT_END_G, DATE, DIET, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS,FASTED_GLU_mg_dL,BASAL_GLU_mg_dL) %>% 
        mutate(
          INTAKE_GR = (FOOD_WEIGHT_START_G - FOOD_WEIGHT_END_G),
          DATE = lubridate::mdy(DATE)
        ) %>% 
        select(ID, INTAKE_GR, DATE, BODY_WEIGHT_G, DIET_FORMULA,COMMENTS,FASTED_GLU_mg_dL,BASAL_GLU_mg_dL) %>% 
        rename(
          BW = BODY_WEIGHT_G
        ) %>% 
        mutate(BW=as.numeric(BW))
    }
  )

cohort_open_filesNZO <- cohort_open_files2 %>%
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>%
  left_join(., META, by="ID") %>%
  filter(COHORT >2 & COHORT <6)


####This works! Retains BASAL_GLU and FASTED_GLU columns and creates a column that combines them####
NZO_glucose_step1 <- cohort_open_filesNZO %>%
  select(ID, INTAKE_GR, DATE, BW, DIET, FASTED_GLU_mg_dL, BASAL_GLU_mg_dL, COHORT, STRAIN, SEX) %>%
  ungroup() %>%
  mutate(Basal = BASAL_GLU_mg_dL) %>%
  mutate(Fasted = FASTED_GLU_mg_dL) %>%
  unite(GLUCOSE_mg_dL, Basal, Fasted, sep = ",", na.rm=TRUE) %>%
  filter(GLUCOSE_mg_dL >0)

NZO_glucose_plot <- NZO_glucose_step1 %>%
  mutate(as.factor(ID)) %>%
  group_by(ID)


#Graph blood glucose distribution
  ggplot(NZO_glucose_plot, aes(x = DATE, y = GLUCOSE_mg_dL)) + 
    geom_point(aes(color=factor(ID))) + facet_wrap("ID")
    theme(legend.position = "right") + labs(x="Date") +labs(y= "Blood glucose (mg/dL)")
  
  
  