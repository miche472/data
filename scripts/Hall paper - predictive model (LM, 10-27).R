
#Metabolic adaptation: Goal is to recreate the Table 1 and Figure 2 from 
#the 2012 paper on the biggest loser contestants 
  #Johannsen et al. 2012 (also Kevin Hall)

#Libraries
library(dplyr)
library(ggplot2)
library(ggrepel)  # for non-overlapping text labels
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)

#use this data frame that I created in Figure7b-RMR_correctedbyLean (LM, 10-21).R
sable_TEE_adj_RMR

#### Table 1 from The Biggest Loser paper (Johannsen et al. 2012) ####
  #This section uses NZO "peak obesity" to make the RMR prediction model. 
  #Rationale: In the paper "baseline" refers to when the contestants were obese. 
  #The equivalent stage in the NZO mice is "peak obesity"
  #(In our experiment "baseline" refers to the state prior to when the mice developed obesity)
  #This section of code intends to mimic Table 1 from Johannsen, so 
  #I use NZO "peak obesity" for RMR prediction model since that mimics the paradigm in the paper
#Nomenclature: For dfs that have only ad lib mice I added a "2"
    #to the end of the df's name to distinguish from dfs with just restricted mice
#----------------------------------------------------------------------------------
#### Restricted mice ####
#This code works
# --- 1. Create the RMR prediction model using restricted, "peak obesity" ---
Hall_peakobesity <- sable_TEE_adj_RMR %>%
  filter(SABLE == "Peak obesity", GROUP == "Restricted")

Hall_peakobesity_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_peakobesity)
summary(Hall_peakobesity_lm)
tidy(Hall_peakobesity_lm)

# --- 2. Apply model to BW loss, BW maintenance, & BW regain restricted mice ---
Hall_BWloss_BWmaint_BWregain_pred <- sable_TEE_adj_RMR %>%
  filter(GROUP == "Restricted", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   
  mutate(
    RMR_pred = predict(Hall_peakobesity_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

# --- 3. Inspect the results ---
Hall_BWloss_BWmaint_BWregain_pred %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)



#Optional:Split results for BWlosss and BW maintenance into three data frames
Hall_BWloss_pred <- Hall_BWloss_BWmaint_BWregain_pred %>% 
  filter(SABLE == "BW loss")
Hall_BWmaint_pred <- Hall_BWloss_BWmaint_BWregain_pred %>% 
  filter(SABLE == "BW maintenance")
Hall_BWregain_pred <- Hall_BWloss_BWmaint_BWregain_pred %>% 
  filter(SABLE == "BW regain")

####New Get predictions from peak obesity as well (which is equivalent to week 0 in Johannsen et al.) ####
####left off here...I don't think this is totally correct. I should use an ANOVA that accounts for the fact
#that there were repeated measures, but 
# --- 1. Create the RMR prediction model using restricted, "peak obesity" ---
Hall_peakobesity <- sable_TEE_adj_RMR %>%
  filter(SABLE == "Peak obesity", GROUP == "Restricted")

Hall_peakobesity_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_peakobesity)
summary(Hall_peakobesity_lm)
tidy(Hall_peakobesity_lm)

# --- 2. Apply model to BW loss, BW maintenance, & BW regain restricted mice ---
Hall_peakobesity_BWloss_BWmaint_BWregain_pred <- sable_TEE_adj_RMR %>%
  filter(GROUP == "Restricted", SABLE %in% c("Peak obesity", "BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   
  mutate(
    RMR_pred = predict(Hall_peakobesity_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

# --- 3. Inspect the results ---
Hall_peakobesity_BWloss_BWmaint_BWregain_pred %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#1. Check assumptions for ANOVA?
# Histogram and QQ plot
ggplot(Hall_peakobesity_BWloss_BWmaint_BWregain_pred, aes(x = residual)) +
  geom_histogram(bins = 20, color = "black", fill = "lightblue") +
  facet_wrap(~SABLE)

# QQ plot by group
ggplot(Hall_peakobesity_BWloss_BWmaint_BWregain_pred, aes(sample = residual)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~SABLE)

#2a. ANOVA to test for differences in mean RMR at different SABLE time points (for restricted mice)
residual_aov <- aov(residual ~ SABLE, data = Hall_peakobesity_BWloss_BWmaint_BWregain_pred)
summary(residual_aov)
#p<0.0001, so do post hoc test

#2b. Post hoc analysis to identify which time points have different mean RMR
TukeyHSD(residual_aov)
#Restricted BW loss and Restricted BW regain, p=0.000127 --> sig diff
#Restricted BW maintenance and Restricted BW regain, p=0.00112 --> sig diff

#3a. ANOVA to test for differences in mean measured RMR at different SABLE time points (restricted)
RMRpred_aov <- aov(RMR_pred ~ SABLE, data = Hall_peakobesity_BWloss_BWmaint_BWregain_pred)
summary(RMRpred_aov)

#3b. Do post hoc test to see which time points there are differences in residuals
TukeyHSD(RMRpred_aov)

#possible way to do repeated measures anova
#For measured RMR
RMR_kcal_day_aov <- aov(RMR_kcal_day ~ SABLE + Error(ID/SABLE), 
                        data=Hall_peakobesity_BWloss_BWmaint_BWregain_pred)
summary(RMR_kcal_day_aov)
pwc <- RMR_kcal_day_aov%>%
pairwise.t.test(
  x = RMR_kcal_day,
  g = SABLE,
  p.adjust.method = 'bonferroni', data=Hall_peakobesity_BWloss_BWmaint_BWregain_pred)


#----------------------------------------------------------------------------------
#### Ad libitum mice ####
#Same code as above, but for ad lib mice
#Uses paradigm from paper (NZO at peak obesity are used to make RMR prediction model)
# --- 1. Create the RMR prediction model using ad lib, "peak obesity" ---
Hall_peakobesity2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "Peak obesity", GROUP == "Ad lib")

Hall_peakobesity2_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_peakobesity2)
summary(Hall_peakobesity2_lm)
tidy(Hall_peakobesity2_lm)

# --- 2. Apply model to BW loss, BW maintenance, & BW regain ad lib mice ---
Hall_BWloss_BWmaint_BWregain_pred2 <- sable_TEE_adj_RMR %>%
  filter(GROUP == "Ad lib", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   
  mutate(
    RMR_pred = predict(Hall_peakobesity2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

# --- 3. Inspect the results ---
Hall_BWloss_BWmaint_BWregain_pred2 %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Optional: Split results for BWlosss and BW maintenance into two data frames (might not be necessary, but helps to check)
Hall_BWloss_pred2 <- Hall_BWloss_BWmaint_BWregain_pred2 %>% 
  filter(SABLE == "BW loss")
Hall_BWmaint_pred2 <- Hall_BWloss_BWmaint_BWregain_pred2 %>% 
  filter(SABLE == "BW maintenance")
Hall_BWregain_pred2 <- Hall_BWloss_BWmaint_BWregain_pred2 %>% 
  filter(SABLE == "BW regain")
#----------------------------------------------------------------------------------
####Combine ad lib and restricted into one table####
  #Remember that the baseline equations used to calculate the predicted 
  #RMR values were different for the two diet restriction groups --> this table 
  #should really be used for stats...just for creating a summary table
Hall_allSABLE_predictions <- bind_rows(
  Hall_BWloss_BWmaint_BWregain_pred,
  Hall_BWloss_BWmaint_BWregain_pred2)

#Confirm that ad lib and restricted mice are both included in this combined df
Hall_allSABLE_predictions %>%
  count(GROUP, SABLE)
#Verify structure
glimpse(Hall_allSABLE_predictions)

#Clean data frame for csv export
Hall_allSABLE_predictions_for_csv <-Hall_allSABLE_predictions %>%
select(-window_start_time, -window_end_time) %>%
  rename(TEE = tee,
         RMR_measured = RMR_kcal_day)

#write_csv
write_csv(x = Hall_allSABLE_predictions_for_csv, "../data/Hall_Table1_pred.csv")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#### Modified Table 1 from The Biggest Loser paper (Johannsen et al. 2012) ####
  #modified to use the baseline RMR values (pre-obesity) to make the regression model
  #Nomenclature: To indicate that this code is modified from the chunk above (it uses
  #baseline NZO data to genereate RMR prediction model) I added a "b" to the end of 
  #the name of each df created. Also, for dfs that have only ad lib mice I added a "2"
  #to the end of the name to distinguish from dfs with just restricted mice
#----------------------------------------------------------------------------------
#### Restricted mice ####
#This code works
# --- 1b. Create the baseline model using restricted, baseline mice ---
Hall_baselineb <- sable_TEE_adj_RMR %>%
  filter(SABLE == "Baseline", GROUP == "Restricted")

Hall_baselineb_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_baselineb)
summary(Hall_baselineb_lm)
tidy(Hall_baselineb_lm)

# --- 2b. Apply model to BW loss, BW maintenance, & BW regain restricted mice ---
Hall_BWloss_BWmaint_BWregain_predb <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   
  mutate(
    RMR_pred = predict(Hall_baselineb_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

# --- 3b. Inspect the results ---
Hall_BWloss_BWmaint_BWregain_predb %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Optional:Split results for BWlosss and BW maintenance into three data frames
Hall_BWloss_predb <- Hall_BWloss_BWmaint_BWregain_predb %>% 
  filter(SABLE == "BW loss")
Hall_BWmaint_predb <- Hall_BWloss_BWmaint_BWregain_predb %>% 
  filter(SABLE == "BW maintenance")
Hall_BWregain_predb <- Hall_BWloss_BWmaint_BWregain_predb %>% 
  filter(SABLE == "BW regain")
#----------------------------------------------------------------------------------
#### Ad libitum mice ####
#Same code as above, but for ad lib mice instead of restricted mice
#I think the model for ad lib mice should predict RMR better (i.e. smaller residuals)
#because they didn't undergo BW loss and therefore shouldn't display
#metabolic adaptation 
# --- 1b. Create the baseline model using ad lib, baseline mice ---
Hall_baselineb2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "baseline", GROUP == "ad lib")

Hall_baselineb2_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_baselineb2)
summary(Hall_baselineb2_lm)
tidy(Hall_baselineb2_lm)

# --- 2b. Apply model to BW loss, BW maintenance, & BW regain ad lib mice ---
Hall_BWloss_BWmaint_BWregain_predb2 <- sable_TEE_adj_RMR %>%
  filter(GROUP == "ad lib", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   
  mutate(
    RMR_pred = predict(Hall_baselineb2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

# --- 3b. Inspect the results ---
Hall_BWloss_BWmaint_BWregain_predb2 %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Optional: Split results for BWlosss and BW maintenance into two data frames (might not be necessary, but helps to check)
Hall_BWloss_predb2 <- Hall_BWloss_BWmaint_BWregain_predb2 %>% 
  filter(SABLE == "BW loss")
Hall_BWmaint_predb2 <- Hall_BWloss_BWmaint_BWregain_predb2 %>% 
  filter(SABLE == "BW maintenance")
Hall_BWregain_predb2 <- Hall_BWloss_BWmaint_BWregain_predb2 %>% 
  filter(SABLE == "BW regain")
#----------------------------------------------------------------------------------
####Combine ad lib and restricted into one table####
#Remember that the baseline equations used to calculate the predicted 
#RMR values were different for the two diet restriction groups --> this table 
#should really be used for stats...just for creating a summary table
Hall_allSABLE_predictionsb <- bind_rows(
  Hall_BWloss_BWmaint_BWregain_predb,
  Hall_BWloss_BWmaint_BWregain_predb2)

#Confirm that ad lib and restricted mice are both included in this combined df
Hall_allSABLE_predictionsb %>%
  count(GROUP, SABLE)
#Verify structure
glimpse(Hall_allSABLE_predictionsb)

#Clean data frame for csv export
Hall_allSABLE_predictionsb_for_csv <-Hall_allSABLE_predictionsb %>%
  select(-window_start_time, -window_end_time) %>%
  rename(TEE = tee,
         RMR_measured = RMR_kcal_day)

#write_csv
write_csv(x = Hall_allSABLE_predictionsb_for_csv, "../data/Hall_Table1_predb.csv")


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#### Figure 2: Adjusted RMR vs lean for NZO mice ####

#Johannsen et al. called RMR at the start of the competition (i.e. obesity) "baseline" and 
    #RMR at the end of the competition (i.e. post weight loss) "wk 30" 
#Equivalent time points in NZO mice were called "peak obesity" and "BW loss", respectively
#Note: the lm used to calculate adj RMR values at peak obesity is not the same as  
      #the line displayed on the graph

#1. Data points,set 1: Calculate adj. RMR for peak obesity timepoint (RMR_pred_step1) 
      # lm(RMR_kcal_day ~ Fat)...only use RMR measured at peak obesity
      #these are y-values to plot for baseline adj RMR on the graph 
      #these y-values will be INPUTS for the regression line in step 2
#2. Regression line: Input RMR_pred_step1 calculated in step 1 into lm(RMR_pred_step1 ~ Lean)
      #fits regression line to display on graph
#3. Data points,set2: Calculate adj. RMR for BW loss (RMR_pred_step3)
      # lm(RMR_kcal_day ~ Fat)...only use RMR measured at BW loss
      #these are y-values to plot for BW loss adj RMR on the graph 
#4. Create graph with 2 sets of data points and 1 regression line:
  #Data points, set 1 --> black dots (y=RMR adj for fat which is predicted in step 1, x=lean mass at peak obesity)
        #y = RMR_pred_step1
        #x = Lean (at BW loss)
  #Data points, set 2 --> white dots (y=RMR adj for fat which is predicted in step 3, x=lean mass at BW loss)
        #y = RMR_pred_step3
        #x = Lean (at BW loss)
  #Regression line --> regression line that is fitted to data points, set 1 (i.e. adj. RMR at peak obesity)
        #(alternatively, use regression line fitted to MEASURED RMR at peak obesity)
  
#----Step 1:Data points, set 1 ----
        #Include only peak of obesity restricted mice
        RMR_pred_step1_R <- sable_TEE_adj_RMR %>%
          filter(SABLE == "Peak obesity", GROUP == "Restricted")
        
        #Create linear model for peak obesity restricted mice using Fat as the independent variable to predict RMR
        RMR_pred_step1_R_lm <- lm(RMR_kcal_day ~ Fat, data = RMR_pred_step1_R)
        summary(RMR_pred_step1_R_lm)
        tidy(RMR_pred_step1_R_lm)
        
        #Apply model to peak obesity restricted mice
        peakobesity_pred_R <- sable_TEE_adj_RMR %>%
          filter(GROUP == "Restricted", SABLE =="Peak obesity") %>%
          ungroup() %>%
          mutate(
            RMR_pred = predict(RMR_pred_step1_R_lm, newdata = .),
            residual = RMR_kcal_day - RMR_pred)
        
        #Inspect the results
        peakobesity_pred_R %>%
          select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
          arrange(SABLE, ID)

#----Step 2:Regression line (Two options)----
        #Option 1: Use Adj RMR --> use peak obesity RMR measured values ADJUSTED for FM
        RMR_step2_R_lm <- lm(RMR_pred ~ Lean, data = peakobesity_pred_R)
        summary(RMR_step2_R_lm)
        tidy(RMR_step2_R_lm)
        
        #Option 2: using unadjusted (i.e. measured values) for RMR at peak obesity
        For_unadjusted <- sable_TEE_adj_RMR %>%
        filter(SABLE == "Peak obesity", GROUP == "Restricted") 
        
        unadjusted_lm <- lm(RMR_kcal_day ~ Lean, data = For_unadjusted)
        summary(unadjusted_lm)
        tidy(unadjusted_lm)

#----Step 3:Data points, set 2 ----
        #Include only BW loss restricted mice
        RMR_pred_step3_R <- sable_TEE_adj_RMR %>%
          filter(SABLE == "BW loss", GROUP == "Restricted")
        
        #Create linear model for BW loss restricted mice using Fat as the independent variable to predict RMR
        RMR_pred_step3_R_lm <- lm(RMR_kcal_day ~ Fat, data = RMR_pred_step3_R)
        summary(RMR_pred_step3_R_lm)
        tidy(RMR_pred_step3_R_lm)
        
        #Apply model to BW loss restricted mice
        BWloss_pred_R <- sable_TEE_adj_RMR %>%
          filter(GROUP == "Restricted", SABLE =="BW loss") %>%
          ungroup() %>%
          mutate(
            RMR_pred = predict(RMR_pred_step3_R_lm, newdata = .),
            residual = RMR_kcal_day - RMR_pred)
        
        #Inspect the results
        BWloss_pred_R %>%
          select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
          arrange(SABLE, ID)

#----Step 4: Create graph ----

#Combine df containing the 2 sets of data points (peak obesity and BW loss)
        plot_data <- bind_rows(
        peakobesity_pred_R,
        BWloss_pred_R)
        
#--- Create the plot --> Choose one regression line only 
        #Option 1: adjusted RMR at peak obesity, from step 1 (RMR_step2_R_lm)
        #Option 2: unadjusted RMR at peak obesity (unadjusted_lm)
        
  #Combine dfs containing the 2 sets of data points (peak obesity and BW loss)
  #Create new variable called Timepoint to allow for easy graphing based on time point
  plot_data <- bind_rows(
  peakobesity_pred_R %>% mutate(Timepoint = "Peak obesity"),
  BWloss_pred_R %>% mutate(Timepoint = "BW loss"))
        
  # Keep only IDs present in both datasets
  ids_in_both <- intersect(peakobesity_pred_R$ID, BWloss_pred_R$ID)
  plot_data_filtered <- plot_data %>% filter(ID %in% ids_in_both)
        
  # Compute axis limits (otherwise R zooms in too much)
  xlims <- range(plot_data_filtered$Lean, na.rm = TRUE)
  ylims <- range(plot_data_filtered$RMR_pred, na.rm = TRUE)
        
  ggplot(plot_data_filtered, aes(x = Lean, y = RMR_pred)) +
    
  #Connecting lines for same mouse IDs at two time points
  geom_line(aes(group = ID), color = "gray60", linewidth = 0.7, alpha = 0.8) +
          
  #Points (fill mapped to Timepoint)
  geom_point(aes(fill = Timepoint), shape = 21, size = 3, color = "black") +
          
  #Mouse ID labels
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = Inf) +
          
  #CHOOSE ONLY ONE: Regression line from RMR_step2_R_lm (peak obesity model)
  stat_function(fun = function(x) coef(RMR_step2_R_lm)[1] + coef(RMR_step2_R_lm)[2] * x,color = "black", linewidth = 1) +
  #stat_function(fun = function(x) coef(unadjusted_lm)[1] + coef(unadjusted_lm)[2] * x,color = "black", linewidth = 1) +      
  
  #Legend: colors + order
  scale_fill_manual(name = "Time point", values = c("Peak obesity" = "black", "BW loss" = "white"),
  breaks = c("Peak obesity", "BW loss")) +
  
  #Axis limits
  coord_cartesian(xlim = xlims, ylim = ylims) +
    
  #Labels and theme
  labs(x = "Lean mass (g)", y = "Adjusted RMR (kcal/day)",
  title = "Peak obesity & BW loss (regression is adjusted RMR at peak obesity)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",legend.title = element_text(size = 12),legend.text = element_text(size = 12))
  
#### Modified Figure 2 to address the EE gap more directly ####
  #Question: Can an EE gap be observed in the NZO mice between baseline and BW regain?
  #Variation of Johannsen 2012 Fig 2 -> addresses baseline vs. BW regain
  #Used the same approach as described in the chunk directly above, but replaced
  #the time points: peak obesity-->baseline and BW loss-->BW regain
  #Nomenclature: added a "b" to each step and df relative to the code chunk above
  
  #----Step 1b:Data points, set 1 (baseline) ----
  #Include only baseline restricted mice
  RMR_pred_step1b_R <- sable_TEE_adj_RMR %>%
    filter(SABLE == "Baseline", GROUP == "Restricted")
  
  #Create linear model for baseline restricted mice using Fat as the independent variable to predict RMR
  RMR_pred_step1b_R_lm <- lm(RMR_kcal_day ~ Fat, data = RMR_pred_step1b_R)
  summary(RMR_pred_step1b_R_lm)
  tidy(RMR_pred_step1b_R_lm)
  
  #Apply model to baseline restricted mice
  baseline_predb_R <- sable_TEE_adj_RMR %>%
    filter(GROUP == "Restricted", SABLE =="Baseline") %>%
    ungroup() %>%
    mutate(
      RMR_pred = predict(RMR_pred_step1b_R_lm, newdata = .),
      residual = RMR_kcal_day - RMR_pred)
  
  #Inspect the results
  baseline_predb_R %>%
    select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
    arrange(SABLE, ID)
  
  #----Step 2b:Regression line (Two options)----
  #Option 1: Use Adj RMR --> use baseline RMR measured values ADJUSTED for Fat
  RMR_step2b_R_lm <- lm(RMR_pred ~ Lean, data = baseline_predb_R)
  summary(RMR_step2b_R_lm)
  tidy(RMR_step2b_R_lm)
  
  #Option 2: using unadjusted (i.e. measured values) for RMR at baseline
  For_unadjustedb <- sable_TEE_adj_RMR %>%
    filter(SABLE == "Baseline", GROUP == "Restricted") 
  
  unadjustedb_lm <- lm(RMR_kcal_day ~ Lean, data = For_unadjustedb)
  summary(unadjustedb_lm)
  tidy(unadjustedb_lm)
  
  #----Step 3b:Data points, set 2 (BW regain) ----
  #Include only BW regain restricted mice
  RMR_pred_step3b_R <- sable_TEE_adj_RMR %>%
    filter(SABLE == "BW regain", GROUP == "Restricted")
  
  #Create linear model for BW regain restricted mice using Fat as the independent variable to predict RMR
  RMR_pred_step3b_R_lm <- lm(RMR_kcal_day ~ Fat, data = RMR_pred_step3b_R)
  summary(RMR_pred_step3b_R_lm)
  tidy(RMR_pred_step3b_R_lm)
  
  #Apply model to BW regain restricted mice
  BWregain_predb_R <- sable_TEE_adj_RMR %>%
    filter(GROUP == "Restricted", SABLE =="BW regain") %>%
    ungroup() %>%
    mutate(
      RMR_pred = predict(RMR_pred_step3b_R_lm, newdata = .),
      residual = RMR_kcal_day - RMR_pred)
  
  #Inspect the results
  BWregain_predb_R %>%
    select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
    arrange(SABLE, ID)
  
  #----Step 4b: Create graph ----
  #Combine df containing the 2 sets of data points (baseine and BW regain)
  plot_datab <- bind_rows(
    baseline_predb_R,
    BWregain_predb_R)
  
  #--- Create the plot --> Choose one regression line only 
  #Option 1: adjusted RMR at baseline, from step 1 (RMR_step2b_R_lm)
  #Option 2: unadjusted RMR at BW regain (unadjustedb_lm)
  
  #Combine dfs containing the 2 sets of data points (baseline and BW regain)
  #Create new variable called Timepoint to allow for easy graphing based on time point
  plot_datab <- bind_rows(
    baseline_predb_R %>% mutate(Timepoint = "Baseline"),
    BWregain_predb_R %>% mutate(Timepoint = "BW regain"))
  
  # Keep only IDs present in both datasets
  ids_in_bothb <- intersect(baseline_predb_R$ID, BWregain_predb_R$ID)
  plot_data_filteredb <- plot_datab %>% filter(ID %in% ids_in_bothb)
  
  # Compute axis limits (otherwise R zooms in too much)
  xlims <- range(plot_data_filteredb$Lean, na.rm = TRUE)
  ylims <- range(plot_data_filteredb$RMR_pred, na.rm = TRUE)
  
  ggplot(plot_data_filteredb, aes(x = Lean, y = RMR_pred)) +
    
    #Connecting lines for same mouse IDs at two time points
    geom_line(aes(group = ID), color = "gray60", linewidth = 0.7, alpha = 0.8) +
    
    #Points (fill mapped to Timepoint)
    geom_point(aes(fill = Timepoint), shape = 21, size = 3, color = "black") +
    
    #Mouse ID labels
    geom_text_repel(aes(label = ID), size = 3, max.overlaps = Inf) +
    
    #CHOOSE ONLY ONE: Regression line from RMR_step2b_R_lm (baseline model)
    stat_function(fun = function(x) coef(RMR_step2b_R_lm)[1] + coef(RMR_step2b_R_lm)[2] * x,color = "black", linewidth = 1) +
    #stat_function(fun = function(x) coef(unadjustedb_lm)[1] + coef(unadjustedb_lm)[2] * x,color = "black", linewidth = 1) +      
    
    #Legend: colors + order
    scale_fill_manual(name = "Time point", values = c("Baseline" = "black", "BW regain" = "white"),
                      breaks = c("Baseline", "BW regain")) +
    
    #Axis limits
    coord_cartesian(xlim = xlims, ylim = ylims) +
    
    #Labels and theme
    labs(x = "Lean mass (g)", y = "Adjusted RMR (kcal/day)",
         title = "Baseline & BW regain (regression is adjusted RMR at baseline)") +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",legend.title = element_text(size = 12),legend.text = element_text(size = 12))
  
  
  
  
  