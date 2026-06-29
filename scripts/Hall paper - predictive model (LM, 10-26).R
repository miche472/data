#New Kevin Hall attempt
#Goal is to recreate the table 1 and figure 2 from the Kevin Hall paper about 
  #metabolic adaptation in the biggest loser contestants

#use this data set that I created in Figure7b-RMR_correctedbyLean
sable_TEE_adj_RMR

#### Table 1 from The Biggest Loser paper (Johannsen et al. 2012) ####
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#### Restricted mice ####
#This code works
# --- 1. Create the baseline model using restricted, baseline mice ---
Hall_baseline <- sable_TEE_adj_RMR %>%
  filter(SABLE == "baseline", GROUP == "restricted")

Hall_baseline_lm <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_baseline)
summary(Hall_baseline_lm)
tidy(Hall_baseline_lm)

# --- 2. Apply model to BW loss, BW maintenance, & BW regain restricted mice ---
Hall_BWloss_BWmaint_BWregain_pred <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   # <-- THIS FIXES THE ERROR
  mutate(
    RMR_pred = predict(Hall_baseline_lm, newdata = .),
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#### Ad libitum mice (model should predict better i think) --> same code as above, but for ad lib mice####
#This works, but I should make sure that everything is named in a logical way
# --- 1. Create the baseline model using ad lib, baseline mice ---
Hall_baseline2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "baseline", GROUP == "ad lib")

Hall_baseline_lm2 <- lm(RMR_kcal_day ~ Lean + Fat, data = Hall_baseline2)
summary(Hall_baseline_lm2)
tidy(Hall_baseline_lm2)

# --- 2. Apply model to BW loss, BW maintenance, & BW regain ad lib mice ---
Hall_BWloss_BWmaint_BWregain_pred2 <- sable_TEE_adj_RMR %>%
  filter(GROUP == "ad lib", SABLE %in% c("BW loss", "BW maintenance", "BW regain")) %>%
  ungroup() %>%   # <-- THIS FIXES THE ERROR
  mutate(
    RMR_pred = predict(Hall_baseline_lm2, newdata = .),
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
#Now try to make graph in figure 2

#Regression model without lean mass
#### Restricted mice ####
# --- 1. Create the baseline model using restricted, baseline mice ---
Hall_baseline_reg_fig2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "baseline", GROUP == "restricted")

Hall_baseline_reg_fig2_lm <- lm(RMR_kcal_day ~ Lean, data = Hall_baseline_reg_fig2)
summary(Hall_baseline_reg_fig2_lm)
tidy(Hall_baseline_reg_fig2_lm)
#I think that this is the regression line that i need to draw on the graph
#{insert plot of the regression line described above}

#Goal: Adjust the measured RMR values by creating predictions using a linear regression with 
    #FM as the independent variable (can add age as an additional independent variable later)
#Approach: Create a regression for each SABLE time point (baseline, BW loss, BW maintenance, BW regain)
    # baseline_points_pred_lm <-lm(RMR_kcal_day ~ Fat, data=baseline_points_pred_2)
    # BWloss_points_pred_lm <-lm(RMR_kcal_day ~ Fat, data=BWloss_points_pred_2)
    # BWmain_points_pred_lm <-lm(RMR_kcal_day ~ Fat, data=BWmain_points_pred_2)
    # BWregain_points_pred_lm <-lm(RMR_kcal_day ~ Fat, data=BWregain_points_pred_2)
#Predicted values of adjusted RMR generated by the models for each mouse can then be 
#paired with the mouse's Lean measured at the time used in the model that produced the predicted value
#For the points on the graph x=Lean of the mouse at the SABLE time point & 
#y= adj. RMR from lm(RMR_kcal_day ~ Fat, data=___)
 
#Data POINTS for Figure 2 BASELINE (Get coordinates for mice at SABLE=baseline)
#Sort data to include only baseline restricted mice
baseline_points_pred_2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "baseline", GROUP == "restricted")

#Create linear model for baseline restricted mice using Fat as the independent variable to predict RMR
baseline_points_pred_2_lm <- lm(RMR_kcal_day ~ Fat, data = baseline_points_pred_2)
summary(baseline_points_pred_2_lm)
tidy(baseline_points_pred_2_lm)

#Apply model to baseline restricted mice
get_baseline_y_values <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE =="baseline") %>%
  ungroup() %>%
  mutate(
    RMR_pred = predict(baseline_points_pred_2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

#Inspect the results
get_baseline_y_values %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Data POINTS for Figure 2 BW LOSS (Get coordinates for mice at SABLE=BW loss)
#Sort data to include only BW loss restricted mice
BWloss_points_pred_2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "BW loss", GROUP == "restricted")

#Create linear model for BW loss restricted mice using Fat as the independent variable to predict RMR
BWloss_points_pred_2_lm <- lm(RMR_kcal_day ~ Fat, data = BWloss_points_pred_2)
summary(BWloss_points_pred_2_lm)
tidy(BWloss_points_pred_2_lm)

#Apply model to BW loss restricted mice
get_BWloss_y_values <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE =="BW loss") %>%
  ungroup() %>%
  mutate(
    RMR_pred = predict(BWloss_points_pred_2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

#Inspect the results
get_BWloss_y_values %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Data POINTS for Figure 2 BW MAINTENANCE (Get coordinates for mice at SABLE=BW maintenance)
#Sort data to include only BW maintenance restricted mice
BWmain_points_pred_2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "BW maintenance", GROUP == "restricted")

#Create linear model for BW maintenance restricted mice using Fat as the independent variable to predict RMR
BWmain_points_pred_2_lm <- lm(RMR_kcal_day ~ Fat, data = BWmain_points_pred_2)
summary(BWmain_points_pred_2_lm)
tidy(BWmain_points_pred_2_lm)

#Apply model to BW maintenance restricted mice
get_BWmain_y_values <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE =="BW maintenance") %>%
  ungroup() %>%
  mutate(
    RMR_pred = predict(BWmain_points_pred_2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

#Inspect the results
get_BWmain_y_values %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

#Data POINTS for Figure 2 BW REGAIN (Get coordinates for mice at SABLE=BW regain)
#Sort data to include only BW regain restricted mice
BWregain_points_pred_2 <- sable_TEE_adj_RMR %>%
  filter(SABLE == "BW regain", GROUP == "restricted")

#Create linear model for BW regain restricted mice using Fat as the independent variable to predict RMR
BWregain_points_pred_2_lm <- lm(RMR_kcal_day ~ Fat, data = BWregain_points_pred_2)
summary(BWregain_points_pred_2_lm)
tidy(BWregain_points_pred_2_lm)

#Apply model to BW regain restricted mice
get_BWregain_y_values <- sable_TEE_adj_RMR %>%
  filter(GROUP == "restricted", SABLE =="BW regain") %>%
  ungroup() %>%
  mutate(
    RMR_pred = predict(BWregain_points_pred_2_lm, newdata = .),
    residual = RMR_kcal_day - RMR_pred)

#Inspect the results
get_BWregain_y_values %>%
  select(ID, SABLE, RMR_kcal_day, RMR_pred, residual) %>%
  arrange(SABLE, ID)

####Prepare for graphing (recreateing figure 2)

#Combine the data frames from the four time points ("RMR adjusted for Fat")####
#The RMR_pred for each SABLE time point were generated using different regression
#equations, so I should probably only use this combined df for easy graphing
All_points_pred_Fig2 <- bind_rows(
  get_baseline_y_values,
  get_BWloss_y_values,
  get_BWmain_y_values,
  get_BWregain_y_values)

#Confirm that All_points_pred_Fig2 has restricted mice from the 4 SABLE times
All_points_pred_Fig2 %>%
  count(GROUP, SABLE)
#Verify structure
glimpse(All_points_pred_Fig2)


#graph y=RMR_pred, y=Lean, and use color or fill to show the different SABLE time points
#also add the regression line that I generated in the first step to the graph
#I should make 3 graphs: 
#baseline and BW loss
#baseline and BW maintenance
#baseline and BW regain

#Prepare data
Plot_All_points_pred_Fig2 <- All_points_pred_Fig2 %>%
  ungroup() %>%
  group_by(ID, SABLE, GROUP) %>%
  filter(GROUP=="restricted", SABLE %in% c("baseline", "BW loss"))

#make graph (still need to add regression line and add line connecting two time points for each ID)
plot_baseline_BWloss <- Plot_All_points_pred_Fig2 %>%
  ggplot(aes(x = Lean, y = RMR_pred, group=ID, fill = SABLE)) + 
  geom_point() +
theme_minimal() +
  labs(y = "Adjusted RMR (kcal/day)", x="Lean (g)", fill = "Time point") 

plot_baseline_BWloss
  
 


