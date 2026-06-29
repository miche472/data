#Looking at EE vs BW using ANCOVA (comparing diet groups at each time point)
  #Described in "A consensus guide to preclinical indirect calorimetry experiments"

#In the paper they do "statistical analysis by ANCOVA for energy expenditure and 
#energy intake, with body mass as a covariate." 
  #Concluded that "In both cases, the groups are statistically different. A mass × group
  #interaction effect was not significant. **P< 0.01; ***P < 0.001."

#We want to compare GROUP (restricted and ad libitum) at a given time point
  #At SABLE=baseline: Does tee differ significantly between GROUPs 
  #when including Weight as a covariate?

#ANCOVA info guide: https://www.datanovia.com/en/lessons/ancova-in-r/

#Libraries
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

####Baseline: TEE vs BW: Baseline ####

#create df that just has data only from SABLE=baseline
baseline_TEEvsBW <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  select(ID,SABLE,GROUP,DRUG,tee,RMR_kcal_day,Lean,Weight)%>%
  filter(SABLE == "baseline")%>%
  mutate_at(c("GROUP"), as.factor)

####Check ANCOVA Assumptions####
    #1. Linearity --> graph of tee vs Weight, R^2, equation
    ggscatter(
      baseline_TEEvsBW, x = "Weight", y = "tee",
      color = "GROUP", add = "reg.line")+
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      labs(
      title = "Baseline: TEE vs. Body weight",
      x = "Body Weight (g)",
      y = "TEE (kcal/day)",
      color = "Restriction group") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = GROUP))
    #Ad lib: y=6.5+0.095x , R^2=0.22
    #Restricted: y=9.1+0.05x, R^2=0.055
    #Based on the flat slope and low R^2, the linearity assumption isn't really met

    #2. Homogeneity of regression slopes
    baseline_TEEvsBW %>% anova_test(tee ~ GROUP*Weight)
    #interaction term wasn't statistically significant, F(1,12)=0.16, p=0.694
    #no significant interaction between the covariate (Weight) & grouping variable (GROUP)

    #3.Normality of residuals
    #Fit the model, the covariate goes first
    model <- lm(tee ~ Weight + GROUP, data = baseline_TEEvsBW)
    # Inspect the model diagnostic metrics
    model.metrics <- augment(model) %>%
    select(-.hat, -.sigma, -.fitted) # Remove details
    head(model.metrics, 3)
    #Shapiro wilk test ()
    shapiro_test(model.metrics$.resid)
    #p=0.976 --> Shapiro test wasn't significant --> assume normality of residuals

    #4. Homogeneity of variance
    model.metrics %>% levene_test(.resid ~ GROUP)
    #p=0.911 --> Levene test wasn't significant --> assume homogeneity of the 
    #residual variances for all groups.

    #5. Check for outliers
    model.metrics %>% 
    filter(abs(.std.resid) > 3) %>%
    as.data.frame()
    #No outliers (no cases with standardized residuals greater than 3 in absolute value)

####ANCOVA Computation####
    
    #Format --> anova_test(y-axis variable ~ covariate + grouping variable)
    res.aov <- baseline_TEEvsBW %>% anova_test(tee ~ Weight + GROUP)
    get_anova_table(res.aov)
    #After controlling for BW, there is a significant difference in TEE between
    #restricted and ad libitum at baseline
    
    #Conclusion: All assumption are met except the linearity assumption. ANCOVA results say that
    #after controlling for BW, there is a sig dif in TEE between restricted and ad lib at baseline
    
    #Does the fact that the linearity assumption isn't met mean that we can't draw conclusions from this?
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

####Peak obesity: TEE vs BW ####
#create df that just has data only from SABLE=peak obesity
peak_TEEvsBW <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  select(ID,SABLE,GROUP,DRUG,tee,RMR_kcal_day,Lean,Weight)%>%
  filter(SABLE == "peak obesity")%>%
  mutate_at(c("GROUP"), as.factor)

####Check ANCOVA Assumptions####
    #1. Linearity --> graph of tee vs Weight, R^2, equation
    ggscatter(
      peak_TEEvsBW, x = "Weight", y = "tee",
      color = "GROUP", add = "reg.line")+
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      labs(
        title = "Peak obesity: TEE vs. Body weight",
        x = "Body Weight (g)",
        y = "TEE (kcal/day)",
        color = "Restriction group") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = GROUP))
    #Ad lib: y=17-0.13x , R^2=0.11
    #Restricted: y=5+0.12x, R^2=0.62
    #Based on the flat slope and low R^2, the linearity assumption isn't really met

    #2. Homogeneity of regression slopes
    peak_TEEvsBW %>% anova_test(tee ~ GROUP*Weight)
    #interaction term wasn't statistically significant, F(1,12)=3.58, p=0.083
    #no significant interaction between the covariate (Weight) & grouping variable (GROUP)

    #3. Normality of residuals
    # Fit the model, the covariate goes first
    model <- lm(tee ~ Weight + GROUP, data = peak_TEEvsBW)
    # Inspect the model diagnostic metrics
    model.metrics <- augment(model) %>%
    select(-.hat, -.sigma, -.fitted) # Remove details
    head(model.metrics, 3)
    #shapiro wilk test
    shapiro_test(model.metrics$.resid)
    #p=0.583 --> Shapiro test wasn't significant --> assume normality of residuals

    #4.Homogeneity of variance
    model.metrics %>% levene_test(.resid ~ GROUP)
    #p=0.0716 --> Levene test wasn't significant --> assume homogeneity of the residual variances for all groups.

    #5.Check for outliers
    model.metrics %>% 
    filter(abs(.std.resid) > 3) %>%
    as.data.frame()
    #No outliers (no cases with standardized residuals greater than 3 in absolute value)

####ANCOVA Computation####
    #Format --> anova_test(y-axis variable ~ covariate + grouping variable)
    res.aov <- peak_TEEvsBW %>% anova_test(tee ~ Weight + GROUP)
    get_anova_table(res.aov)
    #After controlling for BW, there is NOT a significant difference in TEE between
    #restricted and ad libitum at peak obesity F(1,13)=0.937, p=0.351.
    
    #Conclusion: All assumption are met except the linearity assumption. 
    #ANCOVA results say that fter controlling for BW, there is NOT a sig diff 
    #in TEE between restricted and ad libitum NZO mice at peak obesity
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#### BW loss: TEE vs BW ####
    
#create df that just has data only from SABLE=BW loss
    BWloss_TEEvsBW <- sable_TEE_adj_RMR %>%
      ungroup() %>%
      select(ID,SABLE,GROUP,DRUG,tee,RMR_kcal_day,Lean,Weight)%>%
      filter(SABLE == "BW loss")%>%
      mutate_at(c("GROUP"), as.factor)
    
####Check ANCOVA Assumptions####
    #1. Linearity --> graph of tee vs Weight, R^2, equation
    ggscatter(
      BWloss_TEEvsBW, x = "Weight", y = "tee",
      color = "GROUP", add = "reg.line")+
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      labs(
        title = "BW loss: TEE vs. Body weight",
        x = "Body Weight (g)",
        y = "TEE (kcal/day)",
        color = "Restriction group") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      stat_regline_equation(
        aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = GROUP))
    #Ad lib: y=7+0.039x , R^2=0.0067
    #Restricted: y= -0.57+0.19x, R^2=0.33
    #Based on the flat slope and low R^2, the linearity assumption isn't really met
    
    #2. Homogeneity of regression slopes
    BWloss_TEEvsBW %>% anova_test(tee ~ GROUP*Weight)
    #interaction term wasn't statistically significant, F(1,12)=0.468, p=0.507
    #no significant interaction between the covariate (Weight) & grouping variable (GROUP)
    
    #3.Normality of residuals
    #Fit the model, the covariate goes first
    model <- lm(tee ~ Weight + GROUP, data = BWloss_TEEvsBW)
    # Inspect the model diagnostic metrics
    model.metrics <- augment(model) %>%
      select(-.hat, -.sigma, -.fitted) # Remove details
    head(model.metrics, 3)
    #Shapiro wilk test ()
    shapiro_test(model.metrics$.resid)
    #p=0.915 --> Shapiro test wasn't significant --> assume normality of residuals
    
    #4. Homogeneity of variance
    model.metrics %>% levene_test(.resid ~ GROUP)
    #p=0.118 --> Levene test wasn't significant --> assume homogeneity of the 
    #residual variances for all groups.
    
    #5. Check for outliers
    model.metrics %>% 
      filter(abs(.std.resid) > 3) %>%
      as.data.frame()
    #No outliers (no cases with standardized residuals greater than 3 in absolute value)
    
####ANCOVA Computation####
    
    #Format --> anova_test(y-axis variable ~ covariate + grouping variable)
    res.aov <- BWloss_TEEvsBW %>% anova_test(tee ~ Weight + GROUP)
    get_anova_table(res.aov)
    #After controlling for BW, there is a significant difference in TEE between
    #restricted and ad libitum at baseline
    
    #After controlling for BW, there is NOT a significant difference in TEE between
    #restricted and ad libitum after acute BW loss F(1,13)=0.481, p=0.500.
    
    #Does the fact that the linearity assumption mean that we can't draw conclusions from this?
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#### BW maintenance: TEE vs BW ####
    
#create df that just has data only from SABLE=BW loss
    BWmain_TEEvsBW <- sable_TEE_adj_RMR %>%
      ungroup() %>%
      select(ID,SABLE,GROUP,DRUG,tee,RMR_kcal_day,Lean,Weight)%>%
      filter(SABLE == "BW maintenance")%>%
      mutate_at(c("GROUP"), as.factor)
    
####Check ANCOVA Assumptions####
    #1. Linearity --> graph of tee vs Weight, R^2, equation
    ggscatter(
    BWmain_TEEvsBW, x = "Weight", y = "tee",
    color = "GROUP", add = "reg.line")+
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      labs(
        title = "BW maintenance: TEE vs. Body weight",
        x = "Body Weight (g)",
        y = "TEE (kcal/day)",
        color = "Restriction group") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
    stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = GROUP))
    #Ad lib: y=2+0.2x , R^2=0.18
    #Restricted: y=3.3+0.13x, R^2=0.14
    #Based on the flat slope and low R^2, the linearity assumption isn't really met
    
    #2. Homogeneity of regression slopes
    BWmain_TEEvsBW %>% anova_test(tee ~ GROUP*Weight)
    #interaction term wasn't statistically significant, F(1,12)=0.105, p=0.752
    #no significant interaction between the covariate (Weight) & grouping variable (GROUP)
    
    #3.Normality of residuals
    #Fit the model, the covariate goes first
    model <- lm(tee ~ Weight + GROUP, data = BWmain_TEEvsBW)
    # Inspect the model diagnostic metrics
    model.metrics <- augment(model) %>%
      select(-.hat, -.sigma, -.fitted) # Remove details
    head(model.metrics, 3)
    #Shapiro wilk test ()
    shapiro_test(model.metrics$.resid)
    #p=0.945 --> Shapiro test wasn't significant --> assume normality of residuals
    
    #4. Homogeneity of variance
    model.metrics %>% levene_test(.resid ~ GROUP)
    #p=0.388 --> Levene test wasn't significant --> assume homogeneity of the 
    #residual variances for all groups.
    
    #5. Check for outliers
    model.metrics %>% 
      filter(abs(.std.resid) > 3) %>%
      as.data.frame()
    #No outliers (no cases with standardized residuals greater than 3 in absolute value)
    
####ANCOVA Computation####
    
    #Format --> anova_test(y-axis variable ~ covariate + grouping variable)
    res.aov <- BWmain_TEEvsBW %>% anova_test(tee ~ Weight + GROUP)
    get_anova_table(res.aov)
    #After controlling for BW, there is not a significant difference in TEE between
    #restricted and ad libitum at BW maintenance F(1,13)=0.949, p=0.348
    
    #Does the fact that the linearity assumption mean that we can't draw conclusions from this?
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#### BW regain: TEE vs BW ####
    
#Create df that just has data only from SABLE=BW regain
    BWregain_TEEvsBW <- sable_TEE_adj_RMR %>%
    ungroup() %>%
    select(ID,SABLE,GROUP,DRUG,tee,RMR_kcal_day,Lean,Weight)%>%
    filter(SABLE == "BW regain")%>%
    mutate_at(c("GROUP"), as.factor)
    
####Check ANCOVA Assumptions####
    #1. Linearity --> graph of tee vs Weight, R^2, equation
    ggscatter(
      BWregain_TEEvsBW, x = "Weight", y = "tee",
      color = "GROUP", add = "reg.line")+
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      labs(
        title = "BW regain: TEE vs. Body weight",
        x = "Body Weight (g)",
        y = "TEE (kcal/day)",
        color = "Restriction group") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      stat_regline_equation(
        aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = GROUP))
    #Ad lib: y=3.1+0.2x , R^2=0.073
    #Restricted: y=1.5+0.21x, R^2=0.18
    #Based on the flat slope and low R^2, the linearity assumption isn't really met
    
    #2. Homogeneity of regression slopes
    BWregain_TEEvsBW %>% anova_test(tee ~ GROUP*Weight)
    #interaction term wasn't statistically significant, F(1,12)=0.000487, p=0.983
    #no significant interaction between the covariate (Weight) & grouping variable (GROUP)
    
    #3.Normality of residuals
    #Fit the model, the covariate goes first
    model <- lm(tee ~ Weight + GROUP, data = BWregain_TEEvsBW)
    # Inspect the model diagnostic metrics
    model.metrics <- augment(model) %>%
      select(-.hat, -.sigma, -.fitted) # Remove details
    head(model.metrics, 3)
    #Shapiro wilk test ()
    shapiro_test(model.metrics$.resid)
    #p=0.988 --> Shapiro test wasn't significant --> assume normality of residuals
    
    #4. Homogeneity of variance
    model.metrics %>% levene_test(.resid ~ GROUP)
    #p=0.275 --> Levene test wasn't significant --> assume homogeneity of the 
    #residual variances for all groups.
    
    #5. Check for outliers
    model.metrics %>% 
      filter(abs(.std.resid) > 3) %>%
      as.data.frame()
    #No outliers (no cases with standardized residuals greater than 3 in absolute value)
    
####ANCOVA Computation####
    
    #Format --> anova_test(y-axis variable ~ covariate + grouping variable)
    res.aov <- BWregain_TEEvsBW %>% anova_test(tee ~ Weight + GROUP)
    get_anova_table(res.aov)
    #After controlling for BW, there is not a significant difference in TEE between
    #restricted and ad libitum at BW regain F(1,13)=0.378, p=0.549
    
    #Does the fact that the linearity assumption mean that we can't draw conclusions from this?
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    
#### TEEvsBW: All time points ####
#Run ancova for each time point --> couldn't figure out how to do this, so I did it manually
