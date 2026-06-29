#"I measured total energy expenditure (tee) and lean mass (Lean) at 
#five time points (SABLE = Baseline, Peak obesity, BW loss, BW maintenance, BW regain) 
#for 16 mice. Each mouse had a different ID and was part of one diet group 
#(GROUP=Restricted, Ad lib). How can I do an ANOVA to test for a difference in 
#tee between time points for each diet group. I need to include lean mass at each 
#time point as a covariate and account for the fact that each mouse had tee and 
#Lean measured at each SABLE time point. Maybe there is a way to do a repeated 
#measures and/or paired ANOVA"

#ChatGPT suggested that I use a repeated measures ANCOVA 
#this is the code recommended

library(lme4)
library(lmerTest) # for p-values

model <- lmer(tee ~ GROUP * SABLE + Lean + (1 | ID), data = sable_TEE_adj_RMR)

#post hoc:
#To test pairwise differences between SABLE levels within each group, while controlling for Lean:
library(emmeans)

emm <- emmeans(model, ~ SABLE | GROUP, cov.reduce = mean)
pairs(emm, adjust = "tukey")
#visualize adjusted means
plot(emm, comparisons = TRUE)

#Model assumptions
plot(model)        # residual vs fitted
qqnorm(resid(model))
qqline(resid(model))

#Extract and visualize predictions from model
library(emmeans)

emm <- emmeans(model, ~ SABLE * GROUP, cov.reduce = mean)
summary(emm)

#Pairwise comparisons
#Compare adjusted tee at all five time points within each diet group 
pairs(emmeans(model, ~ SABLE | GROUP), adjust = "tukey")

#### for power analysis: this shows that there is a significant difference in adjusted tee at 
#At each time point, compare adjusted tee for the two diet groups
pairs(emmeans(model, ~ GROUP | SABLE), adjust = "tukey")
