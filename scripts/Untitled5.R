ancova_model <- aov(exam ~ technique + current_grade, data = data)



test_weight_cycled <- sable_TEE_adj_RMR %>%
  ungroup() %>%
  select(ID, SABLE, GROUP, DRUG, tee, RMR_kcal_day, Lean, Weight) %>%
  filter(GROUP == "Weight cycled") %>%
  group_by(ID) %>%
n_distinct(test_weight_cycled$ID)
  #mutate_at(c("GROUP"), as.factor) %>%

ancova_model <- aov(tee ~ SABLE + Lean, data = test_weight_cycled)
#view summary of model
Anova(ancova_model, type="III") 

library(multcomp)
install.packages("multcomp")

#define the post hoc comparisons to make
test_postHocs <- glht(ancova_model, linfct = mcp(SABLE = "Tukey"))

#view a summary of the post hoc comparisons
summary(test_postHocs)


#correct TEE for weight cycled mice for BW
ancova_model2 <- aov(tee ~ SABLE + Weight, data = test_weight_cycled)
Anova(ancova_model2, type="III") 
test_postHocs2 <- glht(ancova_model2, linfct = mcp(SABLE = "Tukey"))
summary(test_postHocs2)