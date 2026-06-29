#Using just vehicle treated mice

#### Build linear mixed model for TEE (lean) #### 
sable_TEE_adj_veh <- sable_TEE_adj %>%
group_by(SABLE, ID, GROUP)%>%
filter(DRUG=="Vehicle")%>%
ungroup()

model_TEE_lean_veh <- lmer(tee ~ SABLE * GROUP + Lean + (1 | ID), data = sable_TEE_adj_veh)
summary(model_TEE_lean_veh)

    #Confirm the number of mice in the data frame is correct
    n_distinct(sable_TEE_adj_veh$ID) #n=8 (RTIOXA47 mice removed)
    
    # Calculate estimated marginal means (emmeans) #
    emm_TEE_lean_veh <- emmeans(model_TEE_lean_veh, ~ SABLE * GROUP, cov.reduce = mean)
    emm_TEE_lean_veh_df <- as.data.frame(emm_TEE_lean_veh)
    
    # Pairwise contrasts within each GROUP (ad lib or restricted)
    contrasts_by_group_TEE_lean_veh <- contrast(emm_TEE_lean_veh, method = "pairwise", by = "GROUP")
    
    # Convert to a data frame
    contrasts_TEE_lean_veh_df <- as.data.frame(contrasts_by_group_TEE_lean_veh)
    
    # Filter for restricted group and baseline vs other time points
    restricted_contrast_TEE_lean_veh <- contrasts_TEE_lean_veh_df %>%
      filter(GROUP == "Weight cycled") %>%
      filter(contrast %in% c("Baseline - BW loss", "Baseline - BW maintenance", "Baseline - BW regain"))
    restricted_contrast_TEE_lean_veh
    
    #Filter for significant contrasts
    Sig_contrast_TEE_lean_veh <- contrasts_TEE_lean_veh_df %>%
      filter(p.value <=0.05)
    Sig_contrast_TEE_lean_veh
    
    # Pairwise contrasts within each SABLE (time point)
    contrasts_SABLE_by_group_TEE_lean_veh <- contrast(emm_TEE_lean_veh, method = "pairwise", by = "SABLE")
    
    # Convert to a data frame
    contrasts_SABLE_TEE_lean_veh_df <- as.data.frame(contrasts_SABLE_by_group_TEE_lean_veh)
    contrasts_SABLE_TEE_lean_veh_df
####Conclusion: the only sable time point with a significant difference between ad lib and restircted
      #is BW regain
    #Plot
    ####Scatter plot - Veh only-Graph predicted TEE adjusted for lean mass (NZO) ####
    #Commented out the measured values for TEE since y-axis is adj TEE
    ggplot() +
      #geom_jitter(data = sable_TEE_adj, 
      #aes(x = SABLE, y = tee, color = GROUP),
      #width = 0.2, alpha = 0.4, size = 2) +
      geom_point(data = emm_TEE_lean_veh_df,
                 aes(x = SABLE, y = emmean, color = GROUP),
                 position = position_dodge(0.2), size = 4) +
      geom_line(data = emm_TEE_lean_veh_df,
                aes(x = SABLE, y = emmean, color = GROUP, group = GROUP),
                position = position_dodge(0.2), linewidth = 1.5) + 
      scale_color_manual(values=c('#FAAC41','#5392DB')) +
      geom_errorbar(data = emm_TEE_lean_veh_df,
                    aes(x = SABLE, ymin = emmean - SE, ymax = emmean + SE, color = GROUP),
                    width = 0.15, position = position_dodge(0.2)) +
      theme_minimal(base_size = 14) +
      labs(y = "Adjusted TEE (kcal/day)", x = "Time point",
           color = "Restriction group",
           title = "Vehicle - Total energy expenditure (TEE) adj. for lean mass") +
      format.plot +
      theme(legend.position = "top", 
            plot.title = element_text(hjust=0.5), 
            axis.text = element_text( color="black", size=12),
            axis.text.x = element_text(angle = 45, hjust = 1))

    