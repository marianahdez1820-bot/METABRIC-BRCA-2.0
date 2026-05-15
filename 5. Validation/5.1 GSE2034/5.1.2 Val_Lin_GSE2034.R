
# 6.- Validation ----------------------------------------------------------

# This is the predict phase

gse2034_results <- predict(final_fit, new_data = proof_genes_pt_gse2034, type = "linear_pred") %>%
  bind_cols(proof_genes_pt_gse2034)


# 6.3 To get the p value

validation_test <- coxph(Surv(EVENT_MON,  EVENT_STAT) ~ .pred_linear_pred, data = gse2034_results)

summary(validation_test)

# 6.8 Calculate the actual Concordance Index

c_index_results.2034 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = gse2034_results)

library(survminer)

# Create risk groups based on the median of the predictions

gse2034_results <- gse2034_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))

# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results)

# Plot

ggsurvplot(km_fit, 
           data = gse2034_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("High risk", "Low risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 180),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3, 
           palette = c("#E41A1C", "#377EB8"),
)


gse2034_results <- gse2034_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
gse2034_results$risk_group <- relevel(gse2034_results$risk_group, ref = "Low Risk")


summary_gse2034 <- summary(coxph(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results))


library(timeROC)

# Area under the curve per time

res_auc <- timeROC(T = gse2034_results$EVENT_MON,
                   delta = gse2034_results$ EVENT_STAT,
                   marker = -gse2034_results$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 120), # 3, 5, and 10 years
                   iid = TRUE)

# 6.9.2 View the AUC values

res_auc_gse2034 <- res_auc$AUC %>% 
  as.data.frame()

# View the AUC values

print(res_auc$AUC)


# Multivariate regression cox with clinical data

proof_genes_pt.gse2034.cox <- 
  proof_genes_pt_gse2034 %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("file_name") %>% 
  mutate(SCORE = gse2034_results$.pred_linear_pred
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                SCORE
  ) %>% 
  na.omit()

independent_prog.gse2034 <- coxph(surv_obj ~ SCORE, 
                                  data = proof_genes_pt.gse2034.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)



num_param_compare <- c(9:21)

cat(paste0("The signature got a C-score of ", round(c_index_results.2034$concordance, 2)),
    paste0("an HR of ", round(summary_gse2034$coefficients[2], 2), " (CI 95% of ", round(summary_gse2034$conf.int[3], 2), " - ", round(summary_gse2034$conf.int[4], 2), " pval ", summary_gse2034$coefficients[5], ")"),
    paste0("AUC at 3 years of ", round(res_auc_gse2034[1,], 2), " at 5 years of ", round(res_auc_gse2034[2,], 2), " and at 6 years of "),
    paste0("As an independence factor it has an HR of ", round(independent_prog.gse2034$estimate[independent_prog.gse2034$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.gse2034$conf.low[independent_prog.gse2034$term == "SCORE"], 2), " - ", round(independent_prog.gse2034$conf.high[independent_prog.gse2034$term == "SCORE"], 2), " pval of ", independent_prog.gse2034$p.value[independent_prog.gse2034$term == "SCORE"], ")"),
    sep = ". "
)
