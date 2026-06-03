library(survminer)
library(timeROC)


# 1.- Validation ----------------------------------------------------------

# 1.1 Predict on GSE2034 set

gse2034_results <- predict(final_fit, new_data = proof_genes_pt.gse2034, type = "linear_pred") %>%
  bind_cols(proof_genes_pt.gse2034 %>% 
              rownames_to_column("file_name"))


# 1.2 Calculate the Concordance Index

c_index_results.2034 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = gse2034_results)

# 1.3 Table with confidence interval and z stat and estimated p val

c_index_summary.gse2034 <- data.frame(
  C_Index = c_index_results.2034$concordance,
  SE = sqrt(c_index_results.2034$var)
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

print(c_index_summary.gse2034)



# 1.4 Create risk groups based on the median of the predictions

gse2034_results <- gse2034_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < true_cut$cutpoint$cutpoint[1], "High Risk", "Low Risk"))) # median(.pred_linear_pred) # true_cut$cutpoint$cutpoint[1]

# 1.5 Fit the KM curve

km_fit.gse2034 <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results)

# 1.6 Plot

ggsurvplot(km_fit.gse2034, 
           data = gse2034_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Recurrence validation GSE2034 (39 genes / untreated)",
           font.title = 20,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("High risk", "Low risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 180),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_linedraw(), # ggplot2 theme
           
           linewidth = 3, 
           palette = c("#E41A1C", "#377EB8"),
)


gse2034_results <- gse2034_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# 1.7 Divides time group into less and more than 70 months since if thsi division is not made the proportional hazards assumption is not met

gse2034_split <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ ., 
  data = gse2034_results,
  cut = 70, 
  episode = "time_group",
  id = "patient_id"
)

# 1.8 Run Cox again

gse2034_split$risk_group <- relevel(gse2034_split$risk_group, ref = "Low Risk")


summary_gse2034 <- summary(coxph(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group * strata(time_group), data = gse2034_split))

# 2.- Metric results ------------------------------------------------------

# 2. 1 Area under the curve per time

res_auc.gse2034 <- timeROC(T = gse2034_results$EVENT_MON,
                           delta = gse2034_results$ EVENT_STAT,
                           marker = -gse2034_results$pred_z,
                           cause = 1, # The event code
                           times = c(12, 36, 60, 120), # 3, 5, and 10 years
                           iid = TRUE)

# 2.1.2 Table with confidence interval and z stat and estimated p val

auc_ci.gse2034 <- data.frame(
  AUC  = res_auc.gse2034$AUC,
  SE   = res_auc.gse2034$inference$vect_sd_1,
  time = res_auc.gse2034$times
) %>% 
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

# 2.2 Text

for (i in 1:5) {
  
  cat(paste0(round(auc_ci.gse2034$AUC[i], 2), " (", round(auc_ci.gse2034$conf_int_low95[i], 2), "-", round(auc_ci.gse2034$conf_int_high95[i], 2), "z score ", round(auc_ci.gse2034$z_stat[i], 2),  " pval ", round(auc_ci.gse2034$p_value[i], 4), ")", " at ", auc_ci.gse2034$time[i], " months (",auc_ci.gse2034$time[i] / 12, " years), "))
  
}

# 2.3 View the AUC values
auc_gse2034 <- res_auc.gse2034$AUC 
print(auc_gse2034)


gse2034_results <- 
  gse2034_results %>% 
  mutate(EVENT_STAT = factor(EVENT_STAT))

# 2.4 Object to later plot ROC curves

global_roc.gse2034 <- roc_curve(gse2034_results,
                                EVENT_STAT,
                                .pred_linear_pred
) %>%
  mutate(label = "GSE2034")


# 2.4.2 Combine all time points into a long dataframe

plot_roc.gse2034 <- map_df(c(12, 36, 60, 120), function(i) { # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc.gse2034$FP[, time_label],
    TP = res_auc.gse2034$TP[, time_label],
    Time = factor(i),
    data_set = "GSE2034"
  )
})


# 2.4.3 Labels for faceted plot

facet_labels.gse2034 <- data.frame(
  Time = factor(c(12, 36, 60, 120)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(res_auc.gse2034$AUC[1:4]), 3), "%"),
  data_set = "GSE2034"
)

# 3.4 Plot with facet

ggplot(plot_roc.gse2034, aes(x = FP, y = TP)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Time, ncol = 2) +
  theme_bw() +
  labs(
    title = "ROC Curves by Time Point",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  geom_text(data = facet_labels.gse2034, 
            aes(x = 0.75, y = 0.1, label = AUC_Text), 
            size = 4, fontface = "bold")


# 3.- Cox models  ----------------------------------------------------------

# 3.1 Asigning the correspondant metadata to the tested patients

proof_genes_pt.gse2034.cox <- 
  proof_genes_pt.gse2034 %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("file_name") %>% 
  mutate(SCORE = gse2034_results$.pred_linear_pred
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                SCORE,
                EVENT_MON,
                EVENT_STAT
  ) %>% 
  na.omit()

# 3.1.2 Once again split to test the cox with score as a continuous variable and not as a divided low and high risk

gse2034_split.cox <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ ., 
  data = proof_genes_pt.gse2034.cox,
  cut = 70, 
  episode = "time_group",
  id = "patient_id"
)

# 3.2 Fit the model with the interaction

cox_model.gse2034 <- coxph(
  Surv(tstart, EVENT_MON, EVENT_STAT) ~ SCORE:strata(time_group), 
  data = gse2034_split.cox
)

# 3.3 Tidy

independent_prog.gse2034 <- 
  cox_model.gse2034 %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

# 3.3.2 Clean names for table

supplementary_table.gse2034 <- independent_prog.gse2034 %>%
  mutate(
    Feature = recode_values(term,
                            "SCORE:strata(time_group)time_group=1"              ~ "Time Group 1 (>70 months)",
                            "SCORE:strata(time_group)time_group=2"              ~ "Time Group 2 (<70 months)",
                            default = term
    ),
    
    `Hazard Ratio (HR)` = round(estimate, 3),
    `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
    
    `p-value` = scales::scientific(p.value, digits = 3)
  ) %>%
  dplyr::select(Feature, `Hazard Ratio (HR)`, `95% Confidence Interval`, `p-value`) %>%
  arrange(`Hazard Ratio (HR)`)

# 3.3.3 Table

flextable(supplementary_table.gse2034) %>%
  autofit() 


# 3.4 Forest plot ignoring values that tend to infinite

cox_p_gse2034 <-  independent_prog.gse2034 %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  mutate(
      term = recode_values(
        term,
        "SCORE:strata(time_group)time_group=1"              ~ "Time Group 1 (>70 months)",
        "SCORE:strata(time_group)time_group=2"              ~ "Time Group 2 (<70 months)",
        default = term ),
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbar(orientation = "y", aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_linedraw() + 
  theme_classic(base_size = 15) +
  ggtitle("Multivariate Cox: Recurrance GSE2034") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") + 
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) +
    labs(x = "Hazard Ratio (log scale)", 
         y = "", 
         color = "Significance (p < 0.05)",
         tag = "C")
  



# 5.- Other scores --------------------------------------------------------

# 5.1 Brier score

eval_results.gse2034 <- final_fit %>%
  augment(new_data = proof_genes_pt.gse2034, eval_time = c(36, 60, 120)) 

performance.gse2034 <- eval_results.gse2034 %>%
  brier_survival(truth = surv_obj, .pred)

print(performance.gse2034)

# 5.2 Martingale and Schofeild residuals



cox.zph(cox_model.gse2034)



ggcoxzph(cox.zph(cox_model.gse2034))



ggcoxdiagnostics(cox_model.gse2034, type = "martingale",
                 linear.predictions = FALSE, ggtheme = theme_bw())
