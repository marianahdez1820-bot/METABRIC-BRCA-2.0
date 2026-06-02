
# This file is to validate the signature obtained by lin regression on GSE96058 

library(paletteer)
library(patchwork)


# 1.- Validation Lin Reg ----------------------------------------------------------


# 1.3 Predict

gse96058_results <- predict(final_fit, new_data = proof_genes_pt.gse96058, type = "linear_pred") %>%
  bind_cols(proof_genes_pt.gse96058)

# 1.4 To get the p value

validation_test <- coxph(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, data = gse96058_results)

summary(validation_test)



# 2.- Risk stratification -------------------------------------------------


library(survminer)

# 2.1 Create risk groups based on the median of the predictions


gse96058_results <- gse96058_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < true_cut$cutpoint$cutpoint[1], "High Risk", "Low Risk")))


# 2.2 Fit the KM curve

km_fit.gse96058 <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = gse96058_results)

# 2.3 Plot


ggsurvplot(km_fit.gse96058, 
           data = gse96058_results, 
           pval = TRUE, 
           risk.table = TRUE,
           
           title = "Validation in GSE96058",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           ylim = c(0.7, 1),
           xlim = c(0, 85),         # Zoom in
           break.time.by = 10,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#E41A1C", "#377EB8"), # Custom color palette
           
)


gse96058_results <- gse96058_results %>%
  mutate(pred_z = scale(.pred_linear_pred))


# 2.4 Run Cox again

gse96058_results$risk_group <- relevel(gse96058_results$risk_group, ref = "Low Risk")

summary_gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = gse96058_results))

# 2.5 Calculate the actual Concordance Index

c_index_results.gse96058 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                        data = gse96058_results)


# 2.6 Table with confidence interval and z stat and estimated p val

c_index_summary.gse96058 <- data.frame(
  C_Index  = c_index_results.gse96058$concordance,
  SE       = sqrt(c_index_results.gse96058$var)
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

print(c_index_summary.gse96058)


library(timeROC)

# 2.6 Area under the curve per time

res_auc_gse96058 <- timeROC(T = gse96058_results$EVENT_MON,
                   delta = gse96058_results$EVENT_STAT,
                   marker = -gse96058_results$.pred_linear_pred,
                   cause = 1, # The event code
                   times = c(12, 36, 60, 72), # 3, 5, and 6 years
                   iid = TRUE)

# 2.6.2 Table with confidence interval and z stat and estimated p val

auc_ci.gse96058 <- data.frame(
  AUC  = res_auc_gse96058$AUC,
  SE   = res_auc_gse96058$inference$vect_sd_1,
  time = res_auc_gse96058$times
) %>% 
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

# 2.6.3 Text

for (i in 1:5) {
  
  cat(paste0(round(auc_ci.gse96058$AUC[i], 2), " (", round(auc_ci.gse96058$conf_int_low95[i], 2), "-", round(auc_ci.gse96058$conf_int_high95[i], 2), "z score ", round(auc_ci.gse96058$z_stat[i], 2),  " pval ", round(auc_ci.gse96058$p_value[i], 4), ")", " at ", auc_ci.gse96058$time[i], " months (",auc_ci.gse96058$time[i] / 12, " years), "))
  
}


# 2.7 Combine all time points into a long dataframe

plot_roc.gse96058 <- map_df(c(12, 36, 60, 72), function(i) { # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc_gse96058$FP[, time_label],
    TP = res_auc_gse96058$TP[, time_label],
    Time = factor(i),
    data_set = "GSE96058"
  )
})


# 2.7.5 Global ROC object


gse96058_results <- 
  gse96058_results %>% 
  mutate(EVENT_STAT = factor(EVENT_STAT))

global_roc.gse96058 <- roc_curve(gse96058_results,
                        EVENT_STAT,
                        .pred_linear_pred
                        ) %>%
  mutate(label = "GSE96058")

# 2.8 Extract AUC and SE for each time point

auc_table.gse96058 <- data.frame(
  Time = c(12, 36, 60, 72),
  `AUC (%)` = round(as.numeric(res_auc_gse96058$AUC[c("t=36", "t=60", "t=72", "t=80")]), 4),
  `SE` = res_auc_gse96058$times
)


# 2.8.2 Labels

legend_labels.gse96058 <- paste0("t=", auc_table.gse96058$Time, " (AUC: ", 100 * auc_table.gse96058$AUC..., "%)")

# 2.8.3 Labels for faceted plot

facet_labels.gse96058 <- data.frame(
  Time = factor(c(12, 36, 60, 72)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(res_auc_gse96058$AUC[1:4]), 3), "%"),
  data_set = "GSE96058"
)

# 2.9 Plot with overlap

ggplot(plot_roc.gse96058, aes(x = FP, y = TP, color = Time)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(
    title = "Time-Dependent ROC Curves",
    subtitle = "Comparing model performance across different horizons",
    x = "1 - Specificity (FP)",
    y = "Sensitivity (TP)",
    color = "Time Point"
  ) +
  scale_color_viridis_d(labels = legend_labels.gse96058) 

# 2.9.2 Plot with facet

ggplot(plot_roc.gse96058, aes(x = FP, y = TP)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Time, ncol = 2) +
  theme_bw() +
  labs(
    title = "ROC Curves by Time Point",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  geom_text(data = facet_labels.gse96058, 
            aes(x = 0.75, y = 0.1, label = AUC_Text), 
            size = 4, fontface = "bold")


# 3.- Independence test ---------------------------------------------------

# 3.1 Prepare data for multivariate cox

proof_genes_pt_gse96058.cox <- 
  proof_genes_pt.gse96058 %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title", suffix = c("", ".y")) %>% # Join with metadata
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("title") %>% 
  mutate(HER2 = her2_pred_sgc, 
         LYMPH = lymph_group,
         PAM50 = pam50,
         AGE = as.numeric(age),
         KI67 = ki67_pred_sgc,
         SCORE = gse96058_results$.pred_linear_pred,
         RISK = gse96058_results$risk_group,
         HORMONE = endocrine_tx,
         CHEMO = chemo_tx
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                LYMPH,
                PAM50,
                AGE,
                HER2,
                KI67,
                SCORE,
                RISK,
                EVENT_STAT,
                EVENT_MON,
                HORMONE,
                CHEMO
  ) %>% 
  na.omit()

summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, 
              data = proof_genes_pt_gse96058.cox))

# 3.2 Multivariate cox

cox_model.gse96058 <- coxph(surv_obj ~ PAM50 + KI67 + HER2 + AGE + LYMPH + CHEMO + HORMONE+ SCORE, 
                            data = proof_genes_pt_gse96058.cox)

independent_prog.gse96058 <- cox_model.gse96058 %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

summary(cox_model.gse96058)


# 4.- Results -------------------------------------------------------------

# 4.1 Index numbers to use of the multivariate cox

num_param_compare_gse <- c(2, 3, 5, 7, 8,10, 11, 12, 13, 14)

# 4.2 Concatenate text with desired results

cat(paste0("The signature on GSE96058 got a C-score of ", round(c_index_results.gse96058$concordance, 2)),
    paste0("an HR of ", round(summary_gse96058$coefficients[2], 2), " (CI 95% of ", round(summary_gse96058$conf.int[3], 2), " - ", round(summary_gse96058$conf.int[4], 2), " pval ", summary_gse96058$coefficients[5], ")"),
    paste0("As an independence factor it has an HR of ", round(independent_prog.gse96058$estimate[independent_prog.gse96058$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.gse96058$conf.low[independent_prog.gse96058$term == "SCORE"], 2), " - ", round(independent_prog.gse96058$conf.high[independent_prog.gse96058$term == "SCORE"], 2), " pval of ", independent_prog.gse96058$p.value[independent_prog.gse96058$term == "SCORE"], ")"),
    sep = ". "
)

cat(paste0(independent_prog.gse96058$term[num_param_compare_gse], " with its HR of ", round(independent_prog.gse96058$estimate[num_param_compare_gse], 2), " (CI 95% of ", round(independent_prog.gse96058$conf.low[num_param_compare_gse], 2), " - ", round(independent_prog.gse96058$conf.high[num_param_compare_gse], 2), " pval of ", independent_prog.gse96058$p.value[num_param_compare_gse], ")"),
    sep = ". "
)

# 4.3 Forest Plot

cox_p_gse96058 <- independent_prog.gse96058[num_param_compare_gse, ] %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  filter(!(term == "LYMPHNA")) %>% 
  mutate(
    term = recode_values(
      term,
      "SCORE"              ~ "Signature Score",
      "LYMPH4yoX"              ~ "Lymph nodes >= 4",
      "LYMPHNodeNegative"              ~ "Lymph nodes negative",
      "LYMPHSubMicroMet"              ~ "Lymph micrometastasis",
      "AGE"                ~ "Age at Diagnosis",
      "CHEMO"           ~ "Chemotherapy",
      "HORMONE"         ~ "Hormone Therapy",
      "PAM50LumA"          ~ "Claudin subtype Luminal A",
      "PAM50LumB"          ~ "Claudin subtype Luminal B",
      "PAM50Her2"          ~ "Claudin subtype Her2-enriched",
      "PAM50Normal"        ~ "Claudin subtype Normal-like",
      default = term 
    ),
    term = reorder(term, estimate),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", 
       y = "", 
       color = "Significance (p < 0.05)",
       tag = "C"
       ) +
  theme_classic(base_size = 22) +
  ggtitle("Multivariate Cox: Survival GSE96058") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",) + 
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) 

# 4.4 Boxplot comparing to PAM50


score_subtype_gse96058 <- ggplot(proof_genes_pt_gse96058.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  facet_wrap(~ EVENT_STAT, 
             labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) + 
  theme_classic(base_size = 22) + 
  ggtitle("Score change on subtype: GSE96058") +
  theme(
    strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    y = "", 
    x = "Intrinsic Subtype",
    tag = "C"   ) + 
  scale_fill_paletteer_d("Redmonder::dPBIPuOr") 
  


# 4.5 Wilcoxon test for PAM50


proof_genes_pt_gse96058.cox %>%
  dplyr::select(EVENT_STAT, HORMONE, SCORE) %>% 
  group_by(HORMONE) %>%
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "holm"))



# 4.6 Kruskal wallis

kruskal.test(SCORE ~  CHEMO, data = proof_genes_pt_gse96058.cox)

# 4.6.2 Long format to visualize how the scorescores patient based on event an treatment

proof_genes_long.gse96058 <- proof_genes_pt_gse96058.cox %>%
  mutate(HORMONE = factor(HORMONE),
         CHEMO = factor(CHEMO)) %>% 
  pivot_longer(
    cols = c(HORMONE, CHEMO), 
    names_to = "Parameter",
    values_to = "Value"
  )


# 4.6.3  Plot

score_tx_gse96058 <- ggplot(proof_genes_long.gse96058,
       aes(
         y = SCORE,
         x = Value,
         fill = factor(Value, labels = c("0" = "Untreated", "1" = "Treated"))
       )) +
  geom_boxplot() +
  facet_wrap(
    ~ Parameter + EVENT_STAT,
    scales = "free_x",
    labeller = labeller(EVENT_STAT = as_labeller(c(
      "0" = "Alive", "1" = "Deceased")
    ),
    Parameter  = as_labeller(c("CHEMO" = "Chemotherapy", "HORMONE" = "Hormonal treatment"))
    )
  ) +
  scale_fill_paletteer_d("khroma::iridescent", direction = - 1) + 
  theme_classic(base_size = 15) + 
  ggtitle("Score change on treatment: GSE96058") +
  theme(
    strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  labs(tag = "C",
       y = "",
       x  = "")

summary(coxph(surv_obj ~ SCORE * CHEMO, data = proof_genes_pt_gse96058.cox))



# 6.- Other evaluations ---------------------------------------------------

# 6.1 Brier score

eval_results.gse96058 <- final_fit %>%
  augment(new_data = proof_genes_pt.gse96058, eval_time = c(36, 60, 120)) 

# 6.2 Calculate Brier Score

performance.gse96058 <- eval_results.gse96058 %>%
  brier_survival(truth = surv_obj, .pred)

print(performance.gse96058)

# 6.3 Martingale and Schofeild residuals

cox.zph(cox_model.gse96058)

ggcoxzph(cox.zph(cox_model.gse96058))

ggcoxdiagnostics(cox_model.gse96058, type = "martingale",
                 linear.predictions = FALSE, ggtheme = theme_bw())
