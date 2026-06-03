# In this file we run the external validation for TCGA

library(timeROC)
library(survminer)
library(patchwork)
library(coxphf)

# 1.- EXTERNAL VALIDATION ON TCGA

# 1.1 Predict with the parameters from the fit to the TCGA data



tcga_results <- predict(final_fit, new_data = proof_genes_pt.tcga, type = "linear_pred") %>%
  bind_cols((proof_genes_pt.tcga %>% 
               rownames_to_column("sampleID")))



# 2.- Divide by risk groups -----------------------------------------------


# 2.3 Create risk groups based on the median of the predictions or the cutpoint

tcga_results <-
  tcga_results %>%
  mutate(risk_group = factor(ifelse(.pred_linear_pred < true_cut$cutpoint$cutpoint[1], "High Risk", "Low Risk")))

# 2.3.2 Relevel so as to have low risk as reference

tcga_results$risk_group <- relevel(tcga_results$risk_group, ref = "Low Risk")

# 2.4 Fit the KM curve

km_fit.tcga <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = tcga_results)

# 2.5.2 Plot

ggsurvplot(km_fit.tcga, 
           data = tcga_results, 
           pval = TRUE, 
           risk.table = TRUE,
           
           title = "Validation in TCGA",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 300),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#E7B800", "#2E9FDF"), # Custom color palette
           
)



# 2.6 Run Cox again

summary_cox_tcga <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = tcga_results))

# 2.7 Calculate the Concordance Index

c_index_results.tcga <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = tcga_results)



# 2.8 Table with confidence interval and z stat and estimated p val

c_index_summary.tcga <- data.frame(
  C_Index = c_index_results.tcga$concordance,
  SE = sqrt(c_index_results.tcga$var)
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

print(c_index_summary.tcga)

# 3.- ROC and AUC ---------------------------------------------------------


# 3.1 Area under the curve at 5 time points

res_auc_tcga <- timeROC(T = tcga_results$EVENT_MON,
                        delta = tcga_results$EVENT_STAT,
                        marker = -tcga_results$.pred_linear_pred,
                        cause = 1, # The EVENT code
                        times = c(12, 36, 60, 72, 120), # 3, 5, and 10 years
                        iid = TRUE)

# 3.1.1.2 Table with confidence interval and z stat and estimated p val

auc_ci.tcga <- data.frame(
  AUC  = res_auc_tcga$AUC,
  SE   = res_auc_tcga$inference$vect_sd_1,
  time = res_auc_tcga$times
) %>% 
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
      )

# 3.1.1.3 Text

for (i in 1:5) {
  
 cat(paste0(round(auc_ci.tcga$AUC[i], 2), " (", round(auc_ci.tcga$conf_int_low95[i], 2), "-", round(auc_ci.tcga$conf_int_high95[i], 2), "z score ", round(auc_ci.tcga$z_stat[i], 2),  " pval ", round(auc_ci.tcga$p_value[i], 4), ")", " at ", auc_ci.tcga$time[i], " months (",auc_ci.tcga$time[i] / 12, " years), "))

}


# 3.1.2 View the AUC values

res_auc_res_tcga <- res_auc_tcga$AUC %>% 
  as.data.frame()

# 3.1.3 

tcga_results <- 
  tcga_results %>% 
  mutate(EVENT_STAT = factor(EVENT_STAT))

# 3.1.4 Global ROC curve object

global_roc.tcga <- roc_curve(tcga_results,
                                 EVENT_STAT,
                                 .pred_linear_pred
                             ) %>%
  mutate(label = "TCGA")


# 3.2 Combine all time points into a long dataframe

plot_roc.tcga <- map_df(c(12, 36, 60, 72, 120), function(i) { # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc_tcga$FP[, time_label],
    TP = res_auc_tcga$TP[, time_label],
    Time = factor(i),
    data_set = "TCGA"
  )
})


# 3.3 Labels for faceted plot

facet_labels.tcga <- data.frame(
  Time = factor(c(12, 36, 60, 72, 120)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(res_auc_tcga$AUC[1:5]), 3), "%"),
  data_set = "TCGA"
)

# 3.4 Plot with facet

ggplot(plot_roc.tcga, aes(x = FP, y = TP)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Time, ncol = 2) +
  theme_bw() +
  labs(
    title = "ROC Curves by Time Point",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  geom_text(data = facet_labels.tcga, 
            aes(x = 0.75, y = 0.1, label = AUC_Text), 
            size = 4, fontface = "bold")



# 4.- Multivariate cox ----------------------------------------------------


# 4.1 Prepare objects

# 4.1.1 First merge with score to then merge with the metadata

proof_genes_pt.tcga <- 
  proof_genes_pt.tcga %>% 
  rownames_to_column("sampleID") %>% 
  left_join(tcga_results, by = "sampleID", suffix = c("", ".drop")) %>%
  dplyr::select(-ends_with(".drop"))

# 4.1.2 Merge with metadata and select parameters to evaluate

proof_genes_pt.tcga.cox <-
  proof_genes_pt.tcga %>%
  left_join(refined_data_unique, by = "sampleID") %>%
  as.data.frame() %>%
  column_to_rownames("sampleID") %>%
  mutate(
    SCORE = tcga_results$.pred_linear_pred,
    ) %>%
  dplyr::select(all_of(proof_genes),
                surv_obj,
                LYMPH,
                PAM50,
                AGE,
                SCORE,
                EVENT_STAT,
                EVENT_MON,
                RADIO,
                NEO,
                OTHER_TX,
                TARG_TX,
                SURGERY,
                risk_group,
                HER2,
                HIST,
                MENO,
                INTCLUST)



# 4.2 Actual cox model with parameters to evaluare

cox_model.tcga <- coxph(surv_obj ~ SCORE + AGE + LYMPH + RADIO + NEO + TARG_TX + strata(PAM50) + strata(HER2) + strata(HIST), data = proof_genes_pt.tcga.cox)

summary(coxph(surv_obj ~ AGE + LYMPH + SCORE  ,  data = proof_genes_pt.tcga.cox))

# 4.3 Tidy format

independent_prog.tcga <- cox_model.tcga %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

# 4.3.2 Results

summary(cox_model.tcga)

# 4.4 Results to paragraph 

num_param_compare <- c(2:9)

cat(paste0("The signature on TCGA got a C-score of ", round(c_index_results.tcga$concordance, 2)),
    paste("AUC at 1 year of ", res_auc_res_tcga$.[1], ", at 3 years of ", res_auc_res_tcga$.[2], ", at 5 years of ", res_auc_res_tcga$.[3], ", at 6 years of ", res_auc_res_tcga$.[4], ", at 10 years of ", res_auc_res_tcga$.[5]),
    paste0("an HR of ", round(summary_cox_tcga$coefficients[2], 2), " (CI 95% of ", round(summary_cox_tcga$conf.int[3], 2), " - ", round(summary_cox_tcga$conf.int[4], 2), " pval ", summary_cox_tcga$coefficients[5], ")"),
    paste0("As an independence factor it has an HR of ", round(independent_prog.tcga$estimate[independent_prog.tcga$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.tcga$conf.low[independent_prog.tcga$term == "SCORE"], 2), " - ", round(independent_prog.tcga$conf.high[independent_prog.tcga$term == "SCORE"], 2), " pval of ", independent_prog.tcga$p.value[independent_prog.tcga$term == "SCORE"], ")"),
    sep = ". "
)

cat(paste0(independent_prog.tcga$term[num_param_compare], " with its HR of ", round(independent_prog.tcga$estimate[num_param_compare], 2), " (CI 95% of ", round(independent_prog.tcga$conf.low[num_param_compare], 2), " - ", round(independent_prog.tcga$conf.high[num_param_compare], 2), " pval of ", independent_prog.tcga$p.value[num_param_compare], ")"),
    sep = ". "
)



# 4.5 Forest plot ignoring values that tend to infinite

cox_p_tcga <- independent_prog.tcga %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  mutate(
    
    term = recode_values(
      term,
      "SCORE"              ~ "Signature Score",
      "LYMPH"              ~ "Lymph Node Status",
      "AGE"                ~ "Age at Diagnosis",
      default = term 
    ),
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic(base_size = 22) +
  ggtitle("Multivariate Cox: Survival TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5), 
        legend.box.background = element_rect(color = "black", linewidth = 1),
        legend.key = element_rect(color = "gray80", linewidth = 0.5)) + 
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) + 
  labs(x = "Hazard Ratio (log scale)", 
       y = "Clinical & Molecular Features", 
       color = "Significance (p < 0.05)",
       tag = "B")


# 5.- Other analysis ------------------------------------------------------


# 5.1.2 Boxplot comparing to PAM50 facet wrapped by event status


score_subtype_tcga <- ggplot(proof_genes_pt.tcga.cox %>% drop_na(PAM50), aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  facet_wrap(~ EVENT_STAT,
             scales = "free_x",
             labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) + 
  theme_classic(base_size = 22) + 
  ggtitle("Score change on subtype: TCGA") +
  theme(
    strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  labs(y = "Score", 
       x = "",
       tag = "B") + 
  scale_fill_paletteer_d("Redmonder::dPBIPuOr") 


# 5.2 Wilcoxon test where on each treatment modality or PAM50 group we compare deceased vs live patients (so if we have patients that recieved chemo we compare the score on patients that lived vs died and then on those who did not recieve chemo and so on with every treatment)

proof_genes_pt.tcga.cox %>% 
  group_by(PAM50, EVENT_STAT) %>% 
  filter(!is.na(PAM50)) %>% 
  group_by(PAM50) %>% 
  filter(n_distinct(EVENT_STAT) == 2) %>% 
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "holm"))


# 5.2.2 Similar but for treatment

vars <- c("RADIO", "NEO", "OTHER_TX", "TARG_TX")

# To make data frame do the "loop" with map_dfr

results <- map_dfr(vars, function(var) {
  
  df <- proof_genes_pt.tcga.cox %>%
    group_by(.data[[var]]) %>%
    filter(n_distinct(EVENT_STAT) == 2) %>% 
    group_modify(~ {
      test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
      tidy(test)
    }) %>%
    mutate(level = cur_group()[[1]])  %>%
    ungroup() %>%
    mutate(
      variable = var,
      adj_p_value = p.adjust(p.value, method = "holm")
    ) %>%
    dplyr::select(variable, level, p.value, adj_p_value)
}
)

results

# 5.3 Pivot longer so as to plot the different treatment given scores faceted by event status

proof_genes_pt.tcga_long <- 
  proof_genes_pt.tcga.cox %>% 
  pivot_longer(
    cols = c(NEO, SURGERY, RADIO, TARG_TX), # Add any other parameters here
    names_to = "Parameter",
    values_to = "Value"
  ) %>% 
  mutate(
    Value = case_when(
      toupper(Value) == "YES" ~ "Yes",
      toupper(Value) == "NO"  ~ "No",
      TRUE ~ Value # Keeps things like "Lumpectomy" and "NA" untouched
    ),
    Value = case_when(
      tolower(Value) == "mastectomy nos" ~ "Mastectomy",
      tolower(Value) == "modified radical mastectomy" ~ "Modified Radical Mastectomy",
      TRUE ~ str_to_title(Value) # Capitalizes first letter of other categories like "Lumpectomy", "Other"
    )
  )

# 5.3.2 Plot


score_tx_tcga <- ggplot(data = proof_genes_pt.tcga_long %>% drop_na(Value), aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() +
  facet_wrap( ~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2, labeller = labeller(
    EVENT_STAT = as_labeller(c(
      "0" = "Alive", "1" = "Deceased")
    ),
    Parameter  = as_labeller(c("NEO" = "Neoadjuvant therapy", "RADIO" = "Radiotherapy", "SURGERY" = "Surgery modality", "TARG_TX" = "Targeted treatment"))
  )) +
  scale_fill_paletteer_d("khroma::iridescent", direction = - 1,
                         labels = function(x) str_to_title(x)) + 
  theme_classic(base_size = 15) +
  labs(x = "Treatment modality", 
       y = "",
       fill = "Treated",
       tag = "B") +
  scale_x_discrete(labels = c("0" = "Untreated", "1" = "Treated")) +
  ggtitle("Score change on treatment: TCGA") +
  theme(
    strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 22, hjust = 1, vjust = 1)
  )

# 5.4 Cox based on treatment

for (i in vars) {
  formula <- as.formula(paste("surv_obj ~", i, "* SCORE"))
  print(summary(coxph(formula, data = proof_genes_pt.tcga.cox)))
  
  print(i)
}


# 7.- Other scores --------------------------------------------------------

# 7.1 Brier score

eval_results.tcga <- final_fit %>%
  augment(new_data = proof_genes_pt.tcga, eval_time = c(36, 60, 120)) 

performance.tcga <- eval_results.tcga %>%
  brier_survival(truth = surv_obj, .pred)

print(performance.tcga)

# 7.2 Martingale and Schofeild residuals

cox.zph(cox_model.tcga)



ggcoxzph(cox.zph(cox_model.tcga))



ggcoxdiagnostics(cox_model.tcga, type = "martingale",
                linear.predictions = FALSE, ggtheme = theme_bw())

