# In this file we run the external validation for TCGA

library(timeROC)
library(survminer)

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

km_fit <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = tcga_results)

# 2.5.2 Plot

ggsurvplot(km_fit, 
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

# 2.7 Calculate the actual Concordance Index

c_index_results.tcga <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = tcga_results)



# 3.- ROC and AUC ---------------------------------------------------------


# 3.1 Area under the curve at 5 time points

res_auc_tcga <- timeROC(T = tcga_results$EVENT_MON,
                        delta = tcga_results$EVENT_STAT,
                        marker = -tcga_results$.pred_linear_pred,
                        cause = 1, # The EVENT code
                        times = c(12, 36, 60, 72, 120), # 3, 5, and 10 years
                        iid = TRUE)

# 3.1.2 View the AUC values

res_auc_res_tcga <- res_auc_tcga$AUC %>% 
  as.data.frame()


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
    HER2 = lab_proc_her2_neu_immunohistochemistry_receptor_status,
    LYMPH = as.numeric(lymph_node_examined_count),
    PAM50 = PAM50Call_RNAseq,
    AGE = as.numeric(Age_at_Initial_Pathologic_Diagnosis_nature2012),
    SCORE = tcga_results$.pred_linear_pred,
    RADIO = radiation_therapy,
    SURGERY = factor(breast_carcinoma_primary_surgical_procedure_name),
    NEO = history_of_neoadjuvant_treatment,
    OTHER_TX = additional_pharmaceutical_therapy,
    TARG_TX = targeted_molecular_therapy
  ) %>%
  dplyr::select(all_of(proof_genes),
                surv_obj,
                LYMPH,
                PAM50,
                AGE,
                SCORE,
                EVENT_STAT,
                RADIO,
                NEO,
                OTHER_TX,
                TARG_TX,
                SURGERY,
                risk_group) 

# 4.2 Actual cox model with parameters to evaluare

cox_model.tcga <- coxph(surv_obj ~ SCORE + AGE + LYMPH + PAM50,  data = proof_genes_pt.tcga.cox)

# 4.3 Tidy format

independent_prog.tcga <- cox_model.tcga %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

# 4.3.2 Results

summary(cox_model.tcga)

# 4.4 Results to paragraph 

num_param_compare <- c(2:7)

cat(paste0("The signature on TCGA got a C-score of ", round(c_index_results.tcga$concordance, 2)),
    paste0("an HR of ", round(summary_cox_tcga$coefficients[2], 2), " (CI 95% of ", round(summary_cox_tcga$conf.int[3], 2), " - ", round(summary_cox_tcga$conf.int[4], 2), " pval ", summary_cox_tcga$coefficients[5], ")"),
    paste0("As an independence factor it has an HR of ", round(independent_prog.tcga$estimate[independent_prog.tcga$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.tcga$conf.low[independent_prog.tcga$term == "SCORE"], 2), " - ", round(independent_prog.tcga$conf.high[independent_prog.tcga$term == "SCORE"], 2), " pval of ", independent_prog.tcga$p.value[independent_prog.tcga$term == "SCORE"], ")"),
    sep = ". "
)

cat(paste0(independent_prog.tcga$term[num_param_compare], " with its HR of ", round(independent_prog.tcga$estimate[num_param_compare], 2), " (CI 95% of ", round(independent_prog.tcga$conf.low[num_param_compare], 2), " - ", round(independent_prog.tcga$conf.high[num_param_compare], 2), " pval of ", independent_prog.tcga$p.value[num_param_compare], ")"),
    sep = ". "
)


independent_prog.tcga[num_param_compare, c(1, 2, 5, 6, 7)]

# 4.5 Forest plot ignoring values that tend to infinite

independent_prog.tcga %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  mutate(
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_minimal()


# 5.- Other analysis ------------------------------------------------------


# 5.1 Boxplot comparing to PAM50 independently of evnet status

ggplot(proof_genes_pt.tcga.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material")

# 5.1.2 Boxplot comparing to PAM50 facet wrapped by event status

ggplot(proof_genes_pt.tcga.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material") +
  facet_wrap( ~ EVENT_STAT)


# 5.2 Wilcoxon test where on each treatment modality or PAM50 group we compare deceased vs live patients (so if we have patients that recieved chemo we compare the score on patients that lived vs died and then on those who did not recieve chemo and so on with every treatment)

proof_genes_pt.tcga.cox %>% 
  group_by(PAM50) %>% 
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
})

# 5.3 Pivot longer so as to plot the different treatment given scores faceted by event status

proof_genes_pt.tcga_long <- 
  proof_genes_pt.tcga.cox %>% 
  pivot_longer(
    cols = c(NEO, SURGERY, RADIO, TARG_TX), # Add any other parameters here
    names_to = "Parameter",
    values_to = "Value"
  )

# 5.3.2 Plot


ggplot(data = proof_genes_pt.tcga_long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() + 
  facet_wrap(~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2) + 
  scale_fill_paletteer_d("palettesForR::Pastels") + 
  theme_gray(base_size = 18)


