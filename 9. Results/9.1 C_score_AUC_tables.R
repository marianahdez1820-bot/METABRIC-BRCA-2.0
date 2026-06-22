library(magrittr)


c_index_df <- bind_rows(c_index_summary, c_index_summary.tcga, c_index_summary.gse96058)

c_index_df <- 
  c_index_df %>% 
  relocate(data_set, C_Index, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE)

c_index_df <- 
  c_index_df %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 3)),
    p_value = format(p_value, scientific = TRUE, digits = 3)
  )

flextable(c_index_df) %>%
  autofit()


auc_df <- bind_rows(auc_ci, auc_ci.tcga, auc_ci.gse96058)

auc_df <- 
  auc_df %>% 
  relocate(data_set, time, AUC, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE) %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 3)),
    p_value = format(p_value, scientific = TRUE, digits = 3)
  )
flextable(auc_df) %>%
  autofit() 

summary_median.gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = gse96058_results))

summary_median.tcga <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = tcga_results))

summary_median.metabric <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = test_data))



cox_table <- data.frame(
  "Data base" = c("METABRIC", "TCGA", "GSE96058"),
  
  "HR (CI 95%, p-value )" = c(paste(summary_cox$coefficients[2], " (", summary_cox$conf.int[3], " - ", summary_cox$conf.int[4], ", p-value ", summary_cox$coefficients[5]),
                              paste(summary_cox_tcga$coefficients[2], " (", summary_cox_tcga$conf.int[3], " - ", summary_cox_tcga$conf.int[4], ", p-value ", summary_cox_tcga$coefficients[5]),
                              paste(summary_gse96058$coefficients[2], " (", summary_gse96058$conf.int[3], " - ", summary_gse96058$conf.int[4], ", p-value ", summary_gse96058$coefficients[5])),
  
  "C-Index" = c(summary_cox$concordance[1],
                summary_cox_tcga$concordance[1],
                summary_gse96058$concordance[1]),
  
  
  
  "HR (CI 95%, p-value)" = c(paste(summary_median.metabric$coefficients[2], " (", summary_median.metabric$conf.int[3], " - ", summary_median.metabric$conf.int[4], ", p-value ", summary_median.metabric$coefficients[5]),
                             paste(summary_median.tcga$coefficients[2], " (", summary_median.tcga$conf.int[3], " - ", summary_median.tcga$conf.int[4], ", p-value ", summary_median.tcga$coefficients[5]),
                             paste(summary_median.gse96058$coefficients[2], " (", summary_median.gse96058$conf.int[3], " - ", summary_median.gse96058$conf.int[4], ", p-value ", summary_median.gse96058$coefficients[5])),
  
  "C-Index " = c(summary_median.metabric$concordance[1],
                 summary_median.tcga$concordance[1],
                 summary_median.gse96058$concordance[1])
  
  
)

flextable(cox_table) %>%
  autofit()





