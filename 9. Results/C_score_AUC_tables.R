
c_index_df <- bind_rows(c_index_summary.tcga, c_index_summary.gse96058, c_index_summary)

c_index_df <- 
  c_index_df %>% 
  relocate(data_set, C_Index, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE)

c_index_df <- 
  c_index_df %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 2)),
    across(p_value, \(x) round(x, digits = 12))
  )

flextable(c_index_df) %>%
  autofit() 


auc_df <- bind_rows(auc_ci.tcga, auc_ci.gse96058, auc_ci)

auc_df <- 
  auc_df %>% 
  relocate(data_set, time, AUC, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE) %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 2)),
    across(p_value, \(x) round(x, digits = 12))
  )

flextable(auc_df) %>%
  autofit() 
