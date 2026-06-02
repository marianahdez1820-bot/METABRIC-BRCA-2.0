library(survival)

# In this file we prepare the data for the linear regression model using onl ER+ patients
# and with survival parameters

label <- " For predicting survival in ER positive patients from METABIC cohort"

# 1.- Preparing metadata --------------------------------------------------

er_patients_surv <- metadata.ER_POS_SURV 

ml_metadata <- er_patients_surv

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

proof_genes <- make.names(common_genes_meta.gse96058) 

# 1.3 Object with ER+ patients and expression of only the genes of interest
rownames(counts_data) <- make.names(rownames(counts_data))

proof_genes_pt <- 
  counts_data[proof_genes, er_patients_surv$PATIENT_ID] %>%  # pt means patients
  t()


# 1.4.1 Check that the patients are in the same order

all(rownames(proof_genes_pt) == er_patients_surv$PATIENT_ID)

# 1.4.2 Add a column of EVENT as a binary term for it to be the outcome and the months of survival

proof_genes_pt <- 
  proof_genes_pt %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(er_patients_surv, by = "PATIENT_ID") %>% 
  column_to_rownames("PATIENT_ID") %>%  # Turn to factor for machine learning
  dplyr::select(all_of(proof_genes),
                EVENT_MON,
                EVENT_STAT) %>% 
  filter(EVENT_MON > 0) %>% # Eliminate those with 0 survival months
  drop_na() %>% 
  mutate(surv_obj = Surv(
    time  = EVENT_MON,
    event = EVENT_STAT,
    type  = "right"))



proof_genes_pt <- 
  proof_genes_pt %>% 
  mutate(across(- c(EVENT_MON, EVENT_STAT, surv_obj), scale))

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 proof_genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models

