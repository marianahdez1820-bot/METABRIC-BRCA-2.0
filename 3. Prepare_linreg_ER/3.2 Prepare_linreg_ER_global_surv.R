library(survival)


# In this file we prepare the data for the linear regression model using onl ER+ patients
# and with survival parameters

#boruta_signature <- c("SLC7A2", "CXCL14", "LAD1", "UBE2C", "FCER1A", "INAVA", "COL17A1", "CBX2", "TMEM26", "CDC20", "KIF1A", "ELOVL5", "CDCA5", "PTTG1", "UHRF1", "ZIC2", "WDR72", "BIRC5", "PRC1", "NUSAP1", "RPL26", "AURKA", "TMEM132A", "NKX2.2", "CDC45", "SLC4A8", "CENPF", "ADGRG1", "KIF20A", "ADIPOR2", "CD1E", "PXDNL", "TPX2", "TROAP", "CACNG1", "FAM83D", "EXO1", "ZNF148", "EZR", "SOSTDC1", "FEN1", "COL4A1", "RACGAP1", "KIF2C", "CCNA2", "FOXM1", "CNST", "CKAP2L", "DBN1", "STAT5A", "TDG", "PSMD3", "RRBP1", "RNF24", "TCL1B", "STIP1", "HHEX", "GSK3B", "YBEY", "ATP2A2", "CNIH2", "KIF4A", "ZC3H11A", "TRIM4", "EXT1", "UBE2O", "TBC1D31", "ZMIZ1", "ENO1", "CUL5", "SLC2A4RG", "SHMT2", "AK3", "MICU2", "DIRAS3", "GSTK1", "ACADS", "FLT3", "CCT6B", "C14orf184", "PIGV", "KIF14", "TMEM60", "GTSE1", "SPTBN2", "PSMD2", "GAS2L3", "SPOPL", "SEC62", "MAP2K7", "LARP1", "HMGCR", "NT5M", "HYI", "CNTF", "C11orf71")

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

