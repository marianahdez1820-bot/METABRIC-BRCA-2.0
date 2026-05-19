library(survival)


label <- "For predicting recurrence in ER+ patients from METABRIC cohort"

# boruta_signature_27 <- c("S100P", "CDCA5", "PTTG1", "UHRF1", "CCNB2", "NKX2.2", "CENPF", "KIF20A", "TPX2",
#   "FAM83D", "CTTN", "RACGAP1", "CYB5D2", "TDG", "STAT5A", "SPRY4", "SYNC", "CBX7",
#   "APOBEC3B", "KIF4A", "SHMT2", "PIGV", "GTSE1", "NMRK1", "DNAJC7", "BECN1", "PEX16")
# 
# 
# final_gene_signature <- as.data.frame(read_csv("~/Documents/Servicio_Social/R.project/MyRProject/METABRIC_project/500/final_gene_signature.csv"))



# boruta_signature_49 <- c("LAD1", "UBE2C", "CDC20", "NR2F1", "CDCA5", "PTTG1",
#                          "UHRF1", "CCNB2", "NUSAP1", "ASPM", "SLC4A8", "CENPF",
#                          "LTO1", "TRIB2", "KIF20A", "SAPCD2", "PBX1", "TPX2",
#                          "TROAP", "PPFIA1", "CTTN", "EXO1", "AFMID", "INTS4",
#                          "RACGAP1", "COL4A1", "DBN1", "CKAP2L", "FOXM1", "SMC4",
#                          "DNAJB9", "STAT5A", "STIL", "NT5DC2", "MRPL27", "SPRY4",
#                          "ARPC1B", "SUMF2", "STIP1", "PPP1R14B", "FANCD2", "SLC7A4",
#                          "KIF4A", "GTF3C1", "TRIM37", "HSP90AB1", "SHMT2", "TCF3", "MXD3")

# boruta_signature_35 <- c(
#   "LAD1", "UBE2C", "CDC20", "NR2F1", "CDCA5", "PTTG1", "UHRF1", "CCNB2", "NUSAP1",
#   "ASPM", "CENPF", "LTO1", "TRIB2", "KIF20A", "TPX2", "TROAP", "CTTN", "EXO1",
#   "INTS4", "COL4A1", "DBN1", "FOXM1", "SMC4", "DNAJB9", "STAT5A", "NT5DC2",
#   "MRPL27", "SPRY4", "ARPC1B", "SUMF2", "PPP1R14B", "SLC7A4", "GTF3C1", "HSP90AB1", "SHMT2")



# common_genes_meta.gse2034 <- c("NR2F1", "TRIB2", "STAT5A", "LAD1", "TRIM37", "SPRY4", "COL4A1",
#                                "CTTN", "CCNB2", "GTF3C1", "HSP90AB1", "PPFIA1", "SLC7A4", "RACGAP1",
#                                "PBX1", "DBN1", "NUSAP1", "SLC4A8", "SMC4", "FOXM1", "DNAJB9", "ARPC1B",
#                                "PTTG1", "KIF20A", "NT5DC2", "STIP1", "PPP1R14B", "MXD3", "EXO1", "ASPM",
#                                "STIL", "TROAP", "CDC20", "UBE2C", "TPX2", "CENPF", "TCF3", "SHMT2", "KIF4A")


# 1.- Preparing metadata --------------------------------------------------

er_patients_recu <- metadata.ER_POS_REC  
ml_metadata <- er_patients_recu

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

proof_genes <- make.names(common_genes_meta.gse2034) # common_genes_meta.gse2034  #boruta_signature # #significant_genes #$term # significant_genes$term  #fused_signature #rownames(res_sig) #common_genes
# significant_genes$term %>% 
#   cat(sep = ", ")


# 1.3 Object with all the patients ER + and expression of only the genes of interest
rownames(counts_data) <- make.names(rownames(counts_data))
proof_genes_pt <- counts_data[proof_genes, er_patients_recu$PATIENT_ID]

# 1.4  Scaling is done in the linear regression recipe

proof_genes_pt <- t(proof_genes_pt) 

# 1.5.1 Check that the patients are in the same order


all(rownames(proof_genes_pt) == er_patients_recu$PATIENT_ID)
  
  

# 1.5.2 Add a column called EVENT_STAT and EVENT_MON to create the surv object
# NOTE that this file is recurrence, it is still stored in survival so as to not have to change the 
# linear regression file

proof_genes_pt <- 
  proof_genes_pt %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(er_patients_recu, by = "PATIENT_ID") %>% 
  column_to_rownames("PATIENT_ID") %>%  # Turn to factor for machine learning
  dplyr::select(all_of(proof_genes),
                EVENT_MON,
                EVENT_STAT) %>% 
  drop_na() %>% 
  mutate(surv_obj = Surv(
    time  = EVENT_MON,
    event = EVENT_STAT,
    type  = "right"))

proof_genes_pt <- proof_genes_pt %>% 
  mutate(across(-c(EVENT_MON, EVENT_STAT, surv_obj), scale))

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models
