# In this file we prepare the data for the linear regression model using only ER+ patients
# and with recurrence parameters
boruta_signature <- read.csv("C:/R/METABRIC/Results/final_gene_signature.csv", 
                             check.names = FALSE)
boruta_signature <- as.character(boruta_signature$x)

label <- "For predicting recurrence in ER+ patients from METABRIC cohort"

# 1.- Preparing metadata --------------------------------------------------

er_patients_recu <- metadata %>% 
  filter(ER_IHC == "Positve")

ml_metadata <- er_patients_recu

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

proof_genes <- boruta_signature # common_genes_meta.gse2043  #boruta_signature # #significant_genes #$term # significant_genes$term  #fused_signature #rownames(res_sig) #common_genes
# significant_genes$term %>% 
#   cat(sep = ", ")


# 1.3 Object with all the patients ER + and expression of only the genes of interest
rownames(counts_data) <- make.names(rownames(counts_data))
proof_genes_pt <- counts_data[proof_genes, er_patients_recu$PATIENT_ID]

# 1.4  Scaling is done in the linear regression recipe


proof_genes_pt <- t(proof_genes_pt) 

# 1.5.1 Check that the patients are in the same order


all(rownames(proof_genes_pt) == er_patients_recu$PATIENT_ID)

# 1.5.2 Add a column called SURVIVAL and SURVIVAL_MON to create the surv object
# NOTE that this file is recurrence, it is still stored in survival so as to not have to change the 
# linear regression file

proof_genes_pt <- 
  proof_genes_pt %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(er_patients_recu, by = "PATIENT_ID") %>% 
  column_to_rownames("PATIENT_ID") %>% 
  mutate(EVENT_STAT = as.numeric(RECURR_STAT),
         EVENT_MON = as.numeric(RFS_MONTHS)
  ) %>%  # Turn to factor for machine learning
  dplyr::select(all_of(proof_genes),
                EVENT_MON,
                EVENT_STAT) %>% 
  filter(EVENT_MON > 0) %>% # Eliminate those with 0 survival months
  drop_na() %>% 
  mutate(surv_obj = Surv(
    time  = EVENT_MON,
    event = EVENT_STAT,
    type  = "right"))

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models
