library(survival)

# In this file we prepare the data for the linear regression model using onl ER+ patients
# and with recurrence parameters

label <- "For predicting recurrence in ER+ patients from METABRIC cohort"


# 1.- Preparing metadata --------------------------------------------------

er_patients_recu <- metadata.ER_POS_REC  
ml_metadata <- er_patients_recu

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

proof_genes <- make.names(boruta_signature) 


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

train_data <- 
  proof_genes_pt[rownames(proof_genes_pt) %in% train_rec.id, ]

test_data <- 
  proof_genes_pt[rownames(proof_genes_pt) %in% test_rec.id, ] 



train_data <- train_data %>%
  mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))

test_data <- test_data %>%
  mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))

proof_genes_pt <- proof_genes_pt %>%
  mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models
