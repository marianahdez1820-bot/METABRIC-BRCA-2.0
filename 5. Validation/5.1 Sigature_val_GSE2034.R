# This script is to perform signature validation on GSE 2034

library(oligo)
library(GEOquery)
library(tidyverse)
library(limma)
library(hgu133plus2.db)
library(gridExtra)
library(readr)



# 1.- Load data and metadata_gse.2034  ----------------------------------------------

# 1.1 Get supplementary files

# getGEOSuppFiles("GSE2034", baseDir = "D:/tcga")

# 1.1.2 Untar files
#untar("D:/tcga/GSE2034/GSE2034_RAW.tar", exdir = "D:/tcga/GSE2034/")

# 1.1.3 Listing .cel files

cel_files <- list.celfiles("D:/tcga/GSE2034/", full.names = TRUE, listGzipped = TRUE)

# 1.1.4 Reading in cel files

pre_raw_data <- read.celfiles(cel_files)

raw_data <- pre_raw_data

# 1.2 metadata_GSE2034 

pre_metadata <- read.delim("C:/R/METABRIC/Validation/GSE2034/Metadatos/GSE2034_metadata", sep = ",")
metadata_gse.2034 <- pre_metadata

#> pData of raw_data initially only contains the id for the counts and an index
#> meanwhile metadata_gse.2034  contains the full metadata but the ids differ from the count data ids

# Object that specifies which GSE is in use

gse_obj <- "GSE2034"



# 3.- Preprocessing metadata_gse.2034  --------------------------------------------------

#3.1 Clean metadata_gse.2034 

metadata_gse.2034  <- metadata_gse.2034  %>%
  mutate(
    EVENT_STAT = ifelse(relapse..1.True. == 1, 1, 0), # Object to evaluate event
    EVENT_STAT = as.numeric( EVENT_STAT),
    EVENT_MON = as.numeric(time.to.relapse.or.last.follow.up..months.), # Object to evaluate time to event
    id = GEO.asscession.number,
    ER_STAT = ER.Status,
    BRAIN_REL = Brain.relapses..1.yes..0.no.
  ) %>% 
  dplyr::select(EVENT_STAT,
                EVENT_MON,
                ER_STAT,
                id,
                BRAIN_REL,
                
  )

# 3.2 Object with the names of each file

id <- sampleNames(raw_data)
sample_names <- gsub("\\..*", "", id) # Raking out the .CEL.gz

# 3.2.2 Add the object with the general names to the phenotype data

pData(raw_data)$id <- sample_names

# 3.4 Join metadata_gse.2034  

pData(raw_data) <- 
  pData(raw_data) %>% 
  rownames_to_column("file_name") %>% # To have the full names
  full_join(metadata_gse.2034 , by = "id", keep = FALSE) %>% # Join by names without .CEL.gz
  mutate(comp_file_name = file_name) %>% # Create new column to then add to rownames
  mutate(file_name = comp_file_name) %>% 
  column_to_rownames("comp_file_name") 


metadata_gse.2034  <- pData(raw_data) # Make shure they have the same data this so that metadata also has the identifiers with .CEL.gz to match with the counts


# 4.- Preprocess data -----------------------------------------------------


# Normalize

norm_data <- rma(raw_data) 

# Boxplot after normalization

boxplot(exprs(norm_data), 
        las = 2, 
        main = paste0("RMA normalized - ", gse_obj))

# 4.3 Expression matrix

expr_matrix <- exprs(norm_data)

# 5.- Probe id to symbol --------------------------------------------------

# 5.1 Get mapping to change probe ids to gene names

library(hgu133a.db)

probe_gene <- AnnotationDbi::select(
  hgu133a.db, # The probe identifiers for this affymetrix
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

# 5.2 Delete NA symbols and duplicate symbols and add to rownames

# 5.2.1 Transform matrix to a tidy data frame and add Symbols

gene_expres_matrix <- 
  expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  inner_join(probe_gene, by = "PROBEID") %>%  # Join with the mapping object
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%   # Remove NAs and empty symbols
  mutate(variance = apply(dplyr::select(., -PROBEID, -SYMBOL), 1, var)) %>%   # 5.2.2 Calculate variance for each probe across all samples
  group_by(SYMBOL) %>% 
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # 5.2.3 Keep only the probe with the highest variance per Gene Symbol
  ungroup() %>%
  dplyr::select(-PROBEID, -variance) %>%  # 5.2.4 Remove columns, add symbol to rownames and reformat to matrix
  column_to_rownames("SYMBOL") %>%
  as.matrix()


# 5.3.3 Convert to a data frame so to add clinical columns

proof_genes_pt_gse2034 <- as.data.frame(gene_expres_matrix)

# 5.4 Keep only ER+ patients

metadata.gse.2034_er_pos <- metadata_gse.2034  %>% 
  filter(ER_STAT == "ER+")

# 5.5 Keep oly patients that are found on the counts

proof_genes_pt_gse2034 <- proof_genes_pt_gse2034[,colnames(proof_genes_pt_gse2034) %in% metadata.gse.2034_er_pos$file_name] %>% 
  as.data.frame() %>% 
  t() 

library(survival)

# Keep only genes to analyze and join with the metadata.gse.2034_er_pos that was cleaned earlier

proof_genes_pt_gse2034 <- 
  proof_genes_pt_gse2034[,proof_genes] %>% # Keep only the genes to test
  as.data.frame() %>% 
  rownames_to_column("file_name") %>%
  left_join(metadata.gse.2034_er_pos, by = "file_name") %>% # Join with metadata
  mutate(
    surv_obj = Surv(time = EVENT_MON, event =  EVENT_STAT) # Create the Survival Object inside the dataframe
  ) %>% 
  column_to_rownames("file_name") %>% 
  dplyr::select(- index, # Select only the columns of interest and discard the rest
                - id,
                - ER_STAT,
                - BRAIN_REL)


# Find genes present in both data sets
common_genes_meta.gse2034 <- intersect(proof_genes, colnames(proof_genes_pt_gse2034))

# Check how many are lost
print(paste("Original:", length(proof_genes), "Common:", length(common_genes_meta.gse2034)))

# If there are gene missing, run the lin reg with common_genes as the gene list



##################################################################

# 6.- Validation ----------------------------------------------------------

# 6.1 Extract the trained recipe from the workflow

trained_rec <- extract_recipe(final_fit)

# 6.2 "Bake" the RNA-seq data and by that we mean to apply the same steps of the recipe to the new data

proof_genes_pt_gse2034_baked <- bake(trained_rec, new_data = proof_genes_pt_gse2034)

# This is the predict phase

gse2034_results <- predict(final_fit, new_data = proof_genes_pt_gse2034_baked, type = "linear_pred") %>%
  bind_cols(proof_genes_pt_gse2034)


# 6.3 To get the p value

validation_test <- coxph(Surv(EVENT_MON,  EVENT_STAT) ~ .pred_linear_pred, data = gse2034_results)

summary(validation_test)


library(survminer)

# Create risk groups based on the median of the predictions

gse2034_results <- gse2034_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))

# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results)

# Plot

ggsurvplot(km_fit, 
           data = gse2034_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           palette = c("#E41A1C", "#377EB8"))


gse2034_results <- gse2034_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
gse2034_results$risk_group <- relevel(gse2034_results$risk_group, ref = "Low Risk")


summary_gse2034 <- summary(coxph(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results))


library(timeROC)

# Area under the curve per time

res_auc <- timeROC(T = gse2034_results$EVENT_MON,
                   delta = gse2034_results$ EVENT_STAT,
                   marker = -gse2034_results$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 120), # 3, 5, and 10 years
                   iid = TRUE)

# View the AUC values

print(res_auc$AUC)


# Multivariate regression cox with clinical data

proof_genes_pt.gse2034.cox <- 
  proof_genes_pt_gse2034 %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("file_name") %>% 
  mutate(SCORE = gse2034_results$.pred_linear_pred
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                SCORE
  ) %>% 
  na.omit()

independent_prog.gse2034 <- coxph(surv_obj ~ SCORE, 
                                  data = proof_genes_pt.gse2034.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)
