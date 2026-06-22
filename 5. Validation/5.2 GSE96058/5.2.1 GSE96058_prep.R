#> In this script we obtain and preprocess the metadata.gse96058 and the counts data

library(GEOquery)
library(tidyverse)

# 1.- Download data -----------------------------------------------------------


# 1.1 Download the supplementary file (The actual expression matrix)
# getGEOSuppFiles("GSE96058", baseDir = "D:/GSE96058")

# 1.1.2 Read the specific expression file (SCAN-B typically provides a large .txt or .csv)

raw <- read.csv("D:/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", 
                row.names = 1, check.names = FALSE)
counts_data.gse96058 <- raw
# 1.2 Download metadata.gse96058

gse <- getGEO("GSE96058", GSEMatrix = TRUE)

# 1.2.2 Asign to object

pheno <- pData(gse[[1]])



# 2.- Preprocess metadata.gse96058 -------------------------------------------------

pheno$characteristics_ch1.3 <- gsub("\\D", "", pheno$characteristics_ch1.3)


metadata.gse96058 <-
  pheno %>%
  mutate(
    tissue = source_name_ch1,
    age = as.numeric(`age at diagnosis:ch1`),
    tumor_size = as.numeric(`tumor size:ch1`),
    lymph_group = `lymph node group:ch1`,
    lymph_status =  `lymph node status:ch1`,
    er_status = as.numeric(`er status:ch1`),
    pgr_status = as.numeric(`pgr status:ch1`),
    her2_status = as.numeric(`her2 status:ch1`),
    ki67_status = as.numeric(`ki67 status:ch1`),
    nhg = as.factor(`nhg:ch1`),
    er_pred_mgc = as.numeric(`er prediction mgc:ch1`),
    # This are predictions made by RNA if mgc its Molecular Gene Classifier which is older than SCN which is Single Sample Classifier (SCAN-B) and if SGC its single gene classifier
    er_pred_sgc = as.numeric(`er prediction sgc:ch1`),
    pgr_pred_mfc = as.numeric(`pgr prediction mgc:ch1`),
    pgr_pred_sgc = as.numeric(`pgr prediction sgc:ch1`),
    her2_pred_mfc = as.numeric(`her2 prediction mgc:ch1`),
    her2_pred_sgc = as.numeric(`her2 prediction sgc:ch1`),
    ki67_pred_mfc = as.numeric(`ki67 prediction mgc:ch1`),
    ki67_pred_sgc = as.numeric(`ki67 prediction sgc:ch1`),
    nhg_pred_mgc = as.numeric(`nhg prediction mgc:ch1`),
    pam50 = as.factor(`pam50 subtype:ch1`),
    os_months = as.numeric(`overall survival days:ch1`) / 30.4166667,
    os_status = as.numeric(`overall survival event:ch1`),
    endocrine_tx = as.numeric(`endocrine treated:ch1`),
    chemo_tx = as.numeric(`chemo treated:ch1`)
    
    
  ) %>%
  dplyr::select(
    -c(
      source_name_ch1,
      characteristics_ch1.2,
      characteristics_ch1.3,
      characteristics_ch1.4,
      characteristics_ch1.5,
      characteristics_ch1.6,
      characteristics_ch1.7,
      characteristics_ch1.8,
      characteristics_ch1.9,
      characteristics_ch1.10,
      characteristics_ch1.11,
      characteristics_ch1.12,
      characteristics_ch1.13,
      characteristics_ch1.14,
      characteristics_ch1.15,
      characteristics_ch1.16,
      characteristics_ch1.17,
      characteristics_ch1.18,
      characteristics_ch1.19,
      characteristics_ch1.20,
      characteristics_ch1.21,
      characteristics_ch1.22,
      characteristics_ch1.23,
      characteristics_ch1.24,
      # Each one of the characteristics_ch1. corresponds to its equivalent in the next lines and both correspond in orther to its characteristic in mutate
      `age at diagnosis:ch1`,
      `tumor size:ch1`,
      `lymph node group:ch1`,
      `lymph node status:ch1`,
      `er status:ch1`,
      `pgr status:ch1`,
      `her2 status:ch1`,
      `ki67 status:ch1`,
      `nhg:ch1`,
      `er prediction mgc:ch1`,
      `er prediction sgc:ch1`,
      `pgr prediction mgc:ch1`,
      `pgr prediction sgc:ch1`,
      `her2 prediction mgc:ch1`,
      `her2 prediction sgc:ch1`,
      `ki67 prediction mgc:ch1`,
      `ki67 prediction sgc:ch1`,
      `nhg prediction mgc:ch1`,
      `pam50 subtype:ch1`,
      `overall survival days:ch1`,
      `overall survival event:ch1`,
      `endocrine treated:ch1`,
      `chemo treated:ch1`
      
    )
  )

counts_data.gse96058[1:5,1:5]


metadata.gse96058_er_pos <-
  metadata.gse96058 %>% 
  filter(er_status == 1) %>% 
  rownames_to_column("id") %>% 
  mutate(EVENT_STAT = os_status,
         EVENT_MON = os_months,
         id = NULL) 

# 3.- Preprocess data -----------------------------------------------------

# 3.1 Match it with metadata.gse96058

# 3.1.1 Identify patients in both sets

common_samples <- intersect(colnames(counts_data.gse96058), metadata.gse96058_er_pos$title)

# 3.1.2 Keep the patients in counts data that also have metadata.gse96058

counts_data.gse96058_erpos <- counts_data.gse96058[, common_samples]

counts_data.gse96058_erpos <- t(counts_data.gse96058_erpos)



# 3.2 Find genes present in both data sets

colnames(counts_data.gse96058_erpos) <- make.names(colnames(counts_data.gse96058_erpos))

common_genes_meta.gse96058 <- intersect(proof_genes, colnames(counts_data.gse96058_erpos))


# 4.2 Object with all the patients and expression of only the genes of interest

counts_data.gse96058_erpos <- counts_data.gse96058_erpos[, colnames(counts_data.gse96058_erpos) %in% common_genes_meta.gse96058]

# 4.3 Stop running if there are less genes in GSE96058 than on the signature

if(length(common_genes_meta.gse96058) < length(proof_genes)){
  stop(paste("There are missing genes in GSE96058 relative to the signature, missing ", length(proof_genes) - length(common_genes_meta.gse96058), " gene(s): "),  paste0(proof_genes[!(proof_genes %in% common_genes_meta.gse96058)], sep = ", ")) # Script stops here
}else{
  print("All genes in the signature are on GSE96058")
}


# 4.4 Keep only genes that are to be evaluated

counts_data.gse96058_erpos <- counts_data.gse96058_erpos[ , common_genes_meta.gse96058]

# 4.5 Object with surv object, all of the genes and the event columns

proof_genes_pt.gse96058 <- 
  counts_data.gse96058_erpos %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title") %>% 
  mutate(surv_obj = Surv(time = EVENT_MON, event = EVENT_STAT, type = "right")) %>% 
  dplyr::select(all_of(proof_genes),
                EVENT_STAT,
                EVENT_MON,
                title,
                surv_obj) %>% 
  column_to_rownames("title")


proof_genes_pt.gse96058 <- 
  proof_genes_pt.gse96058 %>% 
  mutate(across(-c(EVENT_MON, EVENT_STAT, surv_obj), ~ as.vector(scale(.x))))
