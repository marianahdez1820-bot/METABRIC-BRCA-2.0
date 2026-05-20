  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(survival)
  library(paletteer)
  
  
  
  # 1.- Loading data --------------------------------------------------------
  
  
  
  # 1.1 Query con base a RNA-Seq y a 3 casos y 3 controles
  
  # tcga_rna <- GDCquery("TCGA-BRCA",
  #                      data.category = "Transcriptome Profiling",
  #                      access = "open",
  #                      experimental.strategy = "RNA-Seq",
  #                      workflow.type = "STAR - Counts"
  # )
  
  # GDCdownload(tcga_rna, method = "api", files.per.chunk = 5)
  
  # 1.2 Prepare data for usage
  
  #tcga_brca_data <- GDCprepare(tcga_rna, directory = "D:/tcga/GDCdata")
  
  # 1.3 Count matrix
  
  # brca_matrix <- assay(tcga_brca_data, "fpkm_unstrand")
  
  # brca_matrix %>%
  #   write.csv("fpkm_unstrand.csv")
  
  brca_matrix <- read.csv("fpkm_unstrand.csv") %>% 
    column_to_rownames("X")
  
  
  # 1.4 Convert to data frame
  
  brca_data <- brca_matrix %>% 
    as.data.frame()
  
  # 1.4.2 Convert the dots to slashes
  
  colnames(brca_data) <- gsub("\\.", "-", colnames(brca_data))
  
  
  # 1.5.- Eliminate duplicates ----------------------------------------------
  
  
  # 1.5 Extract sample type from TCGA barcode
  
  # 1.5.2 Select the 14th - 16th value which correspond to sample type codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  # this to then select only the samples that correspond to primary tumors
  
  sample_type_full <- substr(colnames(brca_data), 14, 16)
  
  # 1.5.3 Maintain only primary tumor samples so as to avoid duplicates
  
  # 1.5.3.1 Keep only counts that correspond to primary tumor
  
  brca_data <- brca_data[, sample_type_full == "01A"]
  
  # 1.5.4 Assign to new object that will be modified to eliminate duplicates
  
  brca_data2 <- brca_data
  
  # 1.5.4.2 Keep the names up until the -01 so as to have it in the same nomenclature as the metadata
  
  colnames(brca_data2) <- substr(colnames(brca_data2), 1, 15)
  
  # 1.5.4.3 We can see that there are 5 patients with 2 samples of the same tumor
  
  names(brca_data)[substr(colnames(brca_data), 1, 15) %in% names(brca_data2)[duplicated(names(brca_data2))]]
  
  # 1.5.5 So we keep only the patient sample with highest variance
  
  # 1.5.5.2 Calculate variance
  
  sample_variance <- apply(brca_data2, 2, var, na.rm = TRUE)
  
  # 1.5.5.3 Match variance with the names 
  
  df <- data.frame(
    sample = colnames(brca_data2), # Column with the colnames
    equivalent = colnames(brca_data), # Column with the full names so as to then reassign to the counts data
    variance = sample_variance # column with the variance
  )
  
  # 1.5.5.3 Group by and keep only the sample with highest variance
  
  selected_samples <- df %>%
    group_by(sample) %>% # Group by the reduced name (TCGA-XX-XXXX-01)
    slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # Keep only the highest variance
    pull(equivalent) # Extract the complete name (TCGA-XX-XXXX-01A-XXX-XXXX-XX)
  
  # 1.6 Keep only the counts of the unique
  
  brca_data_unique <- brca_data[, selected_samples]
  
  # 1.7 Convert back to names compatible with metadata
  
  colnames(brca_data_unique) <- substr(colnames(brca_data_unique), 1, 15)
  
  # 2.- Metadata ------------------------------------------------------------
  
  library(UCSCXenaTools)
  
  
  # 2.1 Generate and Query
  
  data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
    XenaQuery()
  
  # 2.2 Download 
  
  xe_download <- XenaDownload(data_query, destdir = "D:/tcga/GDCdata/Metadata")
  
  # 2.3 Prepare (Load) the data
  
  brca_clinical <- XenaPrepare(xe_download)
  
  # 2.4 Create the Recurrence variables
  
  
  refined_data <- 
    brca_clinical %>%
    mutate(
      SURVIVAL = ifelse(vital_status == "DECEASED", 1, 0),
      SURVIVAL_MON = ifelse(!is.na(days_to_death) & days_to_death >= 0,
                            days_to_death, # If they died they have this parameter
                            days_to_last_followup # Else they have this parameter
      ) / 30.4166667,
      RECURRENCE = ifelse(
        new_neoplasm_event_type %in% c("Locoregional Recurrence", "Distant Metastasis"),
        1,
        0
      ),
      RECURRENCE_MON = ifelse(
        RECURRENCE == 1,
        days_to_new_tumor_event_after_initial_treatment,
        coalesce(days_to_last_followup, days_to_death)
      ) / 30.4166667
    ) %>%
    filter(
      breast_carcinoma_estrogen_receptor_status == "Positive" # Only ER+,
    ) 
  
  
  
  # 2.5 To maintain the samples of the primary tumor and eliminate duplicates as done with counts data
  
  sample_type_full2 <- substr(refined_data$sampleID, 14, 15)
  
  refined_data_unique <- refined_data
  
  refined_data_unique <- refined_data_unique[sample_type_full2 == "01",]
  
  
  # 2.6 Count object that corresponds to the metadata patients
  
  brca_data_filtered <- brca_data_unique[, colnames(brca_data_unique) %in% refined_data_unique$sampleID]
  
  # 2.7 Since there are missing patients since the first data we also filter out patients in metadata that arent on counts 
  
  refined_data_unique <- refined_data_unique[refined_data_unique$sampleID %in% colnames(brca_data_filtered),]
  
  
  
  # 3.- Deleting duplicates and asigning ensembl as rownames ----------------
  
  # 3.1 Deleting the version of ensembl and keeping only the full name
  
  counts_data_duplicates <- 
    brca_data_filtered %>% 
    rownames_to_column("ensembl_version") %>% 
    mutate(ensembl = gsub("\\..*", "", ensembl_version)) %>% 
    column_to_rownames("ensembl_version")
  
  # 3.2 Asign symbol as a column
  
  counts_data_duplicates$symbol <- mapIds(
    org.Hs.eg.db,
    keys = counts_data_duplicates$ensembl, # A donde va a buscar
    column = "SYMBOL",  # Nueva columna con ese formato
    keytype = "ENSEMBL",  # Que formato va a buscar en keys
    multiVals = "first") # Que hacer si hay varios del mismo Key
  
  
  # 3.2 Variance
  
  numeric_data <- counts_data_duplicates %>%
    dplyr::select(where(is.numeric))
  
  counts_data_duplicates$variance <- apply(numeric_data, 1, var)
  
  
  # 3.2.2 Only mantain the version of the gene duplicate with higher variance
  
  counts_data.tcga <- counts_data_duplicates %>% # Initial data
    group_by(symbol) %>% # Group by ensembl
    slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # Order by variance and keep the highest
    ungroup() %>% 
    rownames_to_column("genes") %>% 
    dplyr::select( - variance, 
                   - ensembl,
                   - genes) %>% # Delete variance, and both of the ensembl ids columns
    filter(!is.na(symbol)) %>%  # Delete those that had a NA in symbol
    column_to_rownames("symbol")
  
  
  # 4.- Asigning signature genes --------------------------------------------
  
  # 4.1 List of genes described in the differential expression as being of prognosis for late death
  
  rownames(counts_data.tcga) <- make.names(rownames(counts_data.tcga))
  
  common_genes_meta.tcga <- intersect(proof_genes, rownames(counts_data.tcga))
  
  
  # 4.2 Object with all the patients and expression of only the genes of interest
  
  proof_genes_pt.tcga <- counts_data.tcga[rownames(counts_data.tcga) %in% common_genes_meta.tcga, ]
  
  # 4.3 Stop running if there are less genes in TCGA than on the signature
  
  if(length(common_genes_meta.tcga) < length(proof_genes)){
      stop(paste("There are missing genes in TCGA relative to the signature, missing ", length(proof_genes) - length(common_genes_meta.tcga), " gene(s): "),  paste0(proof_genes[!(proof_genes %in% common_genes_meta.tcga)], sep = ", ")) # Script stops here
  }else{
    print("All genes in the signature are on TCGA")
  }
  
  
  
  # 5.- Prepare object for validation ---------------------------------------
  
  # 5.1 Transpose first so samples are rows, dont scale because the bake in recipe does that 
  
  tcga_transposed <- 
    t(proof_genes_pt.tcga)
  
  # 5.2 Log2 Transform 
  
  tcga_log <- 
    log2(tcga_transposed + 1)
  
  # 5.3 Assign log to the object for ML
  
  proof_genes_pt.tcga <- 
    tcga_log %>% as.data.frame()
  
  
  # 5.4 Check that the patients are in the same order
  
  refined_data_unique <- 
    refined_data_unique[refined_data_unique$sampleID %in% rownames(proof_genes_pt.tcga),]
  
  
  all(rownames(proof_genes_pt.tcga) == refined_data_unique$sampleID)
  
  # 5.5 Add a column of EVENT as a binary term for it to be the outcome
  
  proof_genes_pt.tcga <- 
    proof_genes_pt.tcga %>% 
    rownames_to_column("sampleID") %>% 
    left_join(refined_data_unique, by = "sampleID") %>%  # Join counts with metadata
    column_to_rownames("sampleID") %>% 
    dplyr::select(all_of(proof_genes),  # Keep all the genes to ve evaluated and the oucome variables
                  RECURRENCE_MON, # SURVIVAL_MON for survival and RECURRENCE_MON for recurrence
                  RECURRENCE) %>% # SURVIVAL for survival and RECURRENCE for recurrence
    dplyr::rename(EVENT_STAT = RECURRENCE, # Rename to common term (EVENT_STAT for event and EVENT_MON for time of follow up)
           EVENT_MON = RECURRENCE_MON) %>% 
    mutate(EVENT_STAT = as.numeric(EVENT_STAT),
           EVENT_MON = as.numeric(EVENT_MON)
           ) %>%  
    as.data.frame() %>% 
    mutate(surv_obj =  Surv( # Create survival object
      time  = EVENT_MON,
      event = EVENT_STAT
    )) %>% 
    filter(
      EVENT_MON > 0
    )
  
  
  outcome_analyzed <- "RECURRENCE" # To print at the end so as to ot get confused to what is being analyzed
  
  if(outcome_analyzed == "RECURRENCE"){
    proof_genes_pt.tcga <- 
      proof_genes_pt.tcga %>% 
      filter(EVENT_MON >= 2 & EVENT_MON <= 180)
    print("Recurrene signature, filtered")
  }else{
    "Survival signature, no filter"
}

  
proof_genes_pt.tcga <- 
    proof_genes_pt.tcga %>% 
    mutate(across(- c(EVENT_MON, EVENT_STAT, surv_obj), scale))
  
