library(tidyverse)

# 1.- Data loading --------------------------------------------------

microarray_data <- read.table("D:/brca_metabric/brca_metabric/data_mrna_illumina_microarray.txt",
                              header = TRUE)

# 1.1.2 Metadata loading

clinical_data <- read.delim(
  "D:/brca_metabric/brca_metabric/data_clinical_patient.txt",
  comment.char = "#",
  stringsAsFactors = FALSE
)

# 1.1.3 Check if rownames in clinical data have the same nomenclature as the colname in microarray_data, if they differ by a point or dash (MB.0000 vs MB-0000) use
clinical_data$PATIENT_ID <- gsub("-", ".", clinical_data$PATIENT_ID)

# 1.1.4 Object without ENTREZ and HUGO columns

microarray_data.num <- microarray_data[,-(1:2)] 

# 1.2.- Eliminate duplicates and Hugo IDs -----------------------

# 1.2.1 There are Hugo IDs duplicated, and thus we search for variance and keep the duplicate with higher variance

# 1.2.2 Variance

microarray_data.var <- apply(microarray_data.num, 1, var)

# 1.2.3 Asign new column with variance

microarray_data$variance <- microarray_data.var

# 1.2.4 Only maintain the version of the gene duplicate with higher variance

microarray_data.unique <- microarray_data %>% # Initial data
  group_by(Hugo_Symbol) %>% # Group by HUGO symbol
  slice_max(order_by = variance, n = 1) %>% # Order by variance and keep the highest
  ungroup()


# 1.2.5 Convert from tibble to data frame

microarray_data.unique <- as.data.frame(microarray_data.unique)

# 1.2.6 Asign Hugo symbol as rownames

rownames(microarray_data.unique) <- microarray_data.unique$Hugo_Symbol

# 1.2.7 Delete the columns with hugo, entrez and variance.

microarray_data.unique$Hugo_Symbol <- NULL
microarray_data.unique$Entrez_Gene_Id <- NULL
microarray_data$variance <- NULL

counts_data <- microarray_data.unique
counts_data$variance <- NULL

##############################################

# 2.- METADATA -----------------------


# 2.1 Metadata that contains only the patients of whom we have RNA counts

metadata <- clinical_data[clinical_data$PATIENT_ID %in% colnames(microarray_data.num),]

metadata  <-   metadata %>% 
  mutate(ER_POS = ifelse(ER_IHC == "Positve",
                1,
                0),
         SURVIVAL_STAT = as.numeric(gsub(":.*", "", metadata$OS_STATUS)), # 2.2 Add column of survival state where 1 is dead and 0 is alive
         RECURR_STAT = as.numeric(gsub(":.*", "", metadata$RFS_STATUS))
         )
metadata



# 2.3 Metadata of only the patients that are dead because of breast cancer

metadata_diseased.brca <- metadata[metadata$SURVIVAL_STAT > 0 & metadata$VITAL_STATUS == "Died of Disease",]



# 1.1 Object with live patients and patients diseased by breast cancer

alive_brca.death <- metadata %>% 
  filter(VITAL_STATUS != "Died of Other Causes")
