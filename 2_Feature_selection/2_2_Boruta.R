library(Boruta)
library(survival)
library(parallel)
library(ranger)


# 2.- Metadata Preparation (Filtering for ER+) ------------------

boruta_metadata <- metadata_surv_train  %>% 
  as.data.frame() %>% 
  drop_na(EVENT_STAT, EVENT_MON) %>% # Eliminate rows with NAs
  dplyr::select(PATIENT_ID, EVENT_STAT, EVENT_MON) # Select columns to create outcome and to identify patients

# 3. Align Expression Data & Feature Engineering ------------------

# 3.1 Select count data with only ER+ patients and transpose

boruta_df <- counts_data[, boruta_metadata$PATIENT_ID] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID") %>%
  left_join(boruta_metadata, by = "PATIENT_ID") %>% # Joint with metadata to have the counts, and the objects for survival analysis
  column_to_rownames("PATIENT_ID")


# 3.2 Create the Survival Object

boruta_df$surv_obj <- Surv(time = boruta_df$EVENT_MON, event = boruta_df$EVENT_STAT)

# 3.3 Remove the survival columns that do not correspond to the surv_obj

boruta_df <- boruta_df %>% 
  dplyr::select( - EVENT_STAT, - EVENT_MON) 

# 5. Variance Filtering (Top 3000 genes) ------------------------

# 5.1 Calculate variance

vars <- apply(boruta_df[, setdiff(colnames(boruta_df), "surv_obj")], 2, var)

# 5.2 Sort based on variance and select top 3000 genes

top_genes <- names(sort(vars, decreasing = TRUE))[1:7000]

# 5.2.2 Maintain only genes that correspond to the top 3000 and the surv object

boruta_df_small <- boruta_df[, c(top_genes, "surv_obj")]

# 5.3 Fix any non-standard gene names 

colnames(boruta_df_small) <- make.names(colnames(boruta_df_small))

# 6. Parallel Function  ---------------------------------

impRangerSurv <- function(x, y, ...) {
  temp_df <- as.data.frame(x) # 1. Combine X and Y into one dataframe for ranger, x corresponds to the gene names and y to the surv obj
  temp_df$target_surv <- y
  
  # 2. Use 'dependent.variable.name' instead of 'formula'
  # This is the "fast" interface for high-dimensional data
  
  res <- ranger::ranger(
    dependent.variable.name = "target_surv", 
    data = temp_df, # Data frame created before
    importance = "permutation", 
    num.trees = 500,            
    num.threads = 16,            # Ensure threads are passed here
    ...
  )
  return(res$variable.importance)
}

# 7. Run Parallelized Boruta ---------------------------

set.seed(111)

# Prepare clean X and Y

x_data <- boruta_df_small[, setdiff(colnames(boruta_df_small), "surv_obj")]

y_data <- boruta_df_small$surv_obj

# Run boruta

boruta.signature <- Boruta(
  x = x_data,
  y = y_data,
  getImp = impRangerSurv, 
  doTrace = 3,
  maxRuns = 500
)

#  8. View Results ------------------

print(boruta.signature)
plot(boruta.signature, las = 2, cex.axis = 0.7)

boruta.signature2 <- boruta.signature

# Only the 11 confirmed genes
confirmed_only <- getSelectedAttributes(boruta.signature, withTentative = FALSE)

# The 22 tentative  genes
# (We find these by taking the full list and removing the confirmed ones)

all_selected <- getSelectedAttributes(boruta.signature, withTentative = TRUE)
tentative_only <- setdiff(all_selected, confirmed_only)

cat(confirmed_only, sep = ", ")
cat(tentative_only, sep = ", ")


# Define path
out_path <- "C:/R/METABRIC/Results/"

# 1. Force the decision on those 22 tentative genes
final_boruta_decided <- TentativeRoughFix(boruta.signature)

final_boruta_decided$finalDecision[final_boruta_decided$finalDecision == "Confirmed"]

saveRDS(boruta.signature, paste0(out_path, "boruta_surv_decided.rds"))

# 2. Save the final model object
saveRDS(final_boruta_decided, paste0(out_path, "final_boruta_surv_decided.rds"))

stats <- attStats(boruta.signature)


# 3. Get the names of all selected genes (Confirmed + Fixed Tentatives)
final_gene_names <- getSelectedAttributes(final_boruta_decided)

# 4. Save the gene list as a CSV
write.csv(final_gene_names, paste0(out_path, "final_gene_surv_signature.csv"), row.names = FALSE)

# 5. Export the importance values (the numbers used in your plot)
stats <- attStats(final_boruta_decided)
write.csv(stats, paste0(out_path, "gene_surv_importance_full_stats.csv"))


# 9.- Load data ---------------------------

# Load the final decided Boruta object
final_boruta <- readRDS(paste0(out_path, "final_boruta_decided.rds"))

# Verify it loaded correctly
print(final_boruta)

final_boruta_decided$finalDecision

cat(final_gene_names, sep = ", ")
