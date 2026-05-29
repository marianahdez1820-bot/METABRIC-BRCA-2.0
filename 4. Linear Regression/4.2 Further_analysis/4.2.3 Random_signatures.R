library(tidyverse)
library(survival)
library(foreach)
library(doParallel)


# 1.- Cox on our gene set -------------------------------------------------

# 1.1 Counts set with only patients with metadata transposed and matrix to scale

x_matrix_full <- counts_data[, ml_metadata$PATIENT_ID] %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% rownames(proof_genes_pt)) %>% 
  as.matrix()

# 1.2.1 Scale and center

x_matrix_full <- scale(x_matrix_full)

# 1.2.2 Names of all the genes for later resampling into random signatures

all_genes_pool <- colnames(x_matrix_full)

# 1.3 Add objects needed for survival analysis

full_df_unpenalized <- as.data.frame(x_matrix_full)
full_df_unpenalized$EVENT_MON <- proof_genes_pt$EVENT_MON
full_df_unpenalized$EVENT_STAT <- proof_genes_pt$EVENT_STAT


# 1.4 Construct direct multivariate unpenalized Cox formula

true_formula <- as.formula(paste("Surv(EVENT_MON, EVENT_STAT) ~", paste(proof_genes, collapse = " + ")))

true_cox_fit <- coxph(true_formula, data = full_df_unpenalized)

# 1.4.2 C score of cox model

true_unpenalized_c_index <- summary(true_cox_fit)$concordance[1]

# 2.- Paralelization ------------------------------------------------------

# 2.1 Prepare cores

num_cores <- parallel::detectCores() - 1

cl <- makeCluster(num_cores)

registerDoParallel(cl)

# 2.2 This establishes how many random genes will be samples

sig_size <- length(proof_genes)

# 2.3 Number of random sets

n_permutations <- 1000

# 2.4 Export required clean variables to the processing cores

clusterExport(cl, c("full_df_unpenalized", "all_genes_pool", "sig_size"))

# 2.5 Paralelize 1000 permutations

unpenalized_random_c_indices <- foreach(
  i = 1:n_permutations, 
  .combine = 'c', 
  .packages = c('survival')
) %dopar% {
  
  set.seed(123 + i) # Set different seeds
  
  # Pick number of  random genes that have the same length as the signature evaluated
  
  random_genes <- sample(all_genes_pool, sig_size, replace = FALSE)
  
  # Formulate standard unpenalized multivariate regression
  
  random_formula <- as.formula(paste("Surv(EVENT_MON, EVENT_STAT) ~", paste(random_genes, collapse = " + ")))
  
  # Fit model safely to catch non-convergence instances
  
  c_index <- tryCatch({
    fit_cox_random <- coxph(random_formula, data = full_df_unpenalized)
    summary(fit_cox_random)$concordance[1]
  }, error = function(e) return(NA))
  
  return(c_index)
}

stopCluster(cl)
unpenalized_random_c_indices <- na.omit(unpenalized_random_c_indices)


# 5. Significance estimation

# 5.1 Obtain how many random signatures got a better or equal C index than our gene set

better_or_equal_random <- sum(unpenalized_random_c_indices >= true_unpenalized_c_index)

# 5.1.2 Empirical p dividing the previous number over the number of permutations

empirical_p <- better_or_equal_random / length(unpenalized_random_c_indices)

# CHeck normality

shapiro.test(unpenalized_random_c_indices)

# 5.2 Compute the Z-score and analytical p

z_score <-(true_unpenalized_c_index - mean(unpenalized_random_c_indices)) / sd(unpenalized_random_c_indices)

analytical_p <- 2 * pnorm(-abs(z_score))

cat("Unpenalized True Signature C-index   :", round(true_unpenalized_c_index, 4), "\n")
cat("Median Performance of Random Sets    :", round(median(unpenalized_random_c_indices), 4), "\n")
cat("Empirical P-value vs Random Noise    :", empirical_p, "\n")
