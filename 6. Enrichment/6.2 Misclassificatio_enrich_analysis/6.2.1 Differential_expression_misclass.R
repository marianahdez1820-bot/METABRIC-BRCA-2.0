library(limma)
library(tidyverse)
# 1.- Dividing by BIAS nodes ------------------------------



# 1.1 Generate column corresponding to BIAS nodes, those that have 0 in one group and those with more than 0 in another

col_data <- ml_metadata %>%
  mutate(BIAS = factor(ifelse(PATIENT_ID %in% misclassification_diff_id , 
                1,
                0))) %>% 
  filter(!(BIAS == 1 & EVENT_STAT == 0)) %>%
  # mutate(BIAS = factor(ifelse(EVENT_STAT == 1,
  #                      1,
  #                      0))) %>%
  dplyr::select(PATIENT_ID, BIAS) %>% # Create only the object to use for Limma
  column_to_rownames("PATIENT_ID")




# 2.- Differential expression -----------------------------------------------


# 2.2 Data counts of the patients that had BIAS node information in the metadata

count_data <- counts_data[colnames(counts_data) %in% rownames(col_data)]

# 2.2.2 Making shure they are in the same order

count_data <- count_data[match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))


# 2.3 Generate limma object

library(limma)

# 2.4 Design based on object separating on BIAS nodes

design <- model.matrix(~ 0 + BIAS, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design)) 

# 2.5 Fit

fit <- lmFit(count_data, design)

# 2.5.2 Contrast matrix comparing BIAS 0 to > 0 BIAS

contrast.matrix <- makeContrasts(BIAS1 - BIAS0,
                                 levels = design)

# 2.5.3 Fit based on contrasts

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit)

# 2.6 Results

res <- topTable(fit, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) # 0.1
res_sig

intersect(proof_genes, rownames(res_sig))

