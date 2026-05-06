

library(ggrepel)
library(paletteer)


# 1.- Lasso impact analysis -----------------------------------------------

#> First we do a plot that evaluates how genes are penalized on a Lasso model when lambda increases
#> For that first we extract coefficients from the final fit


tidy_coeffs <- tidy(extract_fit_engine(final_fit))

# 1.2 We then do an inverse log so as that to plot the estimates with relationship with the lambda

plot_data <- tidy_coeffs %>%
  mutate(log_lambda = log(lambda))

# 1.3 Filter for just the endpoints (where log_lambda is at its maximum) so that when we assign labels it only assigns to that point

endpoints <- plot_data %>%
  group_by(term) %>%
  filter(log_lambda == min(log_lambda)) %>%
  ungroup()

# 1.4 Plot

ggplot(plot_data, aes(x = log_lambda, y = estimate, group = term, color = term)) +
  geom_line(alpha = 0.7) + # Make lines slightly transparent so labels pop
  geom_text_repel(
    data = endpoints,
    aes(label = term),
    size = 3,
    hjust = 0,
    nudge_x = - 3,             # Pushes labels further to the right
    direction = "both",           # Forces labels to only move up/down to avoid overlap
    segment.color = "grey70",
    box.padding = 1,         # Reduces space around the text boxes
    max.overlaps = Inf         # Ensures every gene gets a name, even if it's crowded
  ) + 
  theme(plot.margin = margin(r = 500)) +
  scale_color_viridis_d(option = "turbo") + # High-contrast palette for many genes
  theme_minimal() +
  theme(legend.position = "none") + # Legend is redundant since names are on lines
  expand_limits(x = max(plot_data$log_lambda) - 1.5) # Make room for the names



# 2.- Identifying outlier patients ----------------------------------------



# 2.1 Add a row ID so we can find these patients later

id <- train_data %>% 
  rownames_to_column("PATIENT_ID")

# 2.2 Add to the final fit a column of id and keep

model_diagnostics <- augment(
  final_fit, 
  new_data = id, 
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)



# 3.3 Unnest the predictions and find the biggest outliers on a set point and event

outliers <- model_diagnostics %>%
  dplyr::select(PATIENT_ID, EVENT_MON, EVENT_STAT, .pred) %>%
  unnest(.pred) %>%
  filter(.eval_time == 120) %>%
  filter(EVENT_STAT == 1) %>% 
  arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive

print(outliers)

# 3.4 Observe the claudin type of the outliers

for (i in 1:20) {
  
  print(paste0(metadata.ER_POS_SURV$RFS_MONTHS[metadata.ER_POS_SURV$PATIENT_ID == as.character(outliers[i, 1])], " ", metadata.ER_POS_SURV$CLAUDIN_SUBTYPE[metadata.ER_POS_SURV$PATIENT_ID == as.character(outliers[i, 1])], " ", outliers$PATIENT_ID[i]))

  }


# 3.5 Create bias score 

outliers <- outliers %>%
  mutate(
    actual_survival = 1 - EVENT_STAT, # If they lived its 1 if they died its 0
    bias_score = abs(actual_survival - .pred_survival) # If they lived but where scored low then it will have a high bias score (1 - 0.2 = 0.8) and if they died but where scored high they will also have a high absolute bias (abs(0 - 0.9) = 0.9)
  ) %>%
  arrange(desc(bias_score))

# 3.6 IDs of patients identified as higher bias

top_bias_ids <- outliers$PATIENT_ID[1:20]


# 4.- Gene global contribution to score -----------------------------------


# 4.1 Create a data frame that cntains the estimates, if its protectigve or risl and the importance in positive values

coef_df <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  mutate(direction = ifelse(estimate < 0, "Protective", "Risk"),
         importance = estimate ^ 2) 

# 4.2 Plot descending importance

ggplot(coef_df, aes(x = estimate, y = reorder(term, importance), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Risk" = alpha("#b14037", 0.8), "Protective" = alpha("#55a5c4", 0.8))) +
  labs(x = "Coefficient (log hazard)", y = NULL) +
  theme_minimal()



# 5.- How the signature scores different categories -----------------------


# 5.1 Boxplot of the score with respect to PAM50

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material") 

# 5.2 Scatter plot of the score with relation with lymph nodes

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = LYMPH)) +
  geom_point() 
  
  
# 5.3 For plotting the score with respect to diverse treatment modalities
# 5.3.1 Change to pivot longer with respect to the treatment types and add a column with treatment type and another with which treatment did the patient recieve

proof_genes_long <- 
  proof_genes_pt.cox %>%
  pivot_longer(
    cols = c(HORMONE, CHEMO, SURGERY), # Add any other parameters here
    names_to = "Parameter",
    values_to = "Value"
  )

# 5.3.2 Plot using the new Value column for X and Parameter for facets

ggplot(proof_genes_long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() +
  facet_wrap( ~ Parameter + EVENT_STAT, scales = "free_x") + # free_x allows different labels for each box
  scale_fill_paletteer_d("ggsci::deep_purple_material") +
  theme_minimal() +
  labs(x = "Treatment/Parameter Value") +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

# 5.3.3 Kruskal test


kruskal.test(SCORE ~  PAM50, data = proof_genes_pt.cox)

# 5.4 Multivariate cox with treatment types

coxph(surv_obj ~ PAM50 * SCORE, data = proof_genes_pt.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)



# 6.- Resampling ----------------------------------------------------------


set.seed(456) 

# 6.1 Fold 20 times the whole data

final_resamples <- vfold_cv(proof_genes_pt, v = 20, strata = EVENT_STAT)

# 6.2 Apply the final workflow to the resamples

final_resample_results <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = metric_set(concordance_survival)
  )

# 6.2.2 Metrics

collect_metrics(final_resample_results)

# 6.3 Define the metrics with specific time points to resample for AUC

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 6.3 .2 Apply the workflow to the resamples with respect to the AUC

set.seed(123)

res_auc_5y <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = survival_metrics,
    eval_time = c(36, 60, 120) # 5 and 10 years
  )



collect_metrics(res_auc_5y)

# 6.4 Make a df with the resample metrics

resamples <- bind_rows(final_resample_results$.metrics)

# 6.4.2 Again create data frame but this time with the distinct areas under the curve 

resamples_auc <- 
  bind_rows(res_auc_5y$.metrics) %>% 
  na.omit()

# 6.4.2.2 Create object wiuth mean, standard deviation, 95% confidence interval and standard error of the different time points

summary_auc <- resamples_auc %>%
  group_by(.eval_time) %>%
  summarise(
    mean = mean(.estimate),
    sd = sd(.estimate),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

# 6.5 Add the fold numbers

resamples$fold <- final_resample_results$id

# 6.6 Plot The c score obtained at every fold and add a line with the median of the c scores obtained with the 20 folds

ggplot(resamples, aes(x = reorder(fold, .estimate), y = .estimate)) + # Reorder so that it gives an increasing graph
  geom_point() +
  geom_hline(yintercept = mean(resamples$.estimate), linetype = "dashed", color = "red") +
  annotate("text", x = 1.8, y = mean(resamples$.estimate),
           label = paste0("Mean = ", round(mean(resamples$.estimate), 3)),
           vjust = -1) + 
  labs(x = "Folds", y = "Estimate", title = "Refold: C-score")+
  annotate("text", x = 1.8, y = min(resamples$.estimate),
           label = paste0("Min = ", round(min(resamples$.estimate), 3)),
           vjust = -1) +
  annotate("text", x = length(resamples$.estimate) - 0.5, y = max(resamples$.estimate),
           label = paste0("Max = ", round(max(resamples$.estimate), 3)),
           vjust = 2)


# 6.6.2 Plot the mean areas under the curve with its confidence intervals

ggplot(summary_auc, aes(x = factor(.eval_time), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Time", y = "AUC", title = "Refold: Time-dependent AUC (mean ± 95% CI)") +
  annotate("text", x = 1, y = summary_auc$mean[1],
           label = paste0("Mean = ", round(summary_auc$mean[1], 3)),
           vjust = -1) +
  annotate("text", x = 2, y = summary_auc$mean[2],
           label = paste0("Mean = ", round(summary_auc$mean[2], 3)),
           vjust = -1) +
  annotate("text", x = 3, y = summary_auc$mean[3],
           label = paste0("Mean = ", round(summary_auc$mean[3], 3)),
           vjust = -1)


# 6.7 From the fold outlier identify different clinical characteristics

for (i in c(8, 9, 20, 18)) { # Here add the number of the folds to analyze
  
refold <- final_resamples$splits[[i]] # Gives the split of each fold

# Patients used for training in fold i

train_ids <- analysis(refold)

# Patients used for validation in fold i

test_ids <- assessment(refold)

# Uses the metadata and obtains only the patients in i folds test set

meta_testid <- ml_metadata[ml_metadata$PATIENT_ID %in% rownames(test_ids),]

# Print the counts of desired characteristics of patrients in the i fold test data

print(test_ids %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(meta_testid, by = "PATIENT_ID") %>% 
  dplyr::count(CLAUDIN_SUBTYPE, OS_STATUS))


# Aditionally print the names and id of patient identified as high bias of the train set

print(paste0(metadata.ER_POS_SURV$CLAUDIN_SUBTYPE[metadata.ER_POS_SURV$PATIENT_ID %in% top_bias_ids[top_bias_ids %in% rownames(test_ids)]],
       " ",
       top_bias_ids[top_bias_ids %in% rownames(test_ids)])
)

}


# 7.- Observe distributions and shapiro -----------------------------------


# 7.1 Shapiro wilk test of desired distribution

shapiro.test(train_data2$risk_score)

# 7.2 Plot distribution

ggplot(train_data2, aes(x = risk_score)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "#81dcff", size = 1) +
  theme_minimal()


# 7.3 Plot distribution based on event stat

ggplot(test_data, aes(x = risk_score, fill = factor(EVENT_STAT, labels = c("Alive", "Diceased")))) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(fill = "Event",
       x = "Risk score",
       y = "Density",
       title = "Distribution of score with respect to event")

# Convert to factor

proof_genes_pt.cox <- 
  proof_genes_pt.cox %>%
  mutate(EVENT_STAT = factor(EVENT_STAT)) 

# 7.4 Wilcox test of desired data

results_wilcox <- proof_genes_pt.cox %>%
  dplyr::select(EVENT_STAT, PAM50, SCORE) %>% 
  group_by(PAM50) %>%
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "BH"))

print(results_wilcox)

# 7.4.2 Plot the  comparisons

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material") + 
  facet_wrap(~ EVENT_STAT, 
             labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased")))


# 7.4.3 Prepare to plot the different comparisons

proof_genes_pt.long<- 
  proof_genes_pt.cox %>% 
  pivot_longer(
    cols = c(SURGERY, CHEMO), 
    names_to = "Parameter",
    values_to = "Value"
  )

# 7.4.3.2 Plot

ggplot(data = proof_genes_pt.long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() + 
  facet_wrap(~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2) + 
  scale_fill_paletteer_d("palettesForR::Pastels") + 
  theme_gray(base_size = 18)
