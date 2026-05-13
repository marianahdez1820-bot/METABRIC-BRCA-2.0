library(ggrepel)
library(paletteer)
library(patchwork)

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

# 2.2 Utilize the fitted object to make predictions based on time

model_diagnostics <- augment(
  final_fit, 
  new_data = id, 
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)

# 2.2.2 Create list to then add the results

list <- list()


#> 2.3 For loop that at each desired time point calculates the bias scores, extracts patients with high bias scores
#> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time

for (i in c(36, 60, 120)) {
  
  eval_time <- i
  
  # 2.3 Unnest the predictions and find the biggest outliers on a set point and event
  
  outliers <- model_diagnostics %>%
    dplyr::select(PATIENT_ID, EVENT_MON, EVENT_STAT, .pred) %>%
    unnest(.pred) %>%
    filter(.eval_time == eval_time) %>% 
    arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive
  
  
  # 2.4 Object with metadata and score characteristics
  
  outlier_summary <- outliers %>%
    inner_join(metadata.ER_POS_SURV, by = "PATIENT_ID", suffix = c("", ".drop")) %>%
    dplyr::select(
      PATIENT_ID, 
      .pred_survival, 
      EVENT_STAT, 
      EVENT_MON, 
      CLAUDIN_SUBTYPE, 
      LYMPH_NODES_EXAMINED_POSITIVE, 
      THREEGENE,
      INTCLUST,
      .eval_time,
      NPI,
      CELLULARITY
    ) 
  
  
  # 2.5 Create bias score for defined time
  
  outliers_bias <- outlier_summary %>%
    mutate(
      bias_score = (((1 - EVENT_STAT) - .pred_survival) ^ 2) * ((EVENT_MON - .eval_time) * (( 1 - (1 + EVENT_STAT)) + (1 - EVENT_STAT)))
    ) %>%
    arrange(desc(bias_score))
  
  # 2.6 Identify the highest bias patients
  
  # 2.6.1 Obtain mean and sd and then filter baed on patients higher than determined SD
  
  extreme_outliers <- 
    outliers_bias %>% 
    mutate(mean_bias = mean(bias_score),
           sd_bias = sd(bias_score)) %>% 
    filter(bias_score > (mean_bias + 1 * sd_bias))
  
  # 2.6.2 Obtain their IDs
  
  top_bias_ids <- extreme_outliers$PATIENT_ID
  
  print(length(top_bias_ids))
  
  # 2.6.3 Identify different characteristics of patients identified as top bias
  
  print(outliers_bias %>%
          filter(PATIENT_ID %in% top_bias_ids) %>% 
          group_by(INTCLUST, EVENT_STAT) %>%
          summarise(
            count = n(),
            avg_pred_survival = mean(.pred_survival),
            avg_event_time = mean(EVENT_MON)
          ) %>%
          arrange(desc(count))
  )
  
  # 2.7.1 Add a column identifying patients as top bias or not
  
  outliers <- 
    outliers %>% 
    mutate(quadrant = factor(ifelse(PATIENT_ID %in% top_bias_ids,
                                    1,
                                    0)))
  
  print(outliers %>% 
          group_by(EVENT_STAT) %>% 
          dplyr::count(quadrant))
  
  # 2.7.2 Plot
  
  # 2.7.2.1 Plot colored by EVENT_STAT
  
  p1 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(EVENT_STAT), shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
    stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    labs(
      title = "Event Status Distribution",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival",
      color = "Event Stat",
      shape = "Event Stat"
    ) +
    theme_minimal()
  
  # 2.7.2.2 Plot colored by quadrant
  
  p2 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = quadrant, shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) + 
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    labs(
      title = "Quadrant Analysis",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival"
    ) +
    theme_minimal()
  
  # 2.7.3 Combine and stack
  
  print( p1 / p2)
  
  list[[i]] <- top_bias_ids
  
}

# 2.8 Obtain patients found on all of the iterations of the for loop as top bias patients

bias_interesct <- intersect(intersect(list[[36]], list[[60]]) , list[[120]])

# 2.8.2 Similar but all unique so even if they appear once we register them

bias_diff_id <- unique(c(list[[36]], list[[60]], list[[120]]))

# 2.9 Observe characteristics of the bias_intersect patients

outliers_bias %>%
  filter(PATIENT_ID %in% bias_interesct) %>% 
  group_by(INTCLUST, CELLULARITY) %>%
  summarise(
    count = n(),
    avg_pred_survival = mean(.pred_survival),
    avg_event_time = mean(EVENT_MON)
  ) %>%
  arrange(desc(count))

# 3.- Gene global contribution to score -----------------------------------

# 3.1 Create a data frame that cntains the estimates, if its protectigve or risl and the importance in positive values

coef_df <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  mutate(direction = ifelse(estimate < 0, "Protective", "Risk"),
         importance = estimate ^ 2) 

# 3.2 Plot descending importance

ggplot(coef_df, aes(x = estimate, y = reorder(term, importance), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Risk" = alpha("#b14037", 0.8), "Protective" = alpha("#55a5c4", 0.8))) +
  labs(x = "Coefficient (log hazard)", y = NULL) +
  theme_minimal()



# 4.- How the signature scores different categories -----------------------


# 4.1 Boxplot of the score with respect to PAM50

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material") 

# 4.2 Scatter plot of the score with relation with lymph nodes

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = LYMPH)) +
  geom_point() 


# 4.3 For plotting the score with respect to diverse treatment modalities
# 4.3.1 Change to pivot longer with respect to the treatment types and add a column with treatment type and another with which treatment did the patient recieve

proof_genes_long <- 
  proof_genes_pt.cox %>%
  pivot_longer(
    cols = c(HORMONE, CHEMO, SURGERY), # Add any other parameters here
    names_to = "Parameter",
    values_to = "Value"
  )

# 4.3.2 Plot using the new Value column for X and Parameter for facets

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

# 4.3.3 Kruskal test


kruskal.test(SCORE ~  PAM50, data = proof_genes_pt.cox)

# 4.4 Multivariate cox with treatment types

coxph(surv_obj ~ PAM50 * SCORE, data = proof_genes_pt.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)



# 5.- Resampling ----------------------------------------------------------


set.seed(456) 

# 5.1 Fold 20 times the whole data

final_resamples <- mc_cv(
  proof_genes_pt, 
  prop = 0.8, 
  times = 20, 
  strata = EVENT_STAT
)

# 5.2 Apply the final workflow to the resamples

final_resample_results <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = metric_set(concordance_survival)
  )

# 5.2.2 Metrics

collect_metrics(final_resample_results)

# 5.3 Define the metrics with specific time points to resample for AUC

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 5.3 .2 Apply the workflow to the resamples with respect to the AUC

set.seed(123)

res_auc_5y <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = survival_metrics,
    eval_time = c(36, 60, 120) # 5 and 10 years
  )

collect_metrics(res_auc_5y)

# 5.4 Make a df with the resample metrics

resamples <- bind_rows(final_resample_results$.metrics)


# 5.4.2 Again create data frame but this time with the distinct areas under the curve 

resamples_auc <- res_auc_5y %>%
  collect_metrics(summarize = FALSE) %>% 
  na.omit()

# 5.4.2.2 Create object wiuth mean, standard deviation, 95% confidence interval and standard error of the different time points

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

# 5.5 Add the fold numbers

resamples$fold <- final_resample_results$id

# 5.6 Plot The c score obtained at every fold and add a line with the median of the c scores obtained with the 20 folds

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


# 5.6.2 Plot the mean areas under the curve with its confidence intervals

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



ggplot(resamples_auc, aes(x = reorder(id, .estimate), y = .estimate)) + # Reorder so that it gives an increasing graph
  geom_point() +
  labs(x = "Folds", y = "Estimate", title = "Refold: C-score") +
  facet_wrap(~ .eval_time, ncol = 1)



# 5.7 From the fold outlier identify different clinical characteristics

for (i in c(17, 14, 7, 15)) { # Here add the number of the folds to analyze
  
  refold <- final_resamples$splits[[i]] # Gives the split of each fold
  
  # Patients used for training in fold i
  
  train_ids <- analysis(refold)
  
  # Patients used for validation in fold i
  
  test_ids <- assessment(refold)
  
  # Uses the metadata and obtains only the patients in i folds test set
  
  meta_testid <- ml_metadata[ml_metadata$PATIENT_ID %in% rownames(test_ids),]
  
  # Print the counts of desired characteristics of patrients in the i fold test data
  
  print(
    test_ids %>% 
      tibble::rownames_to_column("PATIENT_ID") %>% 
      dplyr::left_join(
        meta_testid,
        by = "PATIENT_ID",
        suffix = c("", ".drop")
      ) %>% 
      dplyr::select(-dplyr::ends_with(".drop")) %>% 
      dplyr::group_by(CLAUDIN_SUBTYPE, EVENT_STAT) %>% 
      dplyr::summarise(
        mean_OS = mean(EVENT_MON, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"   # removes all grouping from the output
      )
  )
  
  
  print(
    test_ids %>% 
      tibble::rownames_to_column("PATIENT_ID") %>% 
      dplyr::left_join(
        meta_testid,
        by = "PATIENT_ID",
        suffix = c("", ".drop")
      ) %>% 
      dplyr::select(-dplyr::ends_with(".drop")) %>% 
      dplyr::group_by(EVENT_STAT) %>% 
      dplyr::summarise(
        mean_OS = mean(EVENT_MON, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"   # removes all grouping from the output
      )
  )
  
  # Aditionally print the names and id of patient identified as high bias of the train set
  
  print(paste0(metadata.ER_POS_SURV$CLAUDIN_SUBTYPE[metadata.ER_POS_SURV$PATIENT_ID %in% top_bias_ids[top_bias_ids %in% rownames(test_ids)]],
               " ",
               top_bias_ids[top_bias_ids %in% rownames(test_ids)])
  )
  
  
  test_pred_refold <- predict(final_fit, new_data =  test_ids, type = "linear_pred")
  test_ids$.pred_linear_pred <- test_pred_refold$.pred_linear_pred
  
  time_roc_refolds <- timeROC(
    T = test_ids$EVENT_MON,
    delta = test_ids$EVENT_STAT,
    marker = - test_ids$.pred_linear_pred,
    cause = 1,
    times = c(36, 60, 120),  # 1y, 3y, 5y, 6y, 10y
    iid = TRUE
  )
  
  plot_roc_bias <- map_df(c(36, 60, 120), function(i){
    
    time_label <- i
    
    data.frame(
      
      FP = time_roc_refolds$FP[,paste0("t=", i)],
      
      TP = time_roc_refolds$TP[,paste0("t=", i)],
      
      Time = factor(i)
    )
    
  })
  
  
  legend_labels_bias <- 
    paste0("t=", time_roc_refolds$times, " (AUC: ", 100 * round(time_roc_refolds$AUC, 3), "%)")
  
  facet_labels_bias <- 
    data.frame(
      Time = factor(c(36, 60, 120)),
      AUC_Text = paste0("AUC: ", round(100 * as.numeric(time_roc_refolds$AUC[1:3]), 3), "%")
    )
  
  p1 <- ggplot(data = plot_roc_bias, aes(x = FP, y = TP, color = Time)) +
    geom_line(linewidth = 1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs( 
      title = paste0("Time-Dependent ROC Curves for fold ", i),
      subtitle = "Comparing model performance across different horizons",
      x = "1 - Specificity (FP)",
      y = "Sensitivity (TP)",
      color = "Time Point"
    ) + 
    scale_color_viridis_d(labels = legend_labels_bias) 
  
  print(p1)
  
  p2 <- ggplot(data = plot_roc_bias, aes(x = FP, y = TP)) +
    geom_line(color = "darkblue", linewidth = 1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ Time, ncol = 2) + 
    theme_bw() +
    labs(
      title = paste0("Time-Dependent ROC Curves for fold ", i),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    geom_text(data = facet_labels_bias, 
              aes(x = 0.75, y = 0.1, label = AUC_Text), 
              size = 4, fontface = "bold")
  
  print(p2)
  
}



# 6.- Observe distributions and shapiro -----------------------------------


# 6.1 Shapiro wilk test of desired distribution

shapiro.test(train_data2$risk_score)

# 6.2 Plot distribution

ggplot(train_data2, aes(x = risk_score)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "#81dcff", size = 1) +
  theme_minimal()


# 6.3 Plot distribution based on event stat

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

# 6.4 Wilcox test of desired data

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

# 6.4.2 Plot the  comparisons

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::deep_purple_material") + 
  facet_wrap(~ EVENT_STAT, 
             labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased")))


# 6.4.3 Prepare to plot the different comparisons

proof_genes_pt.long<- 
  proof_genes_pt.cox %>% 
  pivot_longer(
    cols = c(SURGERY, CHEMO, HORMONE), 
    names_to = "Parameter",
    values_to = "Value"
  )


# 6.4 Wilcox test of desired data

results_wilcox_treat <- proof_genes_pt.long %>%
  dplyr::select(EVENT_STAT, Parameter, Value, SCORE) %>% 
  group_by(Parameter, Value) %>%
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "BH"))

print(results_wilcox_treat$adj_p_value)

# 6.4.3.2 Plot

ggplot(data = proof_genes_pt.long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() + 
  facet_wrap(~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2,  labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) + 
  scale_fill_paletteer_d("palettesForR::Pastels") + 
  theme_gray(base_size = 18)
