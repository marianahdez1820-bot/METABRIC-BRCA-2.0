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
  rownames_to_column("PATIENT_ID") %>%
  mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))

# 2.2 Utilize the fitted object to make predictions based on time

model_diagnostics <- augment(
  final_fit,
  new_data = id,
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)

# 2.2.2 Create list to then add the results

list <- list()

no_img <- 1

#> 2.3 For loop that at each desired time point calculates the bias scores, extracts patients with high bias scores
#> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time

for (i in c(36, 60, 120)) {

  eval_time <- i

  # 2.3 Unnest the predictions and find the biggest outliers on a set point and event

  outliers <- model_diagnostics %>%
    dplyr::select(PATIENT_ID, EVENT_MON, EVENT_STAT, .pred, .pred_time) %>%
    unnest(.pred) %>%
    filter(.eval_time == eval_time) %>%
    arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive


  # 2.4 Object with metadata and score characteristics

  outlier_summary <- outliers %>%
    inner_join(ml_metadata, by = "PATIENT_ID", suffix = c("", ".drop")) %>%
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
      CELLULARITY,
      HISTOLOGICAL_SUBTYPE,
      .pred_time,
      CHEMOTHERAPY,
      RADIO_THERAPY,
      HORMONE_THERAPY,
      BREAST_SURGERY,
      AGE_AT_DIAGNOSIS,
      OS_STATUS,
      OS_MONTHS,
      RFS_MONTHS,
      RFS_STATUS
    )


  # 2.5 Create bias score for defined time

  outliers_bias <- outlier_summary %>%
    mutate(
      bias_score = (((1 - EVENT_STAT) - .pred_survival) ^ 2) * (EVENT_MON - .eval_time) * ( 1 - (2 * EVENT_STAT))
    ) %>%
    arrange(desc(bias_score))

  # 2.6 Identify the highest bias patients

  # 2.6.1 Obtain mean and sd and then filter baed on patients higher than determined SD

  extreme_outliers <-
    outliers_bias %>%
    mutate(mean_bias = mean(bias_score),
           sd_bias = sd(bias_score)) %>%
    filter(bias_score > (mean_bias + (1 * sd_bias)))

  # 2.6.2 Obtain their IDs

  top_bias_ids <- extreme_outliers$PATIENT_ID

  print(length(top_bias_ids))

  # 2.6.3 Identify different characteristics of patients identified as top bias

  print(outliers_bias %>%
          filter(PATIENT_ID %in% top_bias_ids) %>%
          group_by(RADIO_THERAPY) %>%
          summarise(
            count = n(),
            avg_pred_event = mean(.pred_survival),
            avg_event_time = mean(EVENT_MON)
          ) %>%
          arrange(desc(count)),

        n = 200
          )

  # 2.7.1 Add a column identifying patients as top bias or not

  outliers <-
    outliers %>%
    mutate(quadrant = case_when(
             PATIENT_ID %in% top_bias_ids & EVENT_STAT == 0 ~ 2,
             PATIENT_ID %in% top_bias_ids & EVENT_STAT == 1 ~ 1,
             TRUE ~ 0
           ))

  print(outliers %>%
          group_by(EVENT_STAT) %>%
          dplyr::count(quadrant))

  # 2.7.2 Plot
  if(no_img == 0){

  theme_embedded <- theme_linedraw(base_size = 25) +
    theme(
      legend.position = c(0.95, 0.3), # Adjust coordinates (x, y) from 0 to 1
      legend.background = element_rect(fill = alpha("white", 0.5))
    )

  # 2.7.2.1 Plot colored by EVENT_STAT

  p1 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(EVENT_STAT), shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
    stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    scale_color_viridis_d() +
    labs(
      title = "Event Status Distribution",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival",
      color = "Event Stat",
      shape = "Event Stat"
    ) +
    theme_embedded

  # 2.7.2.2 Plot colored by quadrant

  p2 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(quadrant), shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) +
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    scale_color_viridis_d() +
    labs(
      title = "Quadrant Analysis",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival",
      color = "Quadrant",
      shape = "Event Stat"
    ) +
    theme_embedded

  # 2.7.3 Combine and stack

  print( p1 + p2)
  }else{

  }

  list[[i]] <- top_bias_ids

}

# 2.8 Obtain patients found on all of the iterations of the for loop as top bias patients

bias_interesct <- intersect(intersect(list[[36]], list[[60]]) , list[[120]])

# 2.8.2 Similar but all unique so even if they appear once we register them

bias_diff_id <- unique(c(list[[36]], list[[60]], list[[120]]))

# 2.9 Observe characteristics of the bias_diff_id patients

print(outliers_bias %>%
  filter(PATIENT_ID %in% bias_diff_id) %>%
  mutate(intcluster = ifelse(INTCLUST == "3" | INTCLUST == "4ER+" | INTCLUST == "7" | INTCLUST == "8",
                            "Low",
                            "High"),
         lymph_group = case_when(
           LYMPH_NODES_EXAMINED_POSITIVE == 0 ~ 0,
           LYMPH_NODES_EXAMINED_POSITIVE > 0 & LYMPH_NODES_EXAMINED_POSITIVE  < 4 ~ 1,
           LYMPH_NODES_EXAMINED_POSITIVE >= 4 ~ 2
         ),
         any_tx = ifelse(HORMONE_THERAPY == "YES" | RADIO_THERAPY == "YES" | CHEMOTHERAPY == "YES",
                         "Recieved",
                         "Not recieved"),
         non_hm_tx = case_when(
           RADIO_THERAPY == "YES" & CHEMOTHERAPY == "YES" & HORMONE_THERAPY == "YES" ~ "All",
           xor(RADIO_THERAPY == "YES", CHEMOTHERAPY == "YES") & HORMONE_THERAPY == "YES" ~ "Hormone and else",
           RADIO_THERAPY == "NO" & CHEMOTHERAPY == "NO" & HORMONE_THERAPY == "YES" ~ "Only hormone",
           RADIO_THERAPY == "NO" & CHEMOTHERAPY == "NO" & HORMONE_THERAPY == "NO" ~ "Nothing",
           xor(RADIO_THERAPY == "YES", CHEMOTHERAPY == "YES") & HORMONE_THERAPY == "NO" ~ "Else no hormone",
           RADIO_THERAPY == "YES" & CHEMOTHERAPY == "YES" & HORMONE_THERAPY == "NO" ~ "Both else no hormone"
           ),
         age_bin = ifelse(AGE_AT_DIAGNOSIS > 50,
                          ">50",
                          "<50")
         ) %>%
  group_by(EVENT_STAT, HORMONE_THERAPY, RADIO_THERAPY, CHEMOTHERAPY, BREAST_SURGERY, intcluster, age_bin) %>%
  summarise(
    count = n(),
    avg_bias = mean(bias_score),
    avg_pred_event = mean(.pred_survival),
    avg_event_time = mean(EVENT_MON),
    avg_pred_time = mean(.pred_time)
  ) %>%
  arrange(desc(EVENT_STAT)),
n = 200)

# 2.10 Bias scoEVENT_STAT# 2.10 Bias score based on pred time not on tiHISTOLOGICAL_SUBTYPE# 2.10 Bias scoEVENT_STAT# 2.10 Bias score based on pred time not on time evaluated

outlier_pred_time <- model_diagnostics %>% 
  unnest(.pred) %>% 
  dplyr::select(- all_of(proof_genes)) %>% 
  group_by(PATIENT_ID) %>% 
  mutate(bias_score = # Notice the mean in .pred sruvival is the mean of each oatients prediction at all 3 time points since in this case we dont care about time point evaluated but on time point predicted so we also care about the global % predicted 
          (((1 - EVENT_STAT) - mean(.pred_survival)) ^ 2) * (EVENT_MON - .pred_time) * ( 1 - (2 * EVENT_STAT))
  ) %>%
  ungroup() %>% 
  arrange(desc(bias_score))

# 2.10.2 Identify the highest bias patients

# 2.10.3 Obtain mean and sd and then filter based on patients higher than determined SD

extreme_outliers_pred_time <-  
  outlier_pred_time %>% 
  mutate(mean_bias = mean(bias_score),
         sd_bias = sd(bias_score)) %>% 
  filter(bias_score > (mean_bias + 1 * sd_bias))

# 2.10.3 Obtain their IDs

top_bias_ids_pred_time <- extreme_outliers_pred_time$PATIENT_ID


# 2.10.4 Add a column identifying patients as top bias or not

outlier_pred_time <- 
  outlier_pred_time %>% 
  mutate(quadrant = case_when(
    PATIENT_ID %in% top_bias_ids & EVENT_STAT == 0 ~ 2,
    PATIENT_ID %in% top_bias_ids & EVENT_STAT == 1 ~ 1,
    TRUE ~ 0
  ))



# 2.10.5 Plot prepare theme

theme_embedded <- theme_linedraw(base_size = 25) + 
  theme(
    legend.position = c(0.5, 0.1), # Adjust coordinates (x, y) from 0 to 1
    legend.background = element_rect(fill = alpha("white", 0.5))
  )

# 2.10.5.1 Plot colored by EVENT_STAT

p1_global <- ggplot(outlier_pred_time, aes(x = EVENT_MON, y = .pred_survival, color = factor(EVENT_STAT), shape = factor(EVENT_STAT))) +
  geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
  stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
  geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + 
  labs(
    title = "Event Status Distribution",
    x = paste0("Actual Event Time (Months)", eval_time),
    y = "Predicted Survival",
    color = "Event Stat",
    shape = "Event Stat"
  ) +
  theme_embedded

p2_global <- ggplot(outlier_pred_time, aes(x = EVENT_MON, y = .pred_time, color = factor(quadrant), shape = factor(EVENT_STAT))) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) + 
  geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + 
  labs(
    title = "Quadrant Analysis",
    x = paste0("Actual Event Time (Months)", eval_time),
    y = "Predicted time",
    color = "Quadrant",
    shape = "Event Stat"
  ) +
  theme_embedded

p1_global + p2_global

# 2.11 Brier score

model_diagnostics %>% 
  brier_survival(truth = surv_obj, .pred)


# 2.12 Martingale and Schofeild residuals

cox.zph(cox_model)

ggcoxzph(cox.zph(cox_model))

ggcoxdiagnostics(cox_model, type = "martingale",
                 linear.predictions = FALSE, ggtheme = theme_bw())

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
  times = 100, 
  strata = EVENT_STAT
)


# 5.2 Define metrics with specific evaluation times 

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 5.3 Run the resamples (we run for c score and for auc so as to then calculate their individual CI and p val)

# 5.3.1.1 AUC and C score 

final_resample_results <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = survival_metrics,
    eval_time = c(12, 36, 60, 72, 120),  # Required for time-dependent ROC
    control = control_resamples(save_pred = TRUE) # Useful if you need predictions later
  )

# 5.3.1.2 C score

final_resample_cscore <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = metric_set(concordance_survival),
    control = control_resamples(save_pred = TRUE) # Useful if you need predictions later
  )


# 5.3.2 Obtain metrics

collect_metrics(final_resample_results)


# 5.4 Extract unsummarized metrics for every fold

resamples_all <- collect_metrics(final_resample_results, summarize = FALSE)

resamples_cscore <- collect_metrics(final_resample_cscore, summarize = FALSE)

# 5.5.1 Extract the C-scores

resamples_c <- resamples_all %>% 
  filter(.metric == "concordance_survival")

# 5.5.2 Extract the AUC

resamples_auc <- resamples_all %>% 
  filter(.metric == "roc_auc_survival") %>% 
  na.omit()


# 5.6.1 Tables with CI 95%, z statistic and p val for C score

summary_cscore <- resamples_cscore  %>% 
  filter(.metric == "concordance_survival")%>%
  summarise(
    mean = mean(.estimate),
    sd = sd(.estimate),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    z_stat = (mean - 0.5) / se,
    p_value = 1 - pnorm(z_stat)
  )


# 5.6.2 Tables with CI 95%, z statistic and p val for AUC

summary_auc <- resamples_auc %>%
  group_by(.eval_time) %>%
  summarise(
    mean = mean(.estimate),
    sd = sd(.estimate),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    z_stat = (mean - 0.5) / se,
    p_value = 1 - pnorm(z_stat),
    .groups = "drop"
  )

max_fold <- rownames(resamples_c)[resamples_c$.estimate == max(resamples_c$.estimate)]

min_fold <- rownames(resamples_c)[resamples_c$.estimate == min(resamples_c$.estimate)]

# 5.7 From the fold outlier identify different clinical characteristics

for (i in c(as.numeric(min_fold), as.numeric(max_fold))) { # Here add the number of the folds to analyze

  refold <- final_resamples$splits[[i]] # Gives the split of each fold
  
  # Patients used for training in fold i
  
  train_ids <- analysis(refold)
  
  # Patients used for validation in fold i
  
  test_ids <- assessment(refold)
  
  # Uses the metadata and obtains only the patients in i folds test set
  
  meta_testid <- ml_metadata[ml_metadata$PATIENT_ID %in% rownames(test_ids),]
  
  test_ids$score <- predict(final_fit, new_data = test_ids, type = "linear_pred")$.pred_linear_pred
  
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
      dplyr::group_by(HISTOLOGICAL_SUBTYPE, CLAUDIN_SUBTYPE, INTCLUST, EVENT_STAT) %>% 
      dplyr::summarise(
        mean_OS = mean(EVENT_MON, na.rm = TRUE),
        n = dplyr::n(),
        mean_score = mean(score),
        .groups = "drop"   
      ),
    n = 250
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
        .groups = "drop"   
      )
  )
  
  # Aditionally print the names and id of patient identified as high bias of the train set
  
  print(paste0(ml_metadata$CLAUDIN_SUBTYPE[ml_metadata$PATIENT_ID %in% bias_diff_id[bias_diff_id %in% rownames(test_ids)]],
               " ",
               bias_diff_id[bias_diff_id %in% rownames(test_ids)])
  )
  
  
  test_pred_refold <- predict(final_fit, new_data =  test_ids, type = "linear_pred")
  test_ids$.pred_linear_pred <- test_pred_refold$.pred_linear_pred
  
  
}


mc_plot1 <- 
  ggplot(resamples_c, aes(x = reorder(id, .estimate), y = .estimate)) + # Reorder so that it gives an increasing graph
  geom_point() +
  geom_hline(yintercept = mean(resamples_c$.estimate), linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = mean(resamples_c$.estimate),
           label = paste0("Mean = ", round(mean(resamples_c$.estimate), 3)),
           vjust = -1) + 
  labs(x = "Resample split (ordered by performance)", y = "Concordance score", title = "Refold: C-score Distribution (100 iterations)")+
  annotate("text", x = 5.2, y = min(resamples_c$.estimate),
           label = paste0("Min = ", round(min(resamples_c$.estimate), 3)),
           vjust = 0.4) +
  annotate("text", x = length(resamples_c$.estimate) - 4, y = max(resamples_c$.estimate),
           label = paste0("Max = ", round(max(resamples_c$.estimate), 3)),
           vjust = 0) +
  scale_x_discrete(breaks = seq(0, 100, by = 100)) +
  theme_minimal(base_size = 20)


mc_plot2 <- 
  ggplot(summary_auc, aes(x = factor(.eval_time), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(aes(label = paste0("Mean = ", round(mean, 3))), vjust = 0, size = 5, hjust = - 0.2) +
  labs(x = "Time (Months)", y = "AUC", title = "Refold: Time-dependent AUC (mean ± 95% CI)") +
  theme_minimal(base_size = 20)

mc_plot1 / mc_plot2

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
  theme_linedraw() +
  labs(fill = "Event",
       x = "Risk score",
       y = "Density",
       title = "Distribution of score with respect to event")

# Convert to factor

proof_genes_pt.cox <- 
  proof_genes_pt.cox %>%
  mutate(EVENT_STAT = factor(EVENT_STAT)) 

# 6.4 Wilcox test of desired data

proof_genes_pt.cox %>%
  dplyr::select(EVENT_STAT, PAM50, SCORE) %>% 
  group_by(PAM50) %>%
  filter(n_distinct(EVENT_STAT) == 2) %>%
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "holm"))



# 6.4.2 Plot the  comparisons

ggplot(proof_genes_pt.cox, aes(y = SCORE, x = HIST, fill = HIST)) +
  geom_boxplot() +
  scale_fill_paletteer_d("colorBlindness::Blue2Green14Steps") + 
  facet_wrap(~ EVENT_STAT, 
             labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) +
  theme_classic() + 
  geom_vline(xintercept = 8)


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
  filter(n() != 1) %>% 
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "BH"))

print(results_wilcox_treat$adj_p_value)

vars <- c("CHEMO", "HORMONE", "SURGERY")

for (i in vars) {
  formula <- as.formula(paste("surv_obj ~ ", i, " * SCORE"))
  
  proof_genes_pt.txcox <- 
    proof_genes_pt.cox %>% 
    group_by(.data[[i]]) %>%
    filter(n_distinct(EVENT_STAT) == 2,
           !(is.na(.data[[i]])),
           !(.data[[i]] == "")
           ) %>% 
    ungroup()
  cox_sum <- summary(coxph(formula , data = proof_genes_pt.txcox))
  print(cox_sum)
  print(i)
print(paste0("HR of ", round(cox_sum$coefficients[3, 1], 2), " (Ci 95% ", round(cox_sum$conf.int[3, 3], 2), " - ", round(cox_sum$conf.int[3, 4], 2), ", pval ", cox_sum$coefficients[3, 5], ")")
)  
  
  
}

# 6.4.3.2 Plot

ggplot(data = proof_genes_pt.long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() + 
  facet_wrap(~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2,  labeller = labeller(Parameter = c("CHEMO" = "Chemotherapy", "HORMONE" = "Hormone therapy", "SURGERY" = "Surgery Modality"),
                                                                                         EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) + 
  scale_fill_paletteer_d("colorBlindness::Blue2DarkOrange12Steps") +
  theme_classic(base_size = 32)

