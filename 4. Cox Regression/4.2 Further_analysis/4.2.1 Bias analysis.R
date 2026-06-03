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


# 2.- Brier score, Schonfeild and martingale residuals --------------------

# 1 Add a row ID so we can find these patients later

id <- train_data %>%
  rownames_to_column("PATIENT_ID") %>%
  mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))

# 1.2 Utilize the fitted object to make predictions based on time

model_diagnostics <- augment(
  final_fit,
  new_data = id,
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)


# 2.1 Brier score

model_diagnostics %>% 
  brier_survival(truth = surv_obj, .pred)


# 2.2 Martingale and Schofeild residuals

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



# 5.- Resampling ----------------------------------------------------------


set.seed(456) 

# 5.1 Fold 100 times the whole data

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
    eval_time = c(12, 36, 60, 72, 120),  
    control = control_resamples(save_pred = TRUE) 
  )

# 5.3.1.2 C score

final_resample_cscore <- final_wf %>%
  fit_resamples(
    resamples = final_resamples,
    metrics = metric_set(concordance_survival),
    control = control_resamples(save_pred = TRUE) 
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

# 5.6.1 Summary for C-score

summary_cscore <- resamples_cscore %>% 
  filter(.metric == "concordance_survival") %>%
  summarise(
    mean_c    = mean(.estimate),
    median_c  = median(.estimate),
    sd_c      = sd(.estimate),  
    n         = n(),
    q2.5  = quantile(.estimate, probs = 0.025),
    q97.5  = quantile(.estimate, probs = 0.975),
    p_value   = mean(.estimate <= 0.5)
  )

# 5.6.2  Summary for AUC

summary_auc <- resamples_auc %>%
  group_by(.eval_time) %>%
  summarise(
    mean_auc   = mean(.estimate),
    median_auc = median(.estimate),
    sd_auc     = sd(.estimate), 
    n          = n(),
    q2.5   = quantile(.estimate, probs = 0.025),
    q97.5   = quantile(.estimate, probs = 0.975),
    p_value    = mean(.estimate <= 0.5),
    .groups = "drop"
  )

# 5.7 Plot the distribution of the C-index and the AUC at the studied time points with the intervals determined by the quantiles

mc_plot1 <- 
  ggplot(resamples_c, aes(x = reorder(id, .estimate), y = .estimate)) + # Reorder so that it gives an increasing graph
  geom_point() +
  geom_hline(yintercept = mean(resamples_c$.estimate), linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = mean(resamples_c$.estimate),
           label = paste0("Median = ", round(median(resamples_c$.estimate), 3)),
           vjust = -1) + 
  labs(x = "Resample split (ordered by performance)", y = "Concordance score", title = "Refold: C-score Distribution (100 iterations) Recurrence")+
  annotate("text", x = 5.2, y = min(resamples_c$.estimate),
           label = paste0("Min = ", round(min(resamples_c$.estimate), 3)),
           vjust = 0.4) +
  annotate("text", x = length(resamples_c$.estimate) - 4, y = max(resamples_c$.estimate),
           label = paste0("Max = ", round(max(resamples_c$.estimate), 3)),
           vjust = 0) +
  scale_x_discrete(breaks = seq(0, 100, by = 100)) +
  theme_minimal(base_size = 20)


mc_plot2 <- 
  ggplot(summary_auc, aes(x = factor(.eval_time), y = median_auc)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
  geom_text(aes(label = paste0("Median = ", round(median_auc, 3))), vjust = 0, size = 5, hjust = - 0.2) +
  labs(x = "Time (Months)", y = "AUC", title = "Refold: Time-dependent AUC (median and quantile 2.5 and 97.5)") +
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

# 6.4.0 Convert EVENT_sTAT to factor

proof_genes_pt.cox <- 
  proof_genes_pt.cox %>%
  mutate(EVENT_STAT = factor(EVENT_STAT)) 

# 6.4 Wilcox test of Histology, Claudin subtype and INTCLUST data

ch_var <- c("HIST", "INTCLUST", "PAM50")

for (i in ch_var) {

print(proof_genes_pt.cox %>%
  dplyr::select(EVENT_STAT, all_of(i), SCORE) %>% 
  group_by(.data[[i]]) %>%
  filter(n_distinct(EVENT_STAT) == 2) %>%
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "holm"))
)

}

# 6.4.2 Plot the  comparisons


for (i in ch_var) {
  

  score_subtype <- ggplot(
    proof_genes_pt.cox,
    aes(x = .data[[i]], y = SCORE, fill = .data[[i]])
  ) +
    geom_boxplot() +
    scale_fill_paletteer_d("Redmonder::dPBIPuOr") +
    facet_wrap(
      ~ EVENT_STAT,
      labeller = labeller(EVENT_STAT = c("0" = "Alive",
                                         "1" = "Deceased"))
    ) +
    theme_classic(base_size = 22) + 
    ggtitle("Score change on subtype: METABRIC") +
    theme(
      strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
            strip.text = element_text(color = "white", face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    labs(y = "", 
         x = "",
         tag = "A")
  
  print(score_subtype)
}



# 6.4.3 Prepare to plot the different comparisons

proof_genes_pt.long<- 
  proof_genes_pt.cox %>% 
  pivot_longer(
    cols = c(SURGERY, CHEMO, HORMONE), 
    names_to = "Parameter",
    values_to = "Value"
  ) %>% 
  mutate(Value = case_when(
    toupper(Value) == "YES" ~ "Yes",
    toupper(Value) == "NO"  ~ "No",
    Value == "BREAST CONSERVING" ~ "Breast conserving",
    Value == "MASTECTOMY" ~ "Matectomy", 
    TRUE ~ Value 
  ))


# 6.4 Wilcox test between treatment types

proof_genes_pt.long %>%
  dplyr::select(EVENT_STAT, Parameter, Value, SCORE) %>% 
  group_by(Parameter, Value) %>%
  filter(n() != 1) %>% 
  group_modify(~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "holm"),
         )

# 6.5 Cox analysis of interaction between score and treatment 

tx_vars <- c("CHEMO", "HORMONE", "SURGERY")

for (i in tx_vars) {
  
  formula <- as.formula(paste("surv_obj ~ ", i, " * SCORE")) # Establish formula
  
  proof_genes_pt.txcox <- 
    proof_genes_pt.cox %>% 
    group_by(.data[[i]]) %>%
    filter(n_distinct(EVENT_STAT) == 2, # FIlter data with complete separation
           !(is.na(.data[[i]])), # Filter NA
           !(.data[[i]] == "") # Filter fod the unregistered surgery
    ) %>% 
    ungroup()
  
  cox_sum <- summary(coxph(formula , data = proof_genes_pt.txcox))
  print(cox_sum)
  print(i)
  
  
}


# 6.4.3.2 Plot

score_tx <- ggplot(data = proof_genes_pt.long, aes(y = SCORE, x = Value, fill = Value)) +
  geom_boxplot() + 
  facet_wrap(~ Parameter + EVENT_STAT, scales = "free_x", ncol = 2,  labeller = labeller(Parameter = c("CHEMO" = "Chemotherapy", "HORMONE" = "Hormone therapy", "SURGERY" = "Surgery Modality"),
                                                                                         EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) + 
  scale_fill_paletteer_d("khroma::iridescent", direction = - 1) + 
  theme_classic(base_size = 15) + 
  ggtitle("Score change on treatment: METABRIC") +
  theme(
    strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  labs(tag = "A",
       x = "")

