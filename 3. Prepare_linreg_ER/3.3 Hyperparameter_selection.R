library(tidymodels)
library(censored) 
library(doParallel)
library(parallel)

# In this script we do refolds based on Monte Carlo cross validation so as to select the hyperparameters to be used in the regression

# 1.- Initial Tidymoedls workflow -----------------------------------------

# 1.1 Recipe using proof_genes_pt, so all the patients of METABRIC ER+

lr_rec <- recipe(surv_obj ~ ., data = proof_genes_pt) %>% 
  update_role(EVENT_MON, EVENT_STAT, new_role = "non_predictor") %>%
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% 
  step_nzv(all_predictors()) %>% 
  step_normalize(all_predictors()) %>% 
  step_spatialsign()

# 1.2 Model with tune parameters

lr_mod <- proportional_hazards(
  penalty = tune(), 
  mixture = tune()
) %>%
  set_engine("glmnet")

# 1.3 Workflow 

lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


# 2.- Define folds and tune grid ------------------------------------------


# 2.1 Define 20 resamples with Monte Carlo validation (no replacement)

set.seed(123)
resamples_20 <- mc_cv(
  proof_genes_pt, 
  prop = 0.8, 
  times = 20, 
  strata = EVENT_STAT
)

# 2.2 Grid Definition

grid <- grid_regular(
  penalty(range = c( - 4, 1)),   
  mixture(range = c(0, 1)),
  levels = 10
)


# 2.3.1 Set metrics to evaluate with tune grid levels

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 2.3.2 Run grid tuning 

res_nested <- tune_grid(
  lr_wf,
  resamples = resamples_20,
  grid = grid,
  metrics = survival_metrics,
  control = control_grid(save_pred = TRUE),
  eval_time = c(36, 60, 120) 
  
)


# 3.- Tuning results ------------------------------------------------------

# 3.1 Tuning results

final_metrics <- collect_metrics(res_nested)

# 3.2 Select best parameter for c-score

best_params <- select_best(res_nested, metric = "concordance_survival")

# 3.3  Show best parameters for AUC at 3 time points

for (i in c(36, 60, 120)) {
  best_params_auc <- show_best(res_nested, metric = "roc_auc_survival", n = 20, eval_time = i)
  print(best_params_auc)
}



# 4.- Complete model on resamples -----------------------------------------


# 4.1 Finalize workflow with selected best parameter

final_wf <- finalize_workflow(lr_wf, best_params)

# 4.2 Fit to see individual split results

final_resample_results <- final_wf %>%
  fit_resamples(
    resamples = resamples_20,
    metrics = metric_set(concordance_survival)
  )

collect_metrics(final_resample_results)


# 4.3 Define the metrics with specific time points to resample for AUC

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 4.3.2 Apply the workflow to the resamples with respect to the AUC

set.seed(123)

res_auc_5y <- final_wf %>%
  fit_resamples(
    resamples = resamples_20,
    metrics = survival_metrics,
    eval_time = c(36, 60, 120) 
  )

collect_metrics(res_auc_5y)

# 4.4 Make a df with the resample metrics for c score

resamples <- bind_rows(final_resample_results$.metrics)


# 4.4.2 Again create data frame but this time with the distinct areas under the curve 

resamples_auc <- res_auc_5y %>%
  collect_metrics(summarize = FALSE) %>% 
  na.omit()

# 4.4.2.2 Create object wiuth mean, standard deviation, 95% confidence interval and standard error of the different time points

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

# 4.5 Add the fold numbers

resamples$fold <- final_resample_results$id

# 4.6 Plot The c score obtained at every fold and add a line with the median of the c scores obtained with the 20 folds

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


# 4.6.2 Plot the mean areas under the curve with its confidence intervals

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

# 4.6.3 Plot AUC refolds with increasing AUC

ggplot(resamples_auc, aes(x = reorder(id, .estimate), y = .estimate)) + # Reorder so that it gives an increasing graph
  geom_point() +
  labs(x = "Folds", y = "Estimate", title = "Refold: AUC") +
  facet_wrap(~ .eval_time, ncol = 1)




