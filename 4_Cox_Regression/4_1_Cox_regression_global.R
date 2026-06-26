library(tidymodels)
library(censored)
library(survminer)
library(broom)
library(survival)
library(timeROC)
library(flextable)

# InitiaRecurrence# Initial object late_genes.patients created in preprocessing
# Its made so that the modifications have to be done on preprocessing so even if the
# thing to be studied is Recurrence the object will stay as EVENT and EVENT_MON


# 2.- Preparing recipe and model ------------------------------------------

# 2.1 Recipe

lr_rec <- recipe(surv_obj ~ ., data = train_data) %>% # Recurrence object created in preprocessing for linear regression
  update_role(EVENT_MON, EVENT_STAT, new_role = "non_predictor") %>%
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% # Eliminates variables with a single value
  step_nzv(all_predictors()) # Eliminates highly sparse variables


# 2.2 Model

lr_mod <- proportional_hazards(
  penalty = tune(),    # lambda establishes the severith of the penalty
  mixture = tune()     # alpha establishes the type, 1 being lasso, 0 being ridge, and 0.5 being elasticnet
) %>%
  set_engine("glmnet",
             cox.ties = "breslow") # Engine that permits penalizing by elasticnet, ridge, and lasso

# 2.3 Workflow

lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


# 3.- Selecting best penalizing parameters ------------------------------------

if(!exists("best_params")){
  
  set.seed(123)
  
  # 3.1 Folds for evaluating with resamples
  
  folds <- vfold_cv(
    train_data,
    v = 10,
    strata = EVENT_STAT
  )
  
  # 3.2 Grid for penalizing range
  
  grid <- grid_regular(
    penalty(range = c( - 4, 1)),   
    mixture(range = c(0, 1)),
    levels = 10
  )
  
  # 3.3 Running the different penalization methods
  
  res_ml <- tune_grid(
    lr_wf,
    resamples = folds,
    grid = grid,
    metrics = metric_set(concordance_survival), # Evaluates the different penalizing methods by c score
    control = control_grid(save_pred = TRUE)
  )
  
  # 3.3.1 Observe metrics
  
  collect_metrics(res_ml)
  
  # 3.3.2 Object with best parameters for penalizing
  
  best_params <- select_best(res_ml, metric = "concordance_survival")
  
}else{
  print("Parameter selected on MC cross validation")
}


# 4.- Actual training -----------------------------------------------------

# 4.1 Final workflow with the selecting the best parameter tested previously

final_wf <- finalize_workflow(lr_wf, best_params)

# 4.2 Final fit with training data

final_fit <- fit(final_wf, data = train_data)



# 4.2.2 Observing genes that are maintained after penalziation


coef_tbl <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  arrange(desc(abs(estimate)))

cat(coef_tbl$term, sep = ", ")

# train_data.final <- bake(extract_recipe(final_fit), new_data = train_data)

# 4.3 Predictions on train data

train_pred <- predict(final_fit, new_data = train_data, type = "linear_pred")

train_pred <- 
  train_pred %>% 
  as.data.frame()

# 4.3.2 Object with train data and its predicted scores

train_data2 <- 
  train_data %>% 
  mutate(risk_score = train_pred$.pred_linear_pred) 

# 4.4 Using those scores to calculate the cutpoint

true_cut <- 
  surv_cutpoint(
    data = train_data2, 
    time = "EVENT_MON", 
    event = "EVENT_STAT", 
    variables = "risk_score"
  )


median(train_data2$risk_score)

# 5.- Testing -------------------------------------------------------------

#test_data.final <- bake(extract_recipe(final_fit), new_data = test_data)

# 5.1 Predictions on test data

test_pred <- predict(final_fit, new_data = test_data, type = "linear_pred")

# 5.2 Creating groups for Kaplan-Meier curves

# 5.2.1 Creating column on test data with its prediction

test_data$risk_score <- test_pred$.pred_linear_pred


# 5.2.2 Dividing the groups by median so as to establish a high and low risk and create a column


test_data <- 
  test_data %>% 
  mutate(
    risk_group =  factor(ifelse(
      risk_score < true_cut$cutpoint[1,1],
      "High",
      "Low"
    )
    ),
    risk_group_median = factor(ifelse(
      risk_score < median(train_data2$risk_score),
      "High",
      "Low"
    )
    )
  )



test_data$risk_group <- relevel(test_data$risk_group, ref = "Low")

# 5.3 Creating kapan meier curve

# 5.3.1 Based on the Recurrence object compare the risk_groups created previously

fit_km <- survfit(surv_obj ~ risk_group, data = test_data)

# 5.3.2 Plot

ggsurvplot(fit_km,
           data = test_data,
           pval = TRUE, 
           risk.table = TRUE,
           
           title = "Recurrence ER+ METABRIC",
           ylab = "Recurrence probability",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 300),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#c380d3", "#ff89d4"),
)

# 5.4 Tables

survdiff(surv_obj ~ risk_group, data = test_data)

# 5.4.1 Cox regression analysis table with HR

summary_cox <- summary(coxph(surv_obj ~ risk_group, data = test_data))

# 6.- C score

# 6.1 Combine the test and training data to compare concordance

eval_df <- dplyr::bind_cols(
  test_data,
  test_pred
)

# 6.2 Compare concordance

concordance <- concordance(
  surv_obj ~ .pred_linear_pred,
  data = eval_df
)

c_index_summary <- data.frame(
  C_Index  = concordance$concordance,
  SE       = sqrt(concordance$var),
  data_set = "METABRIC"
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )

print(c_index_summary)



# 7.- Area under de curve (AUC) per time --------------------------------------------------------


# 7.1 Creating object on testing data

time_roc <- timeROC(
  T = test_data$EVENT_MON,
  delta = test_data$EVENT_STAT,
  marker = - test_pred$.pred_linear_pred,
  cause = 1,
  times = c(12, 36, 60, 72, 120),  # 1y, 3y, 5y, 6y, 10y
  iid = TRUE
)

# 2.6.2 Table with confidence interval and z stat and estimated p val

auc_ci <- data.frame(
  AUC  = time_roc$AUC,
  SE   = time_roc$inference$vect_sd_1,
  time = time_roc$times,
  data_set = "METABRIC"
) %>% 
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )


test_data <- 
  test_data %>% 
  mutate(EVENT_STAT = factor(EVENT_STAT))

global_roc <- roc_curve(test_data,
                        EVENT_STAT,
                        risk_score
) %>%
  mutate(label = "METABRIC")


# 7.2 Loop that creates data frame with true positive and falsa positives of each time point

plot_roc <- map_df(c(12, 36, 60, 72, 120), function(i){
  
  
  data.frame(
    
    FP = time_roc$FP[,paste0("t=", i)],
    
    TP = time_roc$TP[,paste0("t=", i)],
    
    Time = factor(i),
    
    data_set = "METABRIC"
  )
  
})


# 7.3.2 Object with labels for the plot with the numerical values of the AUCs

# 7.3.3 Same thing but for facet wrap labels

facet_labels <- 
  data.frame(
    Time = factor(c(12, 36, 60, 72, 120)),
    AUC_Text = paste0("AUC: ", round(100 * as.numeric(time_roc$AUC[1:5]), 3), "%"),
    data_set = "METABRIC"
  )



# 7.4.2 Plot of the ROC curve with facet wraping of the different time points

ggplot(data = plot_roc, aes(x = FP, y = TP)) +
  geom_line(color = "darkblue", linewidth = 1) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Time, ncol = 2) + 
  theme_bw() +
  labs(
    title = "ROC Curves by Time Point",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  geom_text(data = facet_labels, 
            aes(x = 0.75, y = 0.1, label = AUC_Text), 
            size = 4, fontface = "bold")


# 8.- Independence test ---------------------------------------------------

# 8.1 Keep only the patients of the test data

proof_genes_pt.cox <- proof_genes_pt[rownames(test_data),]

# 8.2 Create object with parameters to evaluate on the cox model

proof_genes_pt.cox <- 
  proof_genes_pt.cox %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(ml_metadata, by = "PATIENT_ID", suffix = c("", ".y")) %>% # Join with metadata and eliminate duplicates
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("PATIENT_ID") %>% 
  mutate(HER2 = HER2_SNP6, # Create objects for evaluation
         LYMPH = LYMPH_NODES_EXAMINED_POSITIVE,
         AGE = as.numeric(AGE_AT_DIAGNOSIS ),
         MENO = INFERRED_MENOPAUSAL_STATE    ,
         SCORE = test_data$risk_score, # This one is the score of the model
         HORMONE = HORMONE_THERAPY ,
         CHEMO = CHEMOTHERAPY,
         SURGERY = BREAST_SURGERY,
         PAM50 = CLAUDIN_SUBTYPE,
         NPI = NPI,
         HIST = HISTOLOGICAL_SUBTYPE,
         RADIO = RADIO_THERAPY
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                AGE,
                LYMPH,
                HER2,
                MENO,
                SCORE,
                HORMONE,
                CHEMO,
                SURGERY,
                PAM50,
                INTCLUST, 
                EVENT_STAT,
                EVENT_MON,
                NPI,
                HIST,
                RADIO) %>%  
  na.omit()

# 8.3 Multivariate cox comparing the variables including the score to test its independent value


coxph(surv_obj ~ AGE + LYMPH + SCORE, 
      data = proof_genes_pt.cox)

cox_model <- coxph(surv_obj ~ NPI + HORMONE + CHEMO + RADIO + strata(SURGERY) + MENO + strata(HER2) + AGE + SCORE +  PAM50 + strata(INTCLUST) + strata(HIST), 
                   data = proof_genes_pt.cox)

summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, 
              data = proof_genes_pt.cox))


summary(cox_model)

# 8.3.2 Tidy format

independent_prog <- cox_model %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

supplementary_table <- independent_prog %>%
  mutate(
    Feature = recode_values(term,
                            "SCORE"              ~ "Signature Score",
                            "LYMPH"              ~ "Lymph Node Status",
                            "AGE"                ~ "Age at Diagnosis",
                            "NPI"                ~ "Nottingham Prognostic Index",
                            "CHEMOYES"           ~ "Chemotherapy (Yes)",
                            "HORMONEYES"         ~ "Hormone Therapy (Yes)",
                            "MENOPre"            ~ "Menopausal Status (Pre)",
                            "PAM50LumA"          ~ "PAM50 Luminal A",
                            "PAM50LumB"          ~ "PAM50 Luminal B",
                            "PAM50Her2"          ~ "PAM50 Her2-enriched",
                            "PAM50claudin-low"   ~ "PAM50 Claudin-low",
                            "PAM50Normal"        ~ "PAM50 Normal-like",
                            "PAM50NC"        ~ "PAM50 Not classified",
                            "RADIOYES" ~ "Radiotherapy (Yes)",
                            default = term
    ),
    
    `Hazard Ratio (HR)` = round(estimate, 3),
    `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
    
    `p-value` = scales::scientific(p.value, digits = 3)
  ) %>%
  dplyr::select(Feature, `Hazard Ratio (HR)`, `95% Confidence Interval`, `p-value`) %>%
  arrange(`Hazard Ratio (HR)`)

flextable(supplementary_table) %>%
  autofit() 


# 9.- Results -------------------------------------------------------------

# 9.1 Index number of the parameters to print the evaluation and to graph

num_param_compare <- c(1, 2, 3, 4, 6, 7, 10)

# 9.2 Concatenate strings with the respective result

cat(paste0("Signature with ", length(coef_tbl$term), " genes (", paste(coef_tbl$term, collapse = ", "), ")"),
    paste0("The selected parameters were an alpha of ", best_params$mixture, " and a lambda of ", best_params$penalty),
    paste0("The signature got a C-score of ", concordance$concordance),
    paste0("HR of ", round(summary_cox$coefficients[2], 2), " (CI 95% of ", round(summary_cox$conf.int[3], 2), " - ", round(summary_cox$conf.int[4], 2), " pval ", summary_cox$coefficients[5], ")"),
    paste0("AUC at 1 year of ", round(time_roc$AUC[1], 2), " at 3 years of ", round(time_roc$AUC[2], 2), " at 5 years of ", round(time_roc$AUC[3], 2), " at 7 years of ", round(time_roc$AUC[4], 2), " and at 10 years of ", round(time_roc$AUC[5], 2)),
    paste0("As an independence factor it has an HR of ", round(independent_prog$estimate[independent_prog$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog$conf.low[independent_prog$term == "SCORE"], 2), " - ", round(independent_prog$conf.high[independent_prog$term == "SCORE"], 2), " pval of ", independent_prog$p.value[independent_prog$term == "SCORE"], ")"),
    sep = ". "
)

cat(paste0(independent_prog$term[num_param_compare], " with its HR of ", round(independent_prog$estimate[num_param_compare], 2), " (CI 95% of ", round(independent_prog$conf.low[num_param_compare], 2), " - ", round(independent_prog$conf.high[num_param_compare], 2), " pval of ", independent_prog$p.value[num_param_compare], ")"),
    sep = ". "
)

# 9.3 Forest plot of the evaluated parameters


cox_p_metabric <- independent_prog[num_param_compare, ] %>%
  filter(
    estimate > 0.0001,
    conf.high < 100
  ) %>%
  mutate(
    term = recode_values(
      term,
      "SCORE"              ~ "Signature Score",
      "LYMPH"              ~ "Lymph Node Status",
      "AGE"                ~ "Age at Diagnosis",
      "NPI"                ~ "Nottingham Prognostic Index",
      "CHEMOYES"           ~ "Chemotherapy (Yes)",
      "HORMONEYES"         ~ "Hormone Therapy (Yes)",
      "SURGERYMASTECTOMY"  ~ "Mastectomy",
      "PAM50LumA"          ~ "Claudin subtype Luminal A",
      "PAM50LumB"          ~ "Claudin subtype Luminal B",
      "PAM50Her2"          ~ "Claudin subtype Her2-enriched",
      "PAM50claudin-low"   ~ "Claudin subtype Claudin-low",
      "PAM50Normal"        ~ "Claudin subtype Normal-like",
      "HISTMucinous"       ~ "Histology: Mucinous",
      "HISTMixed"          ~ "Histology: Mixed",
      "HISTLobular"        ~ "Histology: Lobular",
      "HISTMedullary"      ~ "Histology: Medullary",
      "INTCLUST2"          ~ "IntClust 2",
      "INTCLUST3"          ~ "IntClust 3",
      "INTCLUST5"          ~ "IntClust 5",
      "INTCLUST6"          ~ "IntClust 6",
      "INTCLUST7"          ~ "IntClust 7",
      "INTCLUST8"          ~ "IntClust 8",
      "INTCLUST9"          ~ "IntClust 9",
      "MENOPre" ~ "Premenopause",
      "RADIOYES" ~ "Radiotherapy (Yes)",
      default = term 
    ),
    term = reorder(term, estimate),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.5,
                linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "Clinical & Molecular Features", color = "Significance (p < 0.05)") +
  theme_classic(base_size = 15) +
  ggtitle("Multivariate Cox: Recurrence METABRIC") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        ) + 
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) + 
  labs(x = "Hazard Ratio (log scale)", 
       y = "", 
       color = "Significance (p < 0.05)",
       tag = "A")


