
# 6.- Validation ----------------------------------------------------------

# This is the predict phase

gse2034_results <- predict(final_fit, new_data = proof_genes_pt_gse2034, type = "linear_pred") %>%
  bind_cols(proof_genes_pt_gse2034 %>% 
              rownames_to_column("file_name"))


# 6.3 To get the p value

validation_test <- coxph(Surv(EVENT_MON,  EVENT_STAT) ~ .pred_linear_pred, data = gse2034_results)

summary(validation_test)

# 6.8 Calculate the actual Concordance Index

c_index_results.2034 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = gse2034_results)

library(survminer)

# Create risk groups based on the median of the predictions

gse2034_results <- gse2034_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))

# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = gse2034_results)

# Plot

ggsurvplot(km_fit, 
           data = gse2034_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("High risk", "Low risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 180),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3, 
           palette = c("#E41A1C", "#377EB8"),
)


gse2034_results <- gse2034_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

gse2034_split <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ ., 
  data = gse2034_results,
  cut = 70, 
  episode = "time_group",
  id = "patient_id"
)

# Run Cox again
gse2034_split$risk_group <- relevel(gse2034_split$risk_group, ref = "Low Risk")


summary_gse2034 <- summary(coxph(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group * strata(time_group), data = gse2034_split))


library(timeROC)

# Area under the curve per time

res_auc.gse2034 <- timeROC(T = gse2034_results$EVENT_MON,
                   delta = gse2034_results$ EVENT_STAT,
                   marker = -gse2034_results$pred_z,
                   cause = 1, # The event code
                   times = c(12, 36, 60, 72, 120), # 3, 5, and 10 years
                   iid = TRUE)

# 6.9.2 View the AUC values

auc_gse2034 <- res_auc.gse2034$AUC %>% 
  as.data.frame()

# View the AUC values

print(auc_gse2034)


gse2034_results <- 
  gse2034_results %>% 
  mutate(EVENT_STAT = factor(EVENT_STAT))

global_roc.gse2034 <- roc_curve(gse2034_results,
                             EVENT_STAT,
                             .pred_linear_pred
) %>%
  mutate(label = "GSE2034")


# 3.2 Combine all time points into a long dataframe

plot_roc.gse2034 <- map_df(c(12, 36, 60, 72, 120), function(i) { # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc.gse2034$FP[, time_label],
    TP = res_auc.gse2034$TP[, time_label],
    Time = factor(i),
    data_set = "GSE2034"
  )
})


# 3.3 Labels for faceted plot

facet_labels.gse2034 <- data.frame(
  Time = factor(c(12, 36, 60, 72, 120)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(res_auc.gse2034$AUC[1:5]), 3), "%"),
  data_set = "GSE2034"
)

# 3.4 Plot with facet

ggplot(plot_roc.gse2034, aes(x = FP, y = TP)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Time, ncol = 2) +
  theme_bw() +
  labs(
    title = "ROC Curves by Time Point",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  geom_text(data = facet_labels.gse2034, 
            aes(x = 0.75, y = 0.1, label = AUC_Text), 
            size = 4, fontface = "bold")




# Multivariate regression cox with clinical data

proof_genes_pt.gse2034.cox <- 
  proof_genes_pt_gse2034 %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("file_name") %>% 
  mutate(SCORE = gse2034_results$.pred_linear_pred
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                SCORE,
                EVENT_MON,
                EVENT_STAT
  ) %>% 
  na.omit()

gse2034_split.cox <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ ., 
  data = proof_genes_pt.gse2034.cox,
  cut = 70, 
  episode = "time_group",
  id = "patient_id"
)

# 2. Fit the model with the interaction

cox_model.gse2034 <- coxph(
  Surv(tstart, EVENT_MON, EVENT_STAT) ~ SCORE:strata(time_group), 
  data = gse2034_split.cox
)


cox_model.gse2034 <- coxph(surv_obj ~ SCORE:strata(time_group), 
                                  data = gse2034_split.cox) 

independent_prog.gse2034 <- 
  cox_model.gse2034 %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)


num_param_compare <- c(9:21)

cat(paste0("The signature got a C-score of ", round(c_index_results.2034$concordance, 2)),
    paste0("an HR of ", round(summary_gse2034$coefficients[2], 2), " (CI 95% of ", round(summary_gse2034$conf.int[3], 2), " - ", round(summary_gse2034$conf.int[4], 2), " pval ", summary_gse2034$coefficients[5], ")"),
    paste0("AUC at 3 years of ", round(auc_gse2034[1,], 2), " at 5 years of ", round(auc_gse2034[2,], 2), " and at 6 years of "),
    paste0("As an independence factor it has an HR of ", round(independent_prog.gse2034$estimate[independent_prog.gse2034$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.gse2034$conf.low[independent_prog.gse2034$term == "SCORE"], 2), " - ", round(independent_prog.gse2034$conf.high[independent_prog.gse2034$term == "SCORE"], 2), " pval of ", independent_prog.gse2034$p.value[independent_prog.gse2034$term == "SCORE"], ")"),
    sep = ". "
)






# 4.5 Forest plot ignoring values that tend to infinite

independent_prog.gse2034 %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  mutate(
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic() +
  scale_y_discrete(labels = c(
    "SCORE:strata(time_group)time_group=1" = "Time Group 1",
    "SCORE:strata(time_group)time_group=2" = "Time Group 2"
  ))


# 6.- Outlier analysis ----------------------------------------------------


# 6.1 Utilize the fitted object to make predictions based on time

# 6.1.1 Change to numeric 

gse2034_results <- 
  gse2034_results %>% 
  mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))

# 6.1.2 Augment on different time points

model_diagnostics <- augment(
  final_fit, 
  new_data = gse2034_results, 
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)

# 6.2 Create list to then add the results

list <- list()


#> 6.3 For loop that at each desired time point calculates the bias scores, extracts patients with high bias scores
#> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time

for (i in c(36, 60, 120)) {
  
  eval_time <- i
  
  # 6.3.1 Unnest the predictions and find the biggest outliers on a set point and event
  
  outliers <- model_diagnostics %>%
    dplyr::select(file_name, EVENT_MON, EVENT_STAT, .pred) %>%
    unnest(.pred) %>%
    filter(.eval_time == eval_time) %>% 
    arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive
  
  
  # 6.3.2 Object with metadata and score characteristics
  
  outlier_summary <- outliers %>%
    inner_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".drop")) %>%
    dplyr::select(
      file_name, 
      .pred_survival, 
      EVENT_STAT, 
      EVENT_MON, 
      BRAIN_REL, 
      .eval_time
    ) 
  
  
  # 6.3.3 Create bias score for defined time
  
  outliers_bias <- outlier_summary %>%
    mutate(
      bias_score = (((1 - EVENT_STAT) - .pred_survival) ^ 2) * ((EVENT_MON - .eval_time) * ( 1 - 2 * (EVENT_STAT)))
    ) %>%
    arrange(desc(bias_score))
  
  # 6.4 Identify the highest bias patients
  
  # 6.4.1 Obtain mean and sd and then filter baed on patients higher than determined SD
  
  extreme_outliers <- 
    outliers_bias %>% 
    mutate(mean_bias = mean(bias_score),
           sd_bias = sd(bias_score)) %>% 
    filter(bias_score > (mean_bias + 1 * sd_bias))
  
  # 6.4.2 Obtain their IDs
  
  top_bias_ids <- extreme_outliers$file_name
  
  print(length(top_bias_ids))
  
 
  # 6.5.1 Add a column identifying patients as top bias or not
  
  outliers <- 
    outliers %>% 
    mutate(quadrant = case_when(
      file_name %in% top_bias_ids & EVENT_STAT == 0 ~ 2,
      file_name %in% top_bias_ids & EVENT_STAT == 1 ~ 1,
      TRUE ~ 0
    ))
  
  print(outliers %>% 
          group_by(EVENT_STAT) %>% 
          dplyr::count(quadrant))
  
  # 6.5.2 Plot
  
  theme_embedded <- theme_classic(base_size = 25) + 
    theme(
      legend.position = c(0.95, 0.3), # Adjust coordinates (x, y) from 0 to 1
      legend.background = element_rect(fill = alpha("white", 0.5))
    )
  
  # 6.5.2.1 Plot colored by EVENT_STAT
  
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
  
  # 6.5.2.2 Plot colored by quadrant
  
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
  
  # 6.5.3 Combine and stack
  
  print( p1 + p2)
  
  list[[i]] <- top_bias_ids
  
}



# 7.- Other scores --------------------------------------------------------

# 7.1 Brier score

eval_results.gse2034 <- final_fit %>%
  augment(new_data = proof_genes_pt_gse2034, eval_time = c(36, 60, 120)) 

performance.gse2034 <- eval_results.gse2034 %>%
  brier_survival(truth = surv_obj, .pred)

print(performance.gse2034)

# 7.2 Martingale and Schofeild residuals



cox.zph(cox_model.gse2034)



ggcoxzph(cox.zph(cox_model.gse2034))



ggcoxdiagnostics(cox_model.gse2034, type = "martingale",
                 linear.predictions = FALSE, ggtheme = theme_bw())

