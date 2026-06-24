library(tidyverse)
library(tidymodels)
library(survival)
library(logistf) 


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

# 1.2.2 Create list to then add the results

outlier_time_list <- list()

no_img <- 0

#> 1.3 For loop that at each desired time point calculates the misclassification scores, extracts patients with high misclassification scores
#> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time

for (i in c(36, 60, 120)) {
  
  eval_time <- i
  
  # 1.3 Unnest the predictions and find the biggest outliers on a set point and event
  
  outliers <- model_diagnostics %>%
    dplyr::select(PATIENT_ID, EVENT_MON, EVENT_STAT, .pred, .pred_time) %>%
    unnest(.pred) %>%
    filter(.eval_time == eval_time) %>%
    arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive
  
  
  # 1.4 Object with metadata and score characteristics
  
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
  
  
  # 1.5 Create misclassification score for defined time
  
  outliers_misclassification <- outlier_summary %>%
    mutate(
      misclassification_score = (((1 - EVENT_STAT) - .pred_survival) ^ 2) * (EVENT_MON - .eval_time) * ( 1 - (2 * EVENT_STAT))
    ) %>%
    arrange(desc(misclassification_score))
  
  # 1.6 Identify the highest misclassification patients
  
  # 1.6.1 Obtain mean and sd and then filter baed on patients higher than determined SD
  
  extreme_outliers <-
    outliers_misclassification %>%
    mutate(mean_misclassification = mean(misclassification_score),
           sd_misclassification = sd(misclassification_score)) %>%
    filter(misclassification_score > (mean_misclassification + (1 * sd_misclassification)))
  
  # 1.6.2 Obtain their IDs
  
  top_misclassification_ids <- extreme_outliers$PATIENT_ID
  
  print(length(top_misclassification_ids))
  
  # 1.6.3 Identify different characteristics of patients identified as top misclassification
  
  print(outliers_misclassification %>%
          filter(PATIENT_ID %in% top_misclassification_ids) %>%
          group_by(RADIO_THERAPY) %>%
          summarise(
            count = n(),
            avg_pred_event = mean(.pred_survival),
            avg_event_time = mean(EVENT_MON)
          ) %>%
          arrange(desc(count)),
        
        n = 200
  )
  
  # 1.7.1 Add a column identifying patients as top misclassification or not
  
  outliers <-
    outliers %>%
    mutate(quadrant = case_when(
      PATIENT_ID %in% top_misclassification_ids & EVENT_STAT == 0 ~ 2,
      PATIENT_ID %in% top_misclassification_ids & EVENT_STAT == 1 ~ 1,
      TRUE ~ 0
    ))
  
  print(outliers %>%
          group_by(EVENT_STAT) %>%
          dplyr::count(quadrant))
  
  # 1.7.2 Plot
  if(no_img == 0){
    
    theme_embedded <- theme_classic(base_size = 25) +
      theme(
        legend.position = c(0.93, 0.2), # Adjust coordinates (x, y) from 0 to 1
        legend.background = element_rect(fill = alpha("white", 0.5))
      )
    
    # 1.7.2.1 Plot colored by EVENT_STAT
    
    p1 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(EVENT_STAT), shape = factor(EVENT_STAT))) +
      geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
      stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
      geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("#DB6D00FF", "#920000FF"),
                         labels = c("Alive", "Deceased")) + 
      scale_shape_manual(
        values = c(16, 17),
        labels = c("Alive", "Deceased")
      ) + 
      labs(
        title = "Event Status Distribution",
        x = paste0("Actual Event Time (Months)", eval_time),
        y = "Predicted Survival",
        color = "Event Stat",
        shape = "Event Stat"
      ) +
      theme_embedded
    
    # 1.7.2.2 Plot colored by quadrant
    
    p2 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(quadrant), shape = factor(EVENT_STAT))) +
      geom_point(size = 2, alpha = 0.7) +
      stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) +
      geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("#006DDBFF","#490092FF",  "#B66DFFFF"),
                         labels = c("Correct", "Early event", "No event")) + 
      scale_shape_manual(
        values = c(16, 17), 
        labels = c("Alive", "Deceased")
      ) + 
        labs(
        title = "Misclassification group Analysis",
        x = paste0("Actual Event Time (Months)", eval_time),
        y = "Predicted Survival",
        color = "Groups",
        shape = "Event Stat"
      ) +
      theme_embedded
    
    # 1.7.3 Combine and stack
    
    print( p1 + p2)
  }else{
    
  }
  
  outlier_time_list[[i]] <- top_misclassification_ids
  
}

# 1.8 Obtain patients found on all of the iterations of the for loop as top misclassification patients

misclassification_interesct <- intersect(intersect(outlier_time_list[[36]], outlier_time_list[[60]]) , outlier_time_list[[120]])

# 1.8.2 Similar but all unique so even if they appear once we register them

misclassification_diff_id <- unique(c(outlier_time_list[[36]], outlier_time_list[[60]], outlier_time_list[[120]]))



# 2.- Identifying cause of misclassification ------------------------------


# 2.1 Mutate characteristics to evaluate in the further analysis

# 2.0 Build a clean baseline dataset entirely independent of loop artifacts
outlier_cause <- id %>% 
  inner_join(ml_metadata, by = "PATIENT_ID", suffix = c("", ".drop")) %>%
  mutate(
    # Binarize IntClust to low and high grade
    intcluster = ifelse(INTCLUST %in% c("3", "4ER+", "7", "8"), "Low", "High"), 
    
    # Assign quadrants based on your combined list of any-timeframe outliers
    quadrant = case_when(
      PATIENT_ID %in% misclassification_diff_id & EVENT_STAT == 1 ~ 1,
      PATIENT_ID %in% misclassification_diff_id & EVENT_STAT == 0 ~ 2,
      TRUE ~ 3
    ),
    
    # Binary targets for your downstream regressions
    unexpected_event = ifelse(quadrant == 1, 1, 0), 
    exceptional      = ifelse(quadrant == 2, 1, 0),
    NPI_bin = ifelse(NPI > (mean(NPI) + 0.11),
                     1,
                     0),
    age_bin = ifelse(AGE_AT_DIAGNOSIS > (mean(AGE_AT_DIAGNOSIS) + 1),
                     1,
                     0),
    treatment = case_when(
      HORMONE_THERAPY == "YES" & RADIO_THERAPY == "NO" & CHEMOTHERAPY == "NO" ~ "Hormone",
      HORMONE_THERAPY == "NO" & (xor(RADIO_THERAPY == "YES", CHEMOTHERAPY == "YES")) ~ "Chemo|Radio",
      HORMONE_THERAPY == "NO" & (RADIO_THERAPY == "YES" & CHEMOTHERAPY == "YES") ~ "Chemo/Radio",
      HORMONE_THERAPY == "YES" & (xor(RADIO_THERAPY == "YES", CHEMOTHERAPY == "YES")) ~ "Hormone-Chemo/Radio",
      HORMONE_THERAPY == "NO" & RADIO_THERAPY == "NO" & CHEMOTHERAPY == "NO" ~ "None",
      HORMONE_THERAPY == "YES" & RADIO_THERAPY == "YES" & CHEMOTHERAPY == "YES" ~ "ALL"
    ),
    treatment_atall = ifelse(HORMONE_THERAPY == "NO" & RADIO_THERAPY == "NO" & CHEMOTHERAPY == "NO",
                             "None",
                             "Some"
)
  )

# 2.1 Comparison between groups and p val adjustment

for (a in c(1, 2)) {
  
  
  # 2.2 Subset data depending on cuadrant analysed
  
  df_misclassification <- subset(outlier_cause, quadrant %in% c(a, 3)) # First Q1 vs Q3 (unexpected vs correct) thenb Q2 vs Q3 (exceptional vs correct)
  
  p_multiple_eval_misclass <- c(
    intclust = chisq.test(table(df_misclassification$quadrant, df_misclassification$INTCLUST), simulate.p.value = TRUE)$p.value, # Many variables unfit for fishers but because of possible sparsity we p simulate
    histology= chisq.test(table(df_misclassification$quadrant, df_misclassification$HISTOLOGICAL_SUBTYPE), simulate.p.value = TRUE)$p.value, # Same
    intcluster = fisher.test(table(df_misclassification$quadrant, df_misclassification$HORMONE_THERAPY))$p.value,
    hormone  = fisher.test(table(df_misclassification$quadrant, df_misclassification$intcluster))$p.value, # Binary withprobaiblity of <5 variables in cells
    chemo    = fisher.test(table(df_misclassification$quadrant, df_misclassification$CHEMOTHERAPY))$p.value,
    radio    = fisher.test(table(df_misclassification$quadrant, df_misclassification$RADIO_THERAPY))$p.value,
    claudin  = chisq.test(table(df_misclassification$quadrant, df_misclassification$CLAUDIN_SUBTYPE), simulate.p.value = TRUE)$p.value,
    nodes    = kruskal.test(LYMPH_NODES_EXAMINED_POSITIVE ~ quadrant, data = df_misclassification)$p.value, # Continuous non parametric data
    age      = kruskal.test(AGE_AT_DIAGNOSIS ~ quadrant, data = df_misclassification)$p.value,
    NPI      = kruskal.test(NPI ~ quadrant, data = df_misclassification)$p.value
  )
  
  # 2.3 Adjust p value with bonferroni holm and print which groups are being comopared
  
  adj_p_misclass_event <- p.adjust(p_multiple_eval_misclass, method = "holm")
  print(paste(" Q", a, " vs Q3 (Standard) Holm Adjusted P-values"))
  options(scipen = 999)
  print(round(adj_p_misclass_event, 5))
  
  
  for (e in c("CLAUDIN_SUBTYPE", "INTCLUST")) {
    
    
    # 2.4 Table with characteristic to evlauate its residuals
    
    misclass_table <- table(df_misclassification$quadrant, df_misclassification[[e]])
    misclass_table <- misclass_table[drop = TRUE]
    
    # 2.5 Run chisq with simulated p values
    
    chisq_obj_exp <- chisq.test(misclass_table, simulate.p.value = TRUE)
    
    # 2.6 Obtain residuals
    
    print(chisq_obj_exp$residuals)
  }
  
}


# 3.- Regressions ---------------------------------------------------------

outlier_cause <- outlier_cause %>% 
  mutate(
    across(CLAUDIN_SUBTYPE, ~ fct_relevel(.x, "LumB")))

# 3.1 Firth penalized logistic regression on unexpectedevent group

clean_unexpected_model_firth <- logistf(
  unexpected_event ~ HORMONE_THERAPY + CLAUDIN_SUBTYPE + CHEMOTHERAPY + RADIO_THERAPY + AGE_AT_DIAGNOSIS + NPI ,
  data = outlier_cause
)

summary(clean_unexpected_model_firth)

# 3.1.2 OR

exp(coef(clean_unexpected_model_firth))

# 3.1.3 CI OR

exp(confint(clean_unexpected_model_firth))


# 3.2 MFirth penalized regression on exceptional group


clean_exceptional_model_firth <- logistf(
  exceptional ~ HORMONE_THERAPY + intcluster + CHEMOTHERAPY + RADIO_THERAPY + AGE_AT_DIAGNOSIS + NPI ,
  data = outlier_cause
)

summary(clean_exceptional_model_firth)

# 3.2.2 OR

exp(coef(clean_exceptional_model_firth))

# 3.2.3 OR with CI

exp(confint(clean_exceptional_model_firth))
