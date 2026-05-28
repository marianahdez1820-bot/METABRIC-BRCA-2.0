library(tidyverse)

# If else block that determines existence of either of the GSE databases ROC data and assigns it to a common object

pot_roc.gse <-
  if(!exists("plot_roc.gse96058") & exists("plot_roc.gse2034")){
    plot_roc.gse2034
  }else if(!exists("plot_roc.gse2034") & exists("plot_roc.gse96058")){
    plot_roc.gse96058
  }else{
    print("None or both objects exist, select 1")
  }

# Same for labels

facet_labels.gse <-
  if(!exists("facet_labels.gse96058") & exists("facet_labels.gse2034")){
    facet_labels.gse2034
    GSE_label <- "GSE2034"
  }else if(!exists("facet_labels.gse2034") & exists("facet_labels.gse96058")){
    
    facet_labels.gse96058
    GSE_label <- "GSE96058"
  }else{
    print("None or both objects exist, select 1")
  }

# Object with all the ROC data from thje different databases

plot_roc_total <- bind_rows(plot_roc, pot_roc.gse, plot_roc.tcga)

# Object with all the label data from all the databases

facet_labels_final <- bind_rows(facet_labels, facet_labels.tcga, facet_labels.gse) %>%
  mutate(y_pos = case_when(
    data_set == "METABRIC" ~ 0.10,
    data_set == "GSE96058" ~ 0.18,
    data_set == "GSE2034" ~ 0.18,
    data_set == "TCGA"     ~ 0.26
  ))

 # Plot

ggplot(plot_roc_total, aes(x = FP, y = TP, color = data_set)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap( ~ Time) +
  scale_color_manual(values = c("METABRIC" = "#c388ff", "GSE2034" = "#DD6868", "TCGA" = "#B8D369")) +
  geom_text(data = facet_labels_final, 
            aes(x = 0.75, y = y_pos, label = AUC_Text, color = data_set), 
            size = 6, fontface = "bold") +
  labs( 
    title = "Time-Dependent ROC Curves",
    subtitle = "Comparing model performance across different horizons",
    x = "1 - Specificity (FP)",
    y = "Sensitivity (TP)",
    color = "Time Point") +
  theme(
    plot.title = element_text(size = 22, face = "bold"),  
    plot.subtitle = element_text(size = 16),              
    axis.title = element_text(size = 16),                 
    axis.text = element_text(size = 12),                  
    strip.text = element_text(size = 14, face = "bold"),  
    legend.title = element_text(size = 14),               
    legend.text = element_text(size = 12)                 
  )

########################################################################################

global_roc_all <- bind_rows(global_roc, global_roc.gse2034, global_roc.tcga)





ggplot(global_roc_all, aes(x = 1 - specificity, y = sensitivity, color = label)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_classic() +
  labs(title = "Global ROC curves")
