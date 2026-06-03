p1 <-  ggsurvplot(fit_km,
           data = test_data,
           risk.table = TRUE,
           risk.table.fontsize = 3,
           
           title = "Recurrence METABRIC",
           ylab = "Recurrence probability",
           font.title = 22,
           legend = "bottom",
           font.legend = 14,
           legend.title = "Risk group",
           font.legend.title = 12,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 12,
           xlab = "Time (months)",
           
           xlim = c(0, 220),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_classic(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#c380d3", "#ff89d4"),
)


p2 <- ggsurvplot(km_fit.gse2034, 
                 data = gse2034_results, 
                 pval = FALSE,              # Changed to FALSE
                 risk.table = TRUE,
                 risk.table.fontsize = 3,
                 
                 title = "Validation in GSE2034",
                 font.title = 22,
                 legend = "bottom",
                 font.legend = 14,
                 legend.title = "Risk group",
                 font.legend.title = 14,
                 legend.labs =  c("High risk", "Low risk"),
                 font.legend.labs = 14,
                 xlab = "Time (months)",
                 ylim = c(0, 1),
                 xlim = c(0, 170),         
                 break.time.by = 50,      
                 ggtheme = theme_classic(), 
                 linewidth = 3,                 
                 palette = c("#d23A1C", "#377EB8")
)


p3 <- ggsurvplot(km_fit.tcga, 
           data = tcga_results, 
           risk.table = TRUE,
           risk.table.fontsize = 3,
           
           
           title = "Validation in TCGA",
           font.title = 22,
           legend = "bottom",
           font.legend = 14,
           legend.title = "Risk group",
           font.legend.title = 14,
           legend.labs = c("Low risk", "High risk") ,
           font.legend.labs = 14,
           xlab = "Time (months)",
           
           xlim = c(0, 170),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_classic(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#E7B800", "#2E9FDF"), # Custom color palette
           
)



# Plots into a list

plot_list <- list(p1, p2, p3)

# Arrange them in a grid 
res <- arrange_ggsurvplots(plot_list, 
                           print = TRUE, 
                           ncol = 3, 
                           nrow = 1, 
                           risk.table.height = 0.22)

