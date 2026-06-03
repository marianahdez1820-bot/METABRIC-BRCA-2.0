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
           
           xlim = c(0, 220),         
           break.time.by = 50,      
           ggtheme = theme_classic(), 
           
           linewidth = 3,             
           palette = c("#c380d3", "#ff89d4"),
)


p2 <- ggsurvplot(km_fit.gse2034, 
                 data = gse2034_results, 
                 pval = FALSE,              
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
           
           xlim = c(0, 170),        
           break.time.by = 50,       
           ggtheme = theme_classic(), 
           
           linewidth = 3,                 
           palette = c("#E7B800", "#2E9FDF"), 
           
)




# list

plot_list <- list(p1, p2, p3)

# add the tags directly to the curve portion of each object

p1$plot <- p1$plot + labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 18))
p2$plot <- p2$plot + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 18))
p3$plot <- p3$plot + labs(tag = "C") + theme(plot.tag = element_text(face = "bold", size = 18))

# tagged objects

plot_list <- list(p1, p2, p3)

#  Arrange and print 
res <- arrange_ggsurvplots(plot_list, 
                           print = TRUE, 
                           ncol = 3, 
                           nrow = 1, 
                           risk.table.height = 0.22)
