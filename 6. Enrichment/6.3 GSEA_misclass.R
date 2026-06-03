library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(enrichplot)


res$entrez <- mapIds(org.Hs.eg.db,
                     keys=rownames(res), # A donde va a buscar
                     column="ENTREZID", # Nueva columna con ese formato
                     keytype="SYMBOL", # Que formato va a buscar en keys
                     multiVals="first") # Que hacer si hay varios del mismo Key
res$name <-   mapIds(org.Hs.eg.db,
                     keys=rownames(res), 
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")





gene_symbols <- rownames(count_data)


fit$genes <- data.frame(Gene.symbol = gene_symbols)

toptable <- topTable(fit, n = Inf)

EnhancedVolcano(toptable,
                lab = toptable$Gene.symbol,
                x = 'logFC',
                y = 'P.Value')




res_lfc <- res$logFC

# 9.1 Asignar al objeto de Log Fold los nombres de ENSEMBL

names(res_lfc) <- rownames(res)

# 9.2 Eliminar los NA

gene_list <- na.omit(res_lfc)


# 9.3 Ordena en orden descendente

gene_list <- sort(gene_list, decreasing = TRUE)


# 9.4 Objeto GSEA de Gene ontology

gse <- gseGO(
  geneList = gene_list, # Objeto creado anteriormente que ordena de manera descendente
  ont = "BP", # Subontologias (BP = Biological process, CC = Celular component, MF = Molecular Function)
  keyType = "SYMBOL", # El codigo de identificación de los genes
  pvalueCutoff = 0.05, 
  OrgDb = "org.Hs.eg.db", # Organismo en este caso Homo sapiens (a diferencia de org.Mm.eg.db que es de raton o org.Sc.eg.db que es de levadura) 
  minGSSize = 10, # Tamaño minimo de ada set para analisis
  maxGSSize = 100 # Tamaño maximo de genes anotados para analizar
)


# 9.4.1 Observar como data frame

gse_df <- as.data.frame(gse)

View(gse_df)


