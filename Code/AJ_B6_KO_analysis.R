#Updated analysis for AJ_B6 RNA-Seq. This data includes all WT and KO. 
#Run DESeq2 on everything, specifying KOs as contrast to see if we get same results as KO only analysis

#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/")

#Load libraries
#Data manipulation
library("tidyverse")
library("dplyr")

#Differential gene expression
library("DESeq2")

#Graphics and visualizations
library("ggplot2")
library("ggrepel")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("pheatmap")
library("EnhancedVolcano")

#GSEA analysis
library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationDbi")

#Read in the raw data
raw <- read.csv("Counts/final_counts.csv", header = TRUE, sep = "\t")

#Obtain gene MGI symbols
raw$Symbols <- mapIds(org.Mm.eg.db, key = raw$Geneid, column = "SYMBOL", 
                      keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no MGI annotation or unknown genes and Rikens
raw <- raw[!is.na(raw$Symbols),]
raw <- raw[!grepl("^Gm\\d+$", raw$Symbol),]

#Format rownames
raw$Gene <- paste(raw$Geneid, raw$Symbols, sep = " - ")
rownames(raw) <- raw$Gene
raw$Gene <- NULL
raw$Geneid <- NULL
raw$Symbols <- NULL


#Create a design table
design <- data.frame(Sample = rep(c("AJ Knockout", "AJ Wildtype", "B6 Knockout", "B6 Wildtype"),
                                  c(12,12,11,12)),
                     Group = rep(c("AJ_KO", "AJ_WT", "B6_KO", "B6_WT"),
                                 c(12,12,11,12)))
rownames(design) <- colnames(raw)

#Sanity check: are all colnames in raw present in the design table?
#are the colnames in raw in the same order as rownames in design?
all(colnames(raw) %in% rownames(design))
all(colnames(raw) == rownames(design))

#Constructing a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = design,
                              design = ~ Group)

#Run DESeq2 algorithm
dds <- DESeq(dds)

#Generate a PCA plot
#Create a summarized experiment
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), colData = colData(dds))

#Plot PCA
PCA <- plotPCA(DESeqTransform(se), intgroup = "Group") +
  ggtitle("AJ vs B6")
PCA

#Get dds results for knockout comparison
KO <- results(dds, contrast = c("Group", "AJ_KO", "B6_KO"))
KO.res <- as.data.frame(KO)
KO.counts <- as.data.frame(counts(dds, normalized = TRUE))
KO_results <- merge(KO.res, KO.counts, by = 0, all = TRUE)
KO_results <- na.omit(KO_results)
KO_results$Symbols <- gsub("^[^-]+-(.*)$", "\\1", KO_results$Row.names)
KO_results <- KO_results[!grepl("^Gm\\d+$", KO_results$Symbols),]

#Filter the results based on padj and log2FoldChange
KO.res <- na.omit(KO.res)
KO.top <- KO.res[KO.res$padj < 0.05 & abs(KO.res$log2FoldChange) > 1, ]
KO.top <- KO.top[order(KO.top$log2FoldChange, decreasing = TRUE), ]

#Determine number up and down
KOup <- KO.top[KO.top$log2FoldChange > 1, ]
KOdown <- KO.top[KO.top$log2FoldChange < -1, ]


###if you're filtering out strain differences, run the following steps
#Get dds results for wildtype comparison (strain differences)
{
  WT <- results(dds, contrast = c("Group", "AJ_WT", "B6_WT"))
  WT.res <- as.data.frame(WT)
  WT.res$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(WT.res))
  
  #Filter the results based on padj and log2FoldChange
  WT.res <- na.omit(WT.res)
  WT.top <- WT.res[WT.res$padj < 0.05 & abs(WT.res$log2FoldChange) > 2, ]
  WT.top <- WT.top[order(WT.top$log2FoldChange, decreasing = TRUE), ]
  
  #Determine number up and down
  WTup <- WT.top[WT.top$log2FoldChange > 1, ]
  WTdown <- WT.top[WT.top$log2FoldChange < -1, ]
  
  #Determine what strain difference genes are present in the KO comparison
  test <- intersect(WT.top$Symbols, KO.top$Symbols)
  
  #Take only the genes that are phenotype-specific genes
  keep.KO <- KO.top[!KO.top$Symbols %in% test,] #65 genes
  keep_order <- rownames(keep.KO)
  keep.KO.list <- keep.KO$Symbols
  
  #Get these genes from the counts
  counts <- as.data.frame(counts(dds, normalized = TRUE))
  counts$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(counts))
  counts.keep.KO <- counts[counts$Symbols %in% keep.KO.list,]
  counts.keep.KO <- counts.keep.KO[keep_order,]
  counts.keep.KO[,c(13:24,36:47)] <- NULL
  counts.keep.KO$Symbols <- NULL
  
  #Write csv
  results <- merge(keep.KO, counts.keep.KO, by ="Symbols", all = TRUE)
  write.csv(results, file = "KOvsKO_phenotype_specific_genes_results.csv")
  
  #Pheatmap
  top.counts.mat <- as.matrix(counts.keep.KO)
  
  #Sample column
  sample_col <- data.frame(Sample = rep(c("AJ-KO", "B6-KO"),
                                        c(12,11)))
  row.names(sample_col) <- colnames(counts.keep.KO)
  
  #Generate heatmap
  heatmap <- pheatmap(log2(top.counts.mat + 1),
                      main = "AJ-KO vs B6-KO", 
                      scale = "row",
                      show_rownames = TRUE,
                      annotation_names_row = TRUE,
                      annotation_col = sample_col,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE,
                      fontsize_row = 5)
  heatmap
  
  gene_symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(keep.KO))
  
  #Volcano
  volcano <- EnhancedVolcano(keep.KO,
                             lab = gene_symbol,
                             boxedLabels = FALSE,
                             x = 'log2FoldChange',
                             y = 'padj',
                             legendPosition = "right",
                             legendLabels = c("NS", "L2FC", "Sig", "Sig & L2FC"),
                             title = "B6-KO vs. AJ-KO")
  volcano
}


#####Re-do analysis in a different way
#Lump together AJ KO and WT (classify as AJ) and do the same for B6
#Filter out things the same way, do GO and KEGG analysis

#Read in the raw data again
data <- read.csv("Counts/final_counts.csv", header = TRUE, sep = "\t")

#Obtain gene MGI symbols
data$Symbols <- mapIds(org.Mm.eg.db, key = data$Geneid, column = "SYMBOL", 
                      keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no MGI annotation or unknown genes and Rikens
data <- data[!is.na(data$Symbols),]
data <- data[!grepl("^Gm\\d+$", data$Symbol),]

#Format rownames
data$Gene <- paste(data$Geneid, data$Symbols, sep = " - ")
rownames(data) <- data$Gene
data$Gene <- NULL
data$Geneid <- NULL
data$Symbols <- NULL

#Create a design table
design2 <- data.frame(Sample = rep(c("AJ","B6"),
                                  c(24,23)),
                     Group = rep(c("AJ", "B6"),
                                 c(24,23)))
rownames(design2) <- colnames(data)

#Sanity check: are all colnames in raw present in the design table?
#are the colnames in raw in the same order as rownames in design?
all(colnames(data) %in% rownames(design2))
all(colnames(data) == rownames(design2))

#Constructing a DESeqDataSet object
dds2 <- DESeqDataSetFromMatrix(countData = data,
                              colData = design2,
                              design = ~ Group)

#Run DESeq2 algorithm
dds2 <- DESeq(dds2)

#Save results
res <- results(dds2, contrast = c("Group", "AJ", "B6"))
res.df <- as.data.frame(res)
res.df <- na.omit(res.df)
counts <- as.data.frame(counts(dds, normalized = TRUE))
res.df$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.df))

res_counts <- merge(res.df, counts, by = 0, all = TRUE)
res_counts <- na.omit(res_counts)

#Filter results
res.top <- res_counts[res_counts$padj < 0.05 & abs(res_counts$log2FoldChange) > 1, ]
res.top <- res.top[order(res.top$log2FoldChange, decreasing = TRUE), ]
res.top.list <- res.top$Symbols
#Get ensembl IDs
res.top$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(res.top))
#Get Entrez IDs
res.top$Entrez <- mapIds(org.Mm.eg.db, key = res.top$Ensembl,
                         column = "ENTREZID", keytype = "ENSEMBL",
                         multiVals = "first")

res.top <- res.top[,c(55,1:54)]
write.csv(res.top, file = "log2FC_1_degs.allVsall.csv")

res.top.order <- rownames(res.top)

#Save results as a csv
res.top.counts <- counts[counts$Symbols %in% res.top.list,]
res.top.counts <- res.top.counts[res.top.order,]
resultsAJB6 <- merge(res.top, res.top.counts, by = "Symbols", all = TRUE)
write.csv(resultsAJB6, file = "AJB6_AllvsAll_top_results.csv")

#Do GO and KEGG analysis
#Extract the ordered log2FoldChange column from top a.deg.top DF
deg.list <- res.top$log2FoldChange




#Name the gene list with the corresponding Entrez IDs
names(deg.list) <- res.top$Entrez

#Run GSEA for GO terms
gse <- gseGO(deg.list,
               ont = "bp",
               OrgDb = "org.Mm.eg.db",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               eps = 1e-300)

#Save GSE results as a DF
gse.res <- as.data.frame(gse)
gse.readable <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gse.readable <- as.data.frame(gse.readable)

#Write csv
write.csv(gse.readable, file = "All_vs_All_GO.csv")

#Plot
gse.dotplot <- dotplot(gse, showCategory = 22, split = ".sign", 
                       font.size = 7.5,
                       label_format = 50,
                       title = "GSEA (GO) AJ vs. B6") +
  facet_grid(.~.sign)
gse.dotplot

#KEGG
kegg <- gseKEGG(deg.list,
                  organism = "mmu",
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  eps = 1e-300)

kegg.readable <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

#Write csv
write.csv(kegg.readable, file = "All_vs_All_KEGG.csv")


#Save KEGG results as a DF
kegg.res <- as.data.frame(kegg.readable)






