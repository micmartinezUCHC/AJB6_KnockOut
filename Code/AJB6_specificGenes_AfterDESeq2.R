#The purpose of this script is to analyze AJ/B6 KO and AJ/B6 WT to look for genes that could drive differences
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

#Read in the KO counts
KOraw <-  read.csv("AJ_B6_KO_counts.csv", header = TRUE, sep = ",")

#Obtain gene MGI symbols
KOraw$Symbols <- mapIds(org.Mm.eg.db, key = KOraw$Geneid, column = "SYMBOL", 
                      keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no MGI annotation or unknown genes and Rikens
KOraw <- KOraw[!is.na(KOraw$Symbols),]
KOraw <- KOraw[!grepl("^Gm\\d+$", KOraw$Symbol),]
KOraw <- KOraw[!grepl("Rik$", KOraw$Symbol),] #Remove things like this: 9930111J21Rik2
KOraw <- KOraw[!grepl("AA\\d+$", KOraw$Symbol),]

#Format rownames
KOraw$Gene <- paste(KOraw$Geneid, KOraw$Symbols, sep = " - ")
rownames(KOraw) <- KOraw$Gene
KOraw$Gene <- NULL
KOraw$Geneid <- NULL
KOraw$Symbols <- NULL

#Create a design table
KOdesign <- data.frame(Sample = rep(c("AJ Knockout", "B6 Knockout"),
                                  c(12,11)),
                     Group = rep(c("AJ_KO","B6_KO"),
                                 c(12,11)))
rownames(KOdesign) <- colnames(KOraw)

#Sanity check: are all colnames in raw present in the design table?
#are the colnames in raw in the same order as rownames in design?
all(colnames(KOraw) %in% rownames(KOdesign))
all(colnames(KOraw) == rownames(KOdesign))

#Create DESeq2 object
KOdds <- DESeqDataSetFromMatrix(countData = KOraw,
                              colData = KOdesign,
                              design = ~ Group)

#Run DESeq2 algorithm
KOdds <- DESeq(KOdds)

#View results
KO <- as.data.frame(results(KOdds))
KO$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(KO))
KO_counts <- as.data.frame(counts(KOdds, normalized = TRUE))
KO_counts$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(KO_counts))
KO.res <- merge(KO, KO_counts, by = "Symbol", all = TRUE)

#Filter the results
KO.res <- na.omit(KO.res)
KO.res.filt <- KO.res[KO.res$padj < 0.05 & abs(KO.res$log2FoldChange) > 1,]


#Write UNFILTERED results as a csv file
write.csv(KO.res.filt, file = "filtered_KO_DEGs.csv")
write.csv(KO.res, file = "unfiltered_KO_DEGs.csv")

#Now we look for the AJ and B6 specific genes
#Find the genes that have low read counts in AJ
minimumCountpergene <- 1
MinSampleWithminimumgeneCounts <- 10

#Low read counts in AJ counts
B6KO.spec <- KO.res.filt[rowSums(data.frame(KO.res.filt[,8:19] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]
AJKO.spec <- KO.res.filt[rowSums(data.frame(KO.res.filt[,20:ncol(KO.res)] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#Write these tables as csv files
write.csv(B6KO.spec, file = "B6KO.specificGenes_filtered.csv")
write.csv(AJKO.spec, file = "AJKO.specificGenes_filtered.csv")

#Read in the heatmap input file
KO_heatmap <- read.csv("KO_heatmap_input.csv", header = TRUE, sep = ",")
rownames(KO_heatmap) <- KO_heatmap$Symbol
KO_heatmap$Symbol <- NULL
KO_heatmap.mat <- as.matrix(KO_heatmap)

#Sample column
sample_col <- data.frame(Sample = rep(c("AJ-KO", "B6-KO"),
                                      c(12,11)))
row.names(sample_col) <- colnames(KO_heatmap)

KOheatmap <- pheatmap(log2(KO_heatmap.mat + 1),
                    main = "AJ-KO vs B6-KO", 
                    scale = "row",
                    show_rownames = TRUE,
                    annotation_names_row = TRUE,
                    annotation_col = sample_col,
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    fontsize_row = 5)
KOheatmap


#Read in the wildtype counts
WTraw <-  read.csv("AJ_B6_WT_counts.csv", header = TRUE, sep = ",")

#Obtain gene MGI symbols
WTraw$Symbols <- mapIds(org.Mm.eg.db, key = WTraw$Geneid, column = "SYMBOL", 
                        keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no MGI annotation or unknown genes and Rikens
WTraw <- WTraw[!is.na(WTraw$Symbols),]
WTraw <- WTraw[!grepl("^Gm\\d+$", WTraw$Symbol),]
WTraw <- WTraw[!grepl("Rik$", WTraw$Symbol),] #Remove things like this: 9930111J21Rik2

#Format rownames
WTraw$Gene <- paste(WTraw$Geneid, WTraw$Symbols, sep = " - ")
rownames(WTraw) <- WTraw$Gene
WTraw$Gene <- NULL
WTraw$Geneid <- NULL
WTraw$Symbols <- NULL

#Create a design table
WTdesign <- data.frame(Sample = rep(c("AJ Wildtype", "B6 Wildtype"),
                                    c(12,12)),
                       Group = rep(c("AJ_WT","B6_WT"),
                                   c(12,12)))
rownames(WTdesign) <- colnames(WTraw)

#Sanity check: are all colnames in raw present in the design table?
#are the colnames in raw in the same order as rownames in design?
all(colnames(WTraw) %in% rownames(WTdesign))
all(colnames(WTraw) == rownames(WTdesign))

#Create DESeq2 object
WTdds <- DESeqDataSetFromMatrix(countData = WTraw,
                                colData = WTdesign,
                                design = ~ Group)

#Run DESeq2 algorithm
WTdds <- DESeq(WTdds)

#View results
WT <- as.data.frame(results(WTdds))
WT$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(WT))
WT_counts <- as.data.frame(counts(WTdds, normalized = TRUE))
WT_counts$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(WT_counts))
WT.res <- merge(WT, WT_counts, by = "Symbol", all = TRUE)

#Filter the results
WT.res <- na.omit(WT.res)
WT.res.filt <- WT.res[WT.res$padj < 0.05 & abs(WT.res$log2FoldChange) > 1,]

#Now we look for the AJ and B6 specific genes
#Find the genes that have low read counts in AJ
minimumCountpergene <- 1
MinSampleWithminimumgeneCounts <- 10

#Low read counts in AJ counts
B6WT.spec <- WT.res.filt[rowSums(data.frame(WT.res.filt[,8:19] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]
AJWT.spec <- WT.res.filt[rowSums(data.frame(WT.res.filt[,20:ncol(WT.res)] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#Write these tables as csv files
write.csv(B6WT.spec, file = "B6WT.specificGenes_filtered.csv")
write.csv(AJWT.spec, file = "AJWT.specificGenes_filtered.csv")

#Read in the heatmap input file
WT_heatmap <- read.csv("WT_heatmap_input.csv", header = TRUE, sep = ",")
rownames(WT_heatmap) <- WT_heatmap$Symbol
WT_heatmap$Symbol <- NULL
WT_heatmap.mat <- as.matrix(WT_heatmap)

#Sample column
sample_col <- data.frame(Sample = rep(c("AJ-WT", "B6-WT"),
                                      c(12,12)))
row.names(sample_col) <- colnames(WT_heatmap)

WTheatmap <- pheatmap(log2(WT_heatmap.mat + 1),
                      main = "AJ-WT vs B6-WT", 
                      scale = "row",
                      show_rownames = TRUE,
                      annotation_names_row = TRUE,
                      annotation_col = sample_col,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE,
                      fontsize_row = 5)
WTheatmap


#Prepare for GSEA analysis
KO <- na.omit(KO)
KO$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(KO))
KO$Entrez <- mapIds(org.Mm.eg.db, key = KO$Ensembl,
                         column = "ENTREZID", keytype = "ENSEMBL",
                         multiVals = "first")
KO <- KO[order(KO$log2FoldChange, decreasing = TRUE),]

KO.degs <- KO$log2FoldChange
names(KO.degs) <- KO$Entrez

gse <- gseGO(KO.degs,
             ont = "bp",
             OrgDb = "org.Mm.eg.db",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             eps = 1e-300)

gse.res <- as.data.frame(gse)
gse.readable <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gse.readable <- as.data.frame(gse.readable)

gse.dotplot <- dotplot(gse, showCategory = 22, split = ".sign", 
                       font.size = 7.5,
                       label_format = 50,
                       title = "GSEA (GO) AJ vs. B6") +
  facet_grid(.~.sign)
gse.dotplot

kegg <- gseKEGG(KO.degs,
                organism = "mmu",
                keyType = "kegg",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 1e-300)

kegg.readable <- as.data.frame(setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))

kegg.dotplot <- dotplot(kegg, showCategory = 15, split = ".sign", 
                       font.size = 7.5,
                       label_format = 50,
                       title = "GSEA (KEGG) AJ-KO vs. B6-KO") +
  facet_grid(.~.sign)
kegg.dotplot




