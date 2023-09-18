setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/")


#Load libraries
suppressWarnings({
  #Data manipulation
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(plyr))
  
  
  #Differential gene expression
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(ashr))
  suppressPackageStartupMessages(library(RNAseqQC))
  suppressPackageStartupMessages(library(vsn))
  
  #Graphics and visualizations
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggh4x))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(patchwork))
  
  
  #GSEA analysis
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(enrichplot))
})

raw <- read.csv("/users/mikemartinez/Desktop/AJB6_KnockOut/Counts/final_counts.tsv", header = TRUE, sep = "\t")
setwd("/Users/mikemartinez/Desktop/AJB6_Knockout/")

#Obtain gene symbols
raw$Symbol <- mapIds(org.Mm.eg.db, key = raw$Gene, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no symbol annotation or unknown genes
raw <- raw[!is.na(raw$Symbol),]

#Omit any gene with no MGI annotation or unknown genes and Rikens
raw <- raw[!grepl("^Gm\\d+$", raw$Symbol),]
raw <- raw[!grepl("Rik$", raw$Symbol),]
raw <- raw[!grepl("Rik\\d+$", raw$Symbol),]
raw <- raw[!grepl("^LOC", raw$Symbol),]
raw <- raw[!grepl("AK\\d+$", raw$Symbol),]
raw <- raw[!grepl("AY\\d+$", raw$Symbol),]

#Format the rownames
raw$geneIDs <- paste(raw$Geneid, raw$Symbol, sep = " - ")
rownames(raw) <- raw$geneIDs
raw$Geneid<- NULL
raw$Symbol <- NULL
raw$geneIDs <- NULL

#Assess the library sizes of each sample
df.m <- reshape2::melt(raw, id.vars =NULL)
colnames(df.m) <- c("Sample", "Counts")

#Plot library sizes violin plot
QC <- ggplot(df.m, aes(Sample,log10(Counts),fill=Sample)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
QC

#Create a design table for DESeq2
design <- data.frame(Sample = rep(c("A/J:KO", "A/J:WT", "B6:KO", "B6:WT"),
                                  c(12,12,11,12)),
                     Group = rep(c("A/J", "A/J", "B6", "B6"),
                                 c(12,12,11,12)),
                     Genotype = rep(c("A/J:KO", "A/J:WT", "B6:KO", "B6:WT"),
                                    c(12,12,11,12)))
rownames(design) <- colnames(raw)

#Check that everything in the design table is in the counts matrix
all(colnames(raw) %in% rownames(design))
all(colnames(raw) == rownames(design))

suppressWarnings({
  dds <- DESeqDataSetFromMatrix(countData = raw,
                                colData = design,
                                design = ~ Group)
})

#Note: whatever you set as the reference group is the denominator according to DESeq2 documentation
dds$Group <- relevel(dds$Group, ref = "A/J")

#QC
total_counts_QC <- plot_total_counts(dds)
total_counts_QC + ggtitle("Total Sample Counts")
library_complexity <- plot_library_complexity(dds)
library_complexity + ggtitle("Library Complexity")
variance_stable <- vst(dds)
MSD_plot <- mean_sd_plot(variance_stable)
MSD_plot + ggtitle("Mean and SD plot")

#Run DESeq2 algorithm
dds <- DESeq(dds)

#Create a summarized experiment
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), colData = colData(dds))

#Plot PCA
PCA <- plotPCA(DESeqTransform(se), intgroup = "Genotype") +
  geom_text_repel(aes(label = rownames(design)), size = 3, max.overlaps = Inf) +
  ggtitle("A/J vs B6") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("Overall_AJB6_all4Groups_PCA.pdf", PCA, width = 12, height = 8)


#Hierarchical clustering
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

#Get correlation values and plot as heatmap
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)

#View the results
res <- results(dds)
res.df <- as.data.frame(res)
res.df <- na.omit(res.df)
res.df$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(res.df))
res.df$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.df))

#Get the normalixed counts
counts <- counts(dds, normalized = TRUE)

#Order and merge with counts
res.ordered <- res.df[order(res.df$log2FoldChange, decreasing = TRUE),]
res.ordered <- merge(res.df, counts, by = 0, all = TRUE)
res.ordered <- na.omit(res.ordered)
write.csv(res.ordered, file = "AJKO_vs_B6KO_All_DEGs_NoUnknowns.csv")

#Number of genes up and down regulated
filtered <- res.ordered[res.ordered$padj < 0.05 & abs(res.ordered$log2FoldChange) > 2, ]
up <- filtered[filtered$log2FoldChange > 0, ]
dn <- filtered[filtered$log2FoldChange < 0, ]

#Make venn-diagram
library("ggvenn")

upGenes <- up$log2FoldChange
dnGenes <- dn$log2FoldChange

#Generate venn list
x <- list("Upregulated" = upGenes,
          "Downregulated" = dnGenes)
ggvenn(x)


#Get ENSEMBL IDs as their own column to map ENTREZ IDs
res.ordered$Entrez <- mapIds(org.Mm.eg.db, key = res.ordered$Ensembl,
                             column = "ENTREZID", keytype = "ENSEMBL",
                             multiVals = "first")

#Get a list of the genes
res.ordered <- res.ordered[order(res.ordered$log2FoldChange, decreasing = TRUE),]
res.ordered.genes <- res.ordered$log2FoldChange

#Assign Entrez IDs as names for the genes
names(res.ordered.genes) <- res.ordered$Entrez

#Remove duplicated Entrez IDs and their corresponding values
unique_entrez_genes <- names(res.ordered.genes[!duplicated(names(res.ordered.genes))])
unique_genes <- res.ordered.genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)

  #Run GSEA for GO terms
  #For GSEA: USE ALL THE DEGS BEFORE FILTERING
  gse <- gseGO(unique_genes,
               ont = "all",
               OrgDb = "org.Mm.eg.db",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               eps = 1e-300,
               verbose = TRUE,
               by = "fgsea")
  gse.readable <- as.data.frame(setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
  
  #Set labels
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = "Enriched in B6", suppressed = "Enriched in AJ")
  )
  
  # Create the dotplot with custom facet labels
  GO.dotplot <- dotplot(gse, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 8, 
                        label_format = 100,
                        title = "AJKO vs B6KO GSEA (GO)") + 
    facet_nested(.~.sign + ONTOLOGY~., labeller = custom_labels, scales = "free")
  
  
  #CODE TO MAKE A CUSTOM BAR PLOT FOR GSEA (GO)
  #Group the results based on BP ontology and separate into positive and negative, arrange by decreasing NES
  BP_terms <- gse.readable[gse.readable$ONTOLOGY == "BP", ]
  posBP_terms <- BP_terms[BP_terms$NES > 0,]
  posBP_terms.orderedNES <- posBP_terms[order(posBP_terms$NES, decreasing = TRUE),]
  
  negBP_terms <- BP_terms[BP_terms$NES < 0,]
  negBP_terms.orderedNES <- negBP_terms[order(negBP_terms$NES, decreasing = FALSE),]
  topnegBP_terms.orderedNES <- negBP_terms.orderedNES[1:15,]
  
  #Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
  BP <- rbind(posBP_terms.orderedNES, topnegBP_terms.orderedNES)
  BP$Sign <- ifelse(BP$NES > 0, "Enriched in B6KO", "Enriched in AJKO")
  
  #Set the GO description to a factor so it is ordrede the same way in the plot
  BP$Description <- factor(BP$Description, levels = BP$Description)
  
  #Group the results based on MF ontology and separate into positive and negative, arrange by decreasing NES
  MF_terms <- gse.readable[gse.readable$ONTOLOGY == "MF", ]
  posMF_terms <- MF_terms[MF_terms$NES > 0, ]
  posMF_terms.orderedNES <- posMF_terms[order(posMF_terms$NES, decreasing = TRUE),]
  
  negMF_terms <- MF_terms[MF_terms$NES < 0,]
  negMF_terms.orderedNES <- negMF_terms[order(negMF_terms$NES, decreasing = FALSE),]
  topnegMF_terms.orderedNES <- negMF_terms.orderedNES[1:15,]
  
  #Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
  MF <- rbind(posMF_terms.orderedNES, topnegMF_terms.orderedNES)
  MF$Sign <- ifelse(MF$NES > 0, "Enriched in B6KO", "Enriched in AJKO")
  
  #Set the GO description to a factor so it is ordrede the same way in the plot
  MF$Description <- factor(MF$Description, levels = MF$Description)
  
  #Group the results based on CC ontology and separate into positive and negative, arrange by decreasing NES
  CC_terms <- gse.readable[gse.readable$ONTOLOGY == "CC", ]
  posCC_terms <- CC_terms[CC_terms$NES > 0, ]
  posCC_terms.orderedNES <- posCC_terms[order(posCC_terms$NES, decreasing = TRUE),]
  
  negCC_terms <- CC_terms[CC_terms$NES < 0,]
  negCC_terms.orderedNES <- negCC_terms[order(negCC_terms$NES, decreasing = FALSE),]
  topnegCC_terms.orderedNES <- negCC_terms.orderedNES[1:2,]
  
  #Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
  CC <- rbind(posCC_terms.orderedNES, topnegCC_terms.orderedNES)
  CC$Sign <- ifelse(CC$NES > 0, "Enriched in B6KO", "Enriched in AJKO")
  
  #Set the GO description to a factor so it is ordrede the same way in the plot
  CC$Description <- factor(CC$Description, levels = CC$Description)
  
  #Rbind BP and MF
  BPMF <- rbind(BP, MF)
  
  #Rbind BPMF and CC
  all <- rbind(BPMF, CC)
  
  #Order all GO terms by decreasing NES and set description as factor
  all <- all[order(all$NES, decreasing = TRUE),]
  all$Description <- factor(all$Description, levels = all$Description)
  
  Pval_gradientPlot <- ggplot(all, aes(x = NES, y = Description, fill = p.adjust, label = setSize)) +
    geom_bar(stat = "identity") +
    geom_text(size = 2, color = "black") +
    facet_grid(ONTOLOGY ~ Sign, scales = "free") +
    scale_fill_gradient(low = "red", high = "blue") +
    scale_x_continuous(breaks = seq(-2,2,1)) +
    theme(axis.text.y = element_text(size = 6)) +
    ggtitle("AJKO vs B6KO: GSEA (GO)") +
    theme_bw()
  Pval_gradientPlot
  
  OntologyPlot <- ggplot(all, aes(x = NES, y = Description, fill = ONTOLOGY, label = setSize)) +
    geom_bar(stat = "identity") +
    geom_text(size = 3, color = "black") +
    facet_wrap(~Sign, scales = "free_x") +
    scale_x_continuous(breaks = seq(-2,2,1)) +
    theme(axis.text.y = element_text(size = 9)) +
    theme_bw() +
    ggtitle("AJKO vs B6KO: GSEA (GO)")
  OntologyPlot
  
  
  #Prepare significant genes
  GOR <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
  core.genes <- str_split(as.data.frame(gse)[,"core_enrichment"], "/")
  my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.5])
  
  #For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
  filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                                paste, collapse = "/")
  gse@result$core_enrichment <- filtered.core.genes
  gse.filtered <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
  
  
  GOCNET <- cnetplot(gse.filtered, node_label = "gene", foldChange = unique_genes, cex_label_gene = 0.6, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 10)
  
  
  #FOR B6#
  #Take the subset of the gse results that have + NES score
  gse.df <- as.data.frame(gse)
  gse.pos <- gse.df[gse.df$NES > 0,]
  gse@result <- gse.pos
  
  #Filter for core genes
  core.genes <- str_split(as.data.frame(gse)[,"core_enrichment"], "/")
  my.selected.genesB6 <- names(unique_genes[abs(unique_genes) > 2.0])
  filtered.core.genesB6 <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genesB6]),
                                paste, collapse = "/")

  #Set the GSEA results to be the filtered B6 results
  gse@result$core_enrichment <- filtered.core.genesB6
  gse.filteredB6 <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
  test.df <- as.data.frame(gse.filteredB6)
 
  #Plot the B6 enrichment results
  GOCNETB6 <- cnetplot(gse.filteredB6, node_label = "gene", foldChange = unique_genes, cex_label_gene = 0.6, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 10)
  
  
  #Plot enrichment plot for AJKO genes
  ES_plot.GO <- enrichplot::gseaplot2(gse, geneSetID = 1:5, title = "GO Results")
  combined_plot.GO <- wrap_plots(ES_plot.GO, nrow = 3, ncol = 1)
  
  
  #Run KEGG enrichment
  kegg <- gseKEGG(unique_genes,
                  organism = "mmu",
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  eps = 1e-300,
                  verbose = TRUE,
                  by = "fgsea")
  kegg.readable <- as.data.frame(setReadable(kegg, "org.Mm.eg.db", "ENTREZID"))
  
  core.genes <- str_split(as.data.frame(kegg)[,"core_enrichment"], "/")
  filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                                paste, collapse = "/")
  kegg@result$core_enrichment <- filtered.core.genes
  kegg.filtered <- setReadable(kegg, "org.Mm.eg.db", "ENTREZID")
  
  KEGGCNET <- cnetplot(kegg.filtered, node_label = "gene", foldChange = unique_genes, cex_label_category = 0.5, cex_label_gene = 0.5, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 4)
  
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = "up in B6KO", suppressed = "up in AJKO")
  )
  
  # Create the dotplot with custom facet labels
  KEGG.dotplot <- dotplot(kegg, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 7.5, 
                        label_format = 50,
                        title = "KEGG Enrichment") + 
    facet_grid(.~.sign, labeller = custom_labels)
  
  ES_plot.KEGG <- enrichplot::gseaplot2(kegg, geneSetID = 1:4, title = "KEGG Results")
  combined_plot.KEGG <- wrap_plots(ES_plot.KEGG, nrow = 3, ncol = 1)
  
  
#Get the AJ and B6 specific genes (from file)
AJspecific <- read.csv("AJKO.specificGenes_filtered.csv", header = TRUE, sep = ",")
AJspecificSymbols <- AJspecific$Symbol

#Filter out the AJ-specific genes from the results dataframe
AJspecGenes <- res.ordered[res.ordered$Symbol %in% AJspecificSymbols, ]

#Get the entezID associated with the AJ-specific genes
AJspecEntrez <- AJspecGenes$Entrez

#Perform ORA
oraGO <- enrichGO(gene = AJspecEntrez,
                  universe = names(unique_genes),
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.20,
                  readable = TRUE)
  
oraGOres <- as.data.frame(oraGO)

oraGO.dotplot <- dotplot(oraGO, 
                      showCategory = 30, 
                      font.size = 7.5, 
                      label_format = 50,
                      title = "ORA: AJKO specific Genes") 


B6specific <- read.csv("B6KO.specificGenes_filtered.csv", header = TRUE, sep = ",")
B6specificSymbols <- B6specific$Symbol

#Filter out the AJ-specific genes from the results dataframe
B6specGenes <- res.ordered[res.ordered$Symbol %in% B6specificSymbols, ]

#Get the entezID associated with the AJ-specific genes
B6specEntrez <- B6specGenes$Entrez

#Perform ORA
oraGOb6 <- enrichGO(gene = B6specEntrez,
                  universe = names(unique_genes),
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.20,
                  readable = TRUE)

oraGOb6res <- as.data.frame(oraGOb6)

oraGOb6.dotplot <- dotplot(oraGOb6, 
                         showCategory = 30, 
                         font.size = 7.5, 
                         label_format = 50,
                         title = "ORA: B6KO specific Genes") 


oraAJkegg <- enrichKEGG(gene = AJspecEntrez,
                      universe = names(unique_genes),
                      organism = "mmu",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.20)
oraAJkegg.res <- as.data.frame(setReadable(oraAJkegg, "org.Mm.eg.db", "ENTREZID"))

orakeggAJ.dotplot <- dotplot(oraAJkegg, 
                           showCategory = 3, 
                           font.size = 7.5, 
                           label_format = 50,
                           title = "ORA KEGG: AJKO specific Genes") 

oraB6kegg <- enrichKEGG(gene = B6specEntrez,
                        universe = names(unique_genes),
                        organism = "mmu",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.20)
oraB6kegg.res <- as.data.frame(setReadable(oraB6kegg, "org.Mm.eg.db", "ENTREZID"))

orakeggB6.dotplot <- dotplot(oraB6kegg, 
                             showCategory = 3, 
                             font.size = 7.5, 
                             label_format = 50,
                             title = "ORA KEGG: B6KO specific Genes") 
orakeggB6.dotplot

                     


m_t2gC7 <- msigdbr(species = "Mus musculus", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

C7_gsea <- GSEA(unique_genes, TERM2GENE = m_t2gC7)
C7_gsea.readable <- as.data.frame(setReadable(C7_gsea, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))

#Filter for core genes
core.genes <- str_split(as.data.frame(C7_gsea)[,"core_enrichment"], "/")
my.selected.genesC7 <- names(unique_genes[abs(unique_genes) > 2.0])
filtered.core.genesC7 <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genesC7]),
                                paste, collapse = "/")

#Set the GSEA results to be the filtered B6 results
C7_gsea@result$core_enrichment <- filtered.core.genesC7
C7_gsea.filtered <- setReadable(C7_gsea, "org.Mm.eg.db", "ENTREZID")
test.df <- as.data.frame(C7_gsea.filtered)

#Plot the B6 enrichment results
GOCNETC7 <- cnetplot(C7_gsea.filtered, node_label = "gene", foldChange = unique_genes, cex_label_gene = 1.0, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 10)
GOCNETC7.title <- GOCNETC7 + ggtitle("AJKO Immunological Signatures Enrichment CNET")
ggsave("ImmunologicalGSEA_AJKO_CNET.pdf", GOCNETC7.title, width = 16, height = 12)












res.top <- res.ordered[res.ordered$padj < 0.05 & abs(res.ordered$log2FoldChange) > 1,]
res.top <- na.omit(res.top)

#Order the genes by log2FC
res.top <- res.top[order(res.top$log2FoldChange, decreasing = TRUE),]


#Merge the tables and save as a csv
top.degs <- res.top
write.csv(res.top, file = "AJKO_vs_B6KO_filteredDEGS_padj005_log2FC1.csv")


#Heatmap fo top 50 and bottom 50 
#Order top.degs df by decreasing log2FC

print("########## BEGINNING HEATMAP GENERATION ##########")


#Check if there are 100 or more DEGs
if (nrow(top.degs) > 100) {
  print("More than 100 DEGs found: Taking top 50 up-regulated and top 50 down-regulated")
  #Take the top 50 and bottom 50
  total_rows <- nrow(res.top)
  top.degs.keep <- c(1:50, (total_rows - 49):total_rows)
  top.degs.subset <- res.top[top.degs.keep,]
  
  #Pull the baseMean column and Log2FC columns as lists
  top.degs.log2fc <- as.matrix(top.degs.subset$log2FoldChange)
  colnames(top.degs.log2fc) <- "Log2FC"
  top.degs.baseMean <- as.matrix(top.degs.subset$baseMean)
  colnames(top.degs.baseMean) <- "BaseMean"
  
  #Extract just the first column and the counts columns
  top.degs.subset <- top.degs.subset[,c(1,10:56)]
  rownames(top.degs.subset) <- top.degs.subset$Row.names
  top.degs.subset$Row.names <- NULL
  
  #Transpose, center, and scale the normalized counts
  scaled <- t(apply(top.degs.subset, 1, scale))
  colnames(scaled) <- colnames(top.degs.subset)
  
  #Map colors to values
  l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
                              max(top.degs.log2fc)),
                            c("white", "red"))
  
  mean.colors <- colorRamp2(c(quantile(top.degs.baseMean)[1],
                              quantile(top.degs.baseMean)[4]),
                            c("white", "red"))
  
  #Isolate gene symbols for labeling
  labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(top.degs.subset))
  top.degs.subset$Symbols <- labels
  
  
  # Create a named vector for the colors
  color_vector <- c("tan", "black")
  names(color_vector) <- c("A/J", "B6")
  
  #Annotation column
  hmAnno <- HeatmapAnnotation(group = design$Group,
                              name = "", 
                              show_annotation_name = FALSE,
                              col = list(group = color_vector))
  
  
  #Set heatmap splitting pattern
  hmSplit <- rep(1:2, c(24,23))
  
  #Heatmap for scaled data
  hmScaled <- Heatmap(scaled,
                      column_labels = colnames(scaled), 
                      name = "Z-score",
                      cluster_rows = TRUE,
                      cluster_columns = TRUE,
                      top_annotation = hmAnno,
                      column_split = hmSplit,
                      column_title = "A/J vs B6")  
  #Heatmap for log2FC values
  hml2FC <- Heatmap(top.degs.log2fc,
                    row_labels = top.degs.subset$Symbols,
                    cluster_rows = FALSE,
                    name = "log2FC",
                    col = l2FC.colors,
                    cell_fun = function(j,i, x, y, w, h, col) {
                      grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                gp = gpar(fontsize = 3, 
                                          col = "black"))})
  #Heatmap for average expression
  hmMean <- Heatmap(top.degs.baseMean,
                    row_labels = top.degs.subset$Symbols,
                    row_names_gp = gpar(fontsize = 5),
                    cluster_rows = FALSE,
                    name = "Avg Expression",
                    col = mean.colors,
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(round(top.degs.baseMean[i, j],1), x, y,
                                gp = gpar(fontsize = 3,
                                          col = "black"))})
  #Draw the final heatmap
  HM <- hmScaled + hml2FC + hmMean
} 

#########
  #########
  #########
  ##########
  #############
  
  #Wildtype

  raw <- read.csv("AJ_B6_WT_counts.csv", header = TRUE, sep = ",")
  
  #Obtain gene symbols
  raw$Symbol <- mapIds(org.Mm.eg.db, key = raw$Gene, column = "SYMBOL",
                       keytype = "ENSEMBL", multiVals = "first")
  
  #Omit any gene with no symbol annotation or unknown genes
  raw <- raw[!is.na(raw$Symbol),]
  
  #Omit any gene with no MGI annotation or unknown genes and Rikens
  raw <- raw[!grepl("^Gm\\d+$", raw$Symbol),]
  raw <- raw[!grepl("Rik$", raw$Symbol),]
  raw <- raw[!grepl("Rik\\d+$", raw$Symbol),]
  raw <- raw[!grepl("^LOC", raw$Symbol),]
  raw <- raw[!grepl("AK\\d+$", raw$Symbol),]
  raw <- raw[!grepl("AY\\d+$", raw$Symbol),]
  
  
  #Format the rownames
  raw$geneIDs <- paste(raw$Geneid, raw$Symbol, sep = " - ")
  rownames(raw) <- raw$geneIDs
  raw$Geneid<- NULL
  raw$Symbol <- NULL
  raw$geneIDs <- NULL
  
  design <- data.frame(Sample = rep(c("AJWT", "B6WT"),
                                    c(12,12)),
                       Group = rep(c("AJWT", "B6WT"),
                                   c(12,12)))
  rownames(design) <- colnames(raw)
  
  all(colnames(raw) %in% rownames(design))
  all(colnames(raw) == rownames(design))
  
  suppressWarnings({
    dds <- DESeqDataSetFromMatrix(countData = raw,
                                  colData = design,
                                  design = ~ Group)
  })
  
  dds$Group <- relevel(dds$Group, ref = "AJWT")
  
  #QC
  total_counts_QC <- plot_total_counts(dds)
  total_counts <- total_counts_QC + ggtitle("Total Sample Counts")
  library_complexity <- plot_library_complexity(dds)
  total_complexity <- library_complexity + ggtitle("Library Complexity")
  
  variance_stable <- vst(dds)
  
  cowplot::plot_grid(total_counts, total_complexity)
  MSD_plot <- mean_sd_plot(variance_stable)
  MSD_plot + ggtitle("Mean and SD plot")
  
  #Create a summarized experiment
  se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), colData = colData(dds))
  
  #Plot PCA
  PCA <- plotPCA(DESeqTransform(se), intgroup = "Group") +
    geom_text_repel(aes(label = rownames(design))) +
    ggtitle("AJWT vs B6WT PCA") +
    stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
    theme_bw()
  ggsave("AJWT_vs_B6WT_PCA.pdf", PCA)
  
  #Run DESeq2 algorithm
  dds <- DESeq(dds)
  
  setwd("/Users/mikemartinez/Desktop/AJB6_WT")
  
  #View the results
  res <- results(dds)
  res.df <- as.data.frame(res)
  res.df <- na.omit(res.df)
  res.df$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(res.df))
  res.df$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.df))
  
  #Get the normalixed counts
  counts <- counts(dds, normalized = TRUE)
  
  #Order and merge with counts
  res.ordered <- res.df[order(res.df$log2FoldChange, decreasing = TRUE),]
  res.ordered <- merge(res.ordered, counts, by = 0, all = TRUE)
  res.ordered <- na.omit(res.ordered)
  write.csv(res.ordered, file = "AJWTvsB6WT_All_DEGs_orderedByLog2FC_decreasing_noUnknowns.csv")
  
  
  #Get ENSEMBL IDs as their own column to map ENTREZ IDs
  res.ordered$Entrez <- mapIds(org.Mm.eg.db, key = res.ordered$Ensembl,
                               column = "ENTREZID", keytype = "ENSEMBL",
                               multiVals = "first")
  
  #Get a list of the genes
  res.ordered <- res.ordered[order(res.ordered$log2FoldChange, decreasing = TRUE),]
  res.ordered.genes <- res.ordered$log2FoldChange
  
  #Assign Entrez IDs as names for the genes
  names(res.ordered.genes) <- res.ordered$Entrez
  
  #Remove duplicated Entrez IDs and their corresponding values
  unique_entrez_genes <- names(res.ordered.genes[!duplicated(names(res.ordered.genes))])
  unique_genes <- res.ordered.genes[unique_entrez_genes]
  unique_genes <- sort(unique_genes, decreasing = TRUE)
  
  #Run GSEA for GO terms
  #For GSEA: USE ALL THE DEGS BEFORE FILTERING
  gse <- gseGO(unique_genes,
               ont = "all",
               OrgDb = "org.Mm.eg.db",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               eps = 1e-300,
               verbose = TRUE,
               by = "fgsea")
  gse.readable <- as.data.frame(setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
  
  #Set labels
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = "Enriched in B6WT", suppressed = "Enriched in AJWT")
  )
  
  # Create the dotplot with custom facet labels
  GO.dotplot <- dotplot(gse, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 7.5, 
                        label_format = 50,
                        title = "GO Enrichment") + 
    facet_grid(.~.sign, labeller = custom_labels)
  
  
  #Prepare significant genes
  GOR <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
  core.genes <- str_split(as.data.frame(gse)[,"core_enrichment"], "/")
  my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.5])
  
  #For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
  filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                                paste, collapse = "/")
  gse@result$core_enrichment <- filtered.core.genes
  gse.filtered <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
  
  
  GOCNET <- cnetplot(gse.filtered, node_label = "gene", foldChange = unique_genes, cex_label_category = 0.5, cex_label_gene = 0.5, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 5)
  
  ES_plot.GO <- enrichplot::gseaplot2(gse, geneSetID = 1:5, title = "GO Results")
  combined_plot.GO <- wrap_plots(ES_plot.GO, nrow = 3, ncol = 1)
  
  
  #Run KEGG enrichment
  kegg <- gseKEGG(unique_genes,
                  organism = "mmu",
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  eps = 1e-300,
                  verbose = TRUE,
                  by = "fgsea")
  kegg.readable <- as.data.frame(setReadable(kegg, "org.Mm.eg.db", "ENTREZID"))
  
  core.genes <- str_split(as.data.frame(kegg)[,"core_enrichment"], "/")
  filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                                paste, collapse = "/")
  kegg@result$core_enrichment <- filtered.core.genes
  kegg.filtered <- setReadable(kegg, "org.Mm.eg.db", "ENTREZID")
  
  KEGGCNET <- cnetplot(kegg.filtered, node_label = "gene", foldChange = unique_genes, cex_label_category = 0.5, cex_label_gene = 0.5, 
                       colorEdge = TRUE, 
                       circular = TRUE, showCategory = 4)
  
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = "Enriched in B6KO", suppressed = "Enriched in AJKO")
  )
  
  # Create the dotplot with custom facet labels
  KEGG.dotplot <- dotplot(kegg, 
                          showCategory = 15, 
                          split = ".sign", 
                          font.size = 7.5, 
                          label_format = 50,
                          title = "KEGG Enrichment") + 
    facet_grid(.~.sign, labeller = custom_labels)
  
  ES_plot.KEGG <- enrichplot::gseaplot2(kegg, geneSetID = 1:4, title = "KEGG Results")
  combined_plot.KEGG <- wrap_plots(ES_plot.KEGG, nrow = 3, ncol = 1)
  
  
  #Get the AJ and B6 specific genes (from file)
  AJspecific <- read.csv("AJKO.specificGenes_filtered.csv", header = TRUE, sep = ",")
  AJspecificSymbols <- AJspecific$Symbol
  
  #Filter out the AJ-specific genes from the results dataframe
  AJspecGenes <- res.ordered[res.ordered$Symbol %in% AJspecificSymbols, ]
  
  #Get the entezID associated with the AJ-specific genes
  AJspecEntrez <- AJspecGenes$Entrez
  
  #Perform ORA
  oraGO <- enrichGO(gene = AJspecEntrez,
                    universe = names(unique_genes),
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.20,
                    readable = TRUE)
  
  oraGOres <- as.data.frame(oraGO)
  
  oraGO.dotplot <- dotplot(oraGO, 
                           showCategory = 30, 
                           font.size = 7.5, 
                           label_format = 50,
                           title = "ORA: AJWT specific Genes") 
  
  
  B6specific <- read.csv("B6KO.specificGenes_filtered.csv", header = TRUE, sep = ",")
  B6specificSymbols <- B6specific$Symbol
  
  #Filter out the AJ-specific genes from the results dataframe
  B6specGenes <- res.ordered[res.ordered$Symbol %in% B6specificSymbols, ]
  
  #Get the entezID associated with the AJ-specific genes
  B6specEntrez <- B6specGenes$Entrez
  
  #Perform ORA
  oraGOb6 <- enrichGO(gene = B6specEntrez,
                      universe = names(unique_genes),
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.20,
                      readable = TRUE)
  
  oraGOb6res <- as.data.frame(oraGOb6)
  
  oraGOb6.dotplot <- dotplot(oraGOb6, 
                             showCategory = 30, 
                             font.size = 7.5, 
                             label_format = 50,
                             title = "ORA: B6WT specific Genes") 
  
  
  oraAJkegg <- enrichKEGG(gene = AJspecEntrez,
                          universe = names(unique_genes),
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.20)
  oraAJkegg.res <- as.data.frame(setReadable(oraAJkegg, "org.Mm.eg.db", "ENTREZID"))
  
  orakeggAJ.dotplot <- dotplot(oraAJkegg, 
                               showCategory = 3, 
                               font.size = 7.5, 
                               label_format = 50,
                               title = "ORA KEGG: AJKO specific Genes") 
  
  oraB6kegg <- enrichKEGG(gene = B6specEntrez,
                          universe = names(unique_genes),
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.20)
  oraB6kegg.res <- as.data.frame(setReadable(oraB6kegg, "org.Mm.eg.db", "ENTREZID"))
  
  orakeggB6.dotplot <- dotplot(oraB6kegg, 
                               showCategory = 3, 
                               font.size = 7.5, 
                               label_format = 50,
                               title = "ORA KEGG: B6KO specific Genes") 
  orakeggB6.dotplot
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  res.top <- res.ordered[res.ordered$padj < 0.05 & abs(res.ordered$log2FoldChange) > 1,]
  res.top <- na.omit(res.top)
  
  #Order the genes by log2FC
  res.top <- res.top[order(res.top$log2FoldChange, decreasing = TRUE),]
  top.degs <- res.top
  
  #Merge the tables and save as a csv
  write.csv(top.degs, file = "AJWTvsB6WT_top_DEGS_padj005_log2FC1_noUnknowns.csv")
  
  
  #Heatmap fo top 50 and bottom 50 
  #Order top.degs df by decreasing log2FC
  
  print("########## BEGINNING HEATMAP GENERATION ##########")
  
  
  #Check if there are 100 or more DEGs
  if (nrow(top.degs) > 100) {
    print("More than 100 DEGs found: Taking top 50 up-regulated and top 50 down-regulated")
    #Take the top 50 and bottom 50
    total_rows <- nrow(res.top)
    top.degs.keep <- c(1:50, (total_rows - 49):total_rows)
    top.degs.subset <- res.top[top.degs.keep,]
    
    #Pull the baseMean column and Log2FC columns as lists
    top.degs.log2fc <- as.matrix(top.degs.subset$log2FoldChange)
    colnames(top.degs.log2fc) <- "Log2FC"
    top.degs.baseMean <- as.matrix(top.degs.subset$baseMean)
    colnames(top.degs.baseMean) <- "BaseMean"
    
    #Extract just the first column and the counts columns
    top.degs.subset <- top.degs.subset[,c(1,10:34)]
    rownames(top.degs.subset) <- top.degs.subset$Row.names
    top.degs.subset$Row.names <- NULL
    top.degs.subset$Entrez <- NULL
    
    #Transpose, center, and scale the normalized counts
    scaled <- t(apply(top.degs.subset, 1, scale))
    colnames(scaled) <- colnames(top.degs.subset)
    
    #Map colors to values
    l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
                                max(top.degs.log2fc)),
                              c("white", "red"))
    
    mean.colors <- colorRamp2(c(quantile(top.degs.baseMean)[1],
                                quantile(top.degs.baseMean)[4]),
                              c("white", "red"))
    
    #Isolate gene symbols for labeling
    labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(top.degs.subset))
    top.degs.subset$Symbols <- labels
    
    
    # Create a named vector for the colors
    color_vector <- c("tan", "black")
    names(color_vector) <- c("AJWT", "B6WT")
    
    #Annotation column
    hmAnno <- HeatmapAnnotation(group = design$Group,
                                name = "", 
                                show_annotation_name = FALSE,
                                col = list(group = color_vector))
    
    
    #Set heatmap splitting pattern
    hmSplit <- rep(1:2, c(12, 12))
    
    #Heatmap for scaled data
    hmScaled <- Heatmap(scaled,
                        column_labels = colnames(scaled), 
                        name = "Z-score",
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        top_annotation = hmAnno,
                        column_split = hmSplit,
                        column_title = "AJWT vs B6WT top Genes")  
    #Heatmap for log2FC values
    hml2FC <- Heatmap(top.degs.log2fc,
                      row_labels = top.degs.subset$Symbols,
                      cluster_rows = FALSE,
                      name = "log2FC",
                      col = l2FC.colors,
                      cell_fun = function(j,i, x, y, w, h, col) {
                        grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                  gp = gpar(fontsize = 3, 
                                            col = "black"))})
    #Heatmap for average expression
    hmMean <- Heatmap(top.degs.baseMean,
                      row_labels = top.degs.subset$Symbols,
                      row_names_gp = gpar(fontsize = 4.5),
                      cluster_rows = FALSE,
                      name = "Avg Expression",
                      col = mean.colors,
                      cell_fun = function(j, i, x, y, w, h, col) {
                        grid.text(round(top.degs.baseMean[i, j],1), x, y,
                                  gp = gpar(fontsize = 3,
                                            col = "black"))})
    #Draw the final heatmap
    HM <- hmScaled + hml2FC + hmMean
    


