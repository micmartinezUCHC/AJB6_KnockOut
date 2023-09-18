#The purpose of this script is to analyze normal to tumor samples for each treatment group
#Mike Martinez v 2.0.1


#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut//")

#Set command line arguments
args <- commandArgs(trailingOnly = TRUE)

######### USERS GUIDE #########
#To run this script, supply the following arguments in this order:
#1: path to counts.csv file
#2: analysis name for generating output files
#3: Group 1 name for generating DESeq2 design table and figures
#4: Number of obsrevations of group 1 in the counts matrix
#5: Group 2 name
#6: Number of observations of group 2 in the counts matrix

###############################

#Check to make sure all arguments are supplied
if (length(args) < 6) {
  stop("At least one argument is missing\n 
  -input counts file\n 
  -analysis basename\n 
  -group 1 name\n 
  -number of group 1 observations\n 
  -group 2 name\n
  -number of group 2 observations")
}

#Load libraries
suppressWarnings({
  #Data manipulation
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(dplyr))
  
  #Differential gene expression
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(ashr))
  
  #Graphics and visualizations
  suppressPackageStartupMessages(library(ggplot2))
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



print("########## BEGINNING PROCESSING ##########")
print(paste("Counts file provided", args[1], sep = " "))

#Read in the data
raw <- read.csv(args[1], header = TRUE, sep = ",")

#Obtain gene symbols
raw$Symbol <- mapIds(org.Mm.eg.db, key = raw$Gene, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no symbol annotation or unknown genes
raw <- raw[!is.na(raw$Symbol),]

#Omit any gene with no MGI annotation or unknown genes and Rikens
raw <- raw[!grepl("^Gm\\d+$", raw$Symbol),]

print("Data successfully loaded")
print(paste("Number of columns in raw data = ", ncol(raw), sep = " "))
print(colnames(raw))

#Format the rownames to join Ensembl ID and symbol
raw$geneIDs <- paste(raw$Geneid, raw$Symbol, sep = " - ")
print(nrow(raw))
rownames(raw) <- raw$geneIDs
raw$X <- NULL
raw$Symbol <- NULL
raw$geneIDs <- NULL
raw$Geneid<- NULL
print("Gene symbols successfully acquired")

#Filter
  #minimumCountpergene <- 10
  #MinSampleWithminimumgeneCounts <- 5

#Filter out low read counts for Normal vs control, Naproxen vs control, etc...
  #raw <- raw[rowSums(data.frame(raw>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

print("########## BEGINNING DESEQ2 ##########")
print(paste("Reference group:", args[3], sep = " "))
print(paste("2nd Group:", args[5], sep = " "))



#Create a design table
print("##### SAMPLE NAMES #####")
print(colnames(raw))

print("GENERATING DESIGN TABLE")
design <- data.frame(Sample = rep(c(args[3], args[5]),
                                  c(args[4],args[6])),
                     Group = rep(c(args[3], args[5]),
                                 c(args[4],args[6])))
rownames(design) <- colnames(raw)
print(design)

#Sanity check: are all colnames in raw present in the design table?
#are the colnames in raw in the same order as rownames in design?
print("Check: all colnames in counts present in design table")
all(colnames(raw) %in% rownames(design))

print("Check; all colnames in counts in same order as design table")
all(colnames(raw) == rownames(design))

suppressWarnings({
  print("Running Deseq2 Algorithm")
  dds <- DESeqDataSetFromMatrix(countData = raw,
                                colData = design,
                                design = ~ Group)
  
  #Set a factor level
  dds$Group <- relevel(dds$Group, ref = args[3])
})

#Run DESeq2 algorithm
dds <- DESeq(dds)
print("DESeq2 successful !!!!!")


#View the results
res <- results(dds)
res.df <- as.data.frame(res)
res.df <- na.omit(res.df)

labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.df))
res.df$Symbols <- labels

res.df <- res.df[!grepl("Rik$", res.df$Symbols),]
res.df <- na.omit(res.df)


#Get the normalixed counts
counts <- counts(dds, normalized = TRUE)

res.ordered <- res.df[order(res.df$log2FoldChange, decreasing = TRUE),]
res.ordered <- merge(res.ordered, counts, by = 0, all = TRUE)
res.ordered$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(res.ordered))
res.ordered$Symbol <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.ordered))
write.csv(res.ordered, file = paste(args[2], "All_DEGs_orderedByLog2FC_decreasing.csv", sep = "_"))


#Principal Component Analysis

print("########## PLOTTING PRINCIPAL COMPONENTS ##########")

#Create a summarized experiment
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), colData = colData(dds))

#Plot PCA
PCA <- plotPCA(DESeqTransform(se), intgroup = "Group") +
  geom_text_repel(aes(label = rownames(design))) +
  ggtitle(args[2]) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group))
ggsave(paste(args[2], "PCA.pdf", sep = "_"), PCA)


#GO and KEGG analysis

print("########## BEGINNING GSEA: GO AND KEGG ##########")
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

suppressWarnings({
  #Run GSEA for GO terms
  #For GSEA: USE ALL THE DEGS BEFORE FILTERING
  gse <- gseGO(unique_genes,
               ont = "bp",
               OrgDb = "org.Mm.eg.db",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               eps = 1e-300,
               verbose = TRUE,
               by = "fgsea")
  
  if (length(gse@result$Description > 0)) {
    gse.readable <- as.data.frame(setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
    write.csv(gse.readable, file = paste(args[2], "GO_enrichment_results.csv", sep = "_"))
    
    # Custom facet label mapping
    custom_labels <- labeller(
      .default = label_value,
      .sign = c(activated = paste("up in", args[5], sep = " "), suppressed = paste("up in", args[3], sep = " "))
    )
    
    # Create the dotplot with custom facet labels
    GO.dotplot <- dotplot(gse, 
                          showCategory = 15, 
                          split = ".sign", 
                          font.size = 7.5, 
                          label_format = 50,
                          title = paste(args[2], "KEGG Enrichment", sep = " ")) + 
      facet_grid(.~.sign, labeller = custom_labels)
    ggsave(paste(args[2], "GO_enrichment_dotplot.pdf", sep = "_"), GO.dotplot)
    
    #Plot biological network plot
    GOR <- setReadable(gse, "org.Mm.eg.db", "ENTREZID")
    GOCNET <- cnetplot(GOR, node_label = "all", foldChange = unique_genes, cex_label_category = 0.2, cex_label_gene = 0.2, colorEdge = TRUE,
                       circular = FALSE)
    ggsave(paste(args[2], "GO_cnet.pdf", sep = "_"), GOCNET)
    
    #Plot Running enrichment score plot
    ES_plot.GO <- enrichplot::gseaplot2(gse, geneSetID = 1:5, title = "GO Results")
    combined_plot.GO <- wrap_plots(ES_plot.GO, nrow = 3, ncol = 1)
    ggsave(paste(args[2], "GO_ES_plot.pdf", sep = "_"), combined_plot.GO)
    
  } else {
    print("No GO terms enriched")
  }
  
  #Run KEGG enrichment
  kegg <- gseKEGG(unique_genes,
                  organism = "mmu",
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  eps = 1e-300,
                  verbose = TRUE,
                  by = "fgsea")
  
  
  
  
  if (length(kegg@result$Description > 0)) {
    kegg.readable <- as.data.frame(setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
    write.csv(kegg.readable, file = paste(args[2], "KEGG_enrichment_results.csv", sep = "_"))
    
    # Custom facet label mapping
    custom_labels <- labeller(
      .default = label_value,
      .sign = c(activated = paste("up in", args[5], sep = " "), suppressed = paste("up in", args[3], sep = " "))
    )
    
    # Create the dotplot with custom facet labels
    kegg.dotplot <- dotplot(kegg, 
                            showCategory = 15, 
                            split = ".sign", 
                            font.size = 7.5, 
                            label_format = 50,
                            title = paste(args[2], "KEGG Enrichment", sep = " ")) + 
      facet_grid(.~.sign, labeller = custom_labels)
    ggsave(paste(args[2], "KEGG_enrichment_dotplot.pdf", sep = "_"), kegg.dotplot)
    
    #Plot biological network plot
    keggR <- setReadable(kegg, "org.Mm.eg.db", "ENTREZID")
    keggCNET <- cnetplot(keggR, node_label = "all", foldChange = unique_genes, cex_label_category = 0.2, cex_label_gene = 0.2, colorEdge = TRUE,
                         circular = FALSE)
    ggsave(paste(args[2], "KEGG_cnet.pdf", sep = "_"), keggCNET)
    
    #Plot Running enrichment score plot
    ES_plot.kegg <- enrichplot::gseaplot2(kegg, geneSetID = 1:5, title = "KEGG Results")
    combined_plot.kegg <- wrap_plots(ES_plot.kegg, nrow = 3, ncol = 1)
    ggsave(paste(args[2], "KEGG_ES_plot.pdf", sep = "_"), combined_plot.kegg)
    
    
  } else {
    print("No Kegg terms enriched")
  }
  
  
  
  
  
  
  
  
  
  #######ADDITIONAL UNIVERSAL GENE SET ANALYSES
  #H: HALLMARK GENE SETS
  #C1: POSITIONAL GENE SETS
  #C2: CURATED GENE SETS
  #C3: MOTIF GENE SETS
  #C4: COMPUTATIOANL GENE SETS
  #C5: GO GENE SETS
  #C6: ONCOGENIC SIGNATUES
  #C7: IMMUNOLOGIC SIGNATUTES
  
  #Specify the term to gene mapping and category for C6 Oncology signatures and C7 Immune Signatures
  m_t2gC6 <- msigdbr(species = "Mus musculus", category = "C6") %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2gC7 <- msigdbr(species = "Mus musculus", category = "C7") %>%
    dplyr::select(gs_name, entrez_gene)
  
  print("########## BEGINNING GSEA: ONCOGENIC SIGNATURES ##########")
  C6_gsea <- GSEA(unique_genes, TERM2GENE = m_t2gC6)
  C6_gsea.readable <- as.data.frame(setReadable(C6_gsea, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
  
  if (length(C6_gsea@result$Description > 0)) {
    write.csv(C6_gsea.readable, file = paste(args[2], "C6_gsea_results.csv", sep = "_"))
    
    #Custom labels
    custom_labels <- labeller(
      .default = label_value,
      .sign = c(activated = paste("up in", args[5], sep = " "), suppressed = paste("up in", args[3], sep = " "))
    )
    
    # Create the dotplot with custom facet labels
    C6.dotplot <- dotplot(C6_gsea, 
                          showCategory = 15, 
                          split = ".sign", 
                          font.size = 7.5, 
                          label_format = 50,
                          title = paste(args[2], "Oncological Enrichment", sep = " ")) + 
      facet_grid(.~.sign, labeller = custom_labels)
    ggsave(paste(args[2], "C6_enrichment_dotplot.pdf", sep = "_"), C6.dotplot)
  } else {
    print("No Oncogenic Signatures enriched")
  }
  
  print("########## BEGINNING GSEA: IMMUNOLOGIC SIGNATURES ##########")
  C7_gsea <- GSEA(unique_genes, TERM2GENE = m_t2gC7)
  C7_gsea.readable <- as.data.frame(setReadable(C7_gsea, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
  
  if (length(C7_gsea@result$Description > 0)) {
    write.csv(C7_gsea.readable, file = paste(args[2], "C7_gsea_results.csv", sep = "_"))
    
    # Custom facet label mapping
    custom_labels <- labeller(
      .default = label_value,
      .sign = c(activated = paste("up in", args[5], sep = " "), suppressed = paste("up in", args[3], sep = " "))
    )
    
    # Create the dotplot with custom facet labels
    C7.dotplot <- dotplot(C7_gsea, 
                          showCategory = 15, 
                          split = ".sign", 
                          font.size = 7.5, 
                          label_format = 50,
                          title = paste(args[2], "Immunological Enrichment", sep = " ")) + 
      facet_grid(.~.sign, labeller = custom_labels)
    ggsave(paste(args[2], "C7_enrichment_dotplot.pdf", sep = "_"), C7.dotplot)
    
  } else {
    print("No Immunologic Signatures enriched")
  }
}) #####END OF SUPPRESSWARNINGS BLOCK

#Filter the data frame based on padj and log2FC
res.top <- res.df[res.df$padj < 0.05 & abs(res.df$log2FoldChange) > 1,]
res.top <- na.omit(res.top)

#Order the genes by log2FC
res.top <- res.top[order(res.top$log2FoldChange, decreasing = TRUE),]

#Get the normalized counts associated with the top degs
counts.keep <- rownames(res.top)
deg.counts <- counts[counts.keep,]

#Merge the tables and save as a csv
top.degs <- merge(res.top, deg.counts, by = 0, all = TRUE)
write.csv(top.degs, file = paste(args[2],"top_DEGS.csv", sep = "_"))


#Heatmap fo top 50 and bottom 50 
#Order top.degs df by decreasing log2FC
top.degs <- top.degs[order(top.degs$log2FoldChange, decreasing = TRUE),]

print("########## BEGINNING HEATMAP GENERATION ##########")


#Check if there are 100 or more DEGs
if (nrow(top.degs) > 100) {
  print("More than 100 DEGs found: Taking top 50 up-regulated and top 50 down-regulated")
  #Take the top 50 and bottom 50
  total_rows <- nrow(top.degs)
  top.degs.keep <- c(1:50, (total_rows - 49):total_rows)
  top.degs.subset <- top.degs[top.degs.keep,]
  
  #Pull the baseMean column and Log2FC columns as lists
  top.degs.log2fc <- as.matrix(top.degs.subset$log2FoldChange)
  colnames(top.degs.log2fc) <- "Log2FC"
  top.degs.baseMean <- as.matrix(top.degs.subset$baseMean)
  colnames(top.degs.baseMean) <- "BaseMean"
  
  #Extract just the first column and the counts columns
  top.degs.subset <- top.degs.subset[,c(1,8:ncol(top.degs.subset))]
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
  names(color_vector) <- c(args[3], args[5])
  
  #Annotation column
  hmAnno <- HeatmapAnnotation(group = design$Group,
                              name = "", 
                              show_annotation_name = FALSE,
                              col = list(group = color_vector))
  
  
  #Set heatmap splitting pattern
  hmSplit <- rep(1:2, c(args[4], args[6]))
  
  #Heatmap for scaled data
  hmScaled <- Heatmap(scaled,
                      column_labels = colnames(scaled), 
                      name = "Z-score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      top_annotation = hmAnno,
                      column_split = hmSplit,
                      column_title = paste(args[2],"heatmap", sep = " "))  
  #Heatmap for log2FC values
  hml2FC <- Heatmap(top.degs.log2fc,
                    row_labels = rownames(top.degs.subset),
                    cluster_rows = FALSE,
                    name = "log2FC",
                    col = l2FC.colors,
                    cell_fun = function(j,i, x, y, w, h, col) {
                      grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                gp = gpar(fontsize = 3, 
                                          col = "black"))})
  #Heatmap for average expression
  hmMean <- Heatmap(top.degs.baseMean,
                    row_labels = rownames(top.degs.subset),
                    row_names_gp = gpar(fontsize = 3),
                    cluster_rows = FALSE,
                    name = "Avg Expression",
                    col = mean.colors,
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(round(top.degs.baseMean[i, j],1), x, y,
                                gp = gpar(fontsize = 3,
                                          col = "black"))})
  #Draw the final heatmap
  HM <- hmScaled + hml2FC + hmMean
  pdf(file = paste(args[2], "heatmap.pdf", sep = "."))
  
  #Draw the complex heatmap
  draw(HM)
  #Close the PDF device
  dev.off()
} else {
  print("There is less than 100 DEGs present! Using all DEGs to generate heatmap")
  
  #Redefine top degs
  top.degs <- merge(res.top, deg.counts, by = 0, all = TRUE)
  top.degs <- top.degs[order(top.degs$log2FoldChange, decreasing = TRUE),]
  print(paste("There are", nrow(top.degs), "DEGs present after filtering", sep = " "))
  print(top.degs$Row.names)
  
  #Pull the baseMean column and Log2FC columns as lists
  top.degs.log2fc <- as.matrix(top.degs$log2FoldChange)
  colnames(top.degs.log2fc) <- "Log2FC"
  top.degs.baseMean <- as.matrix(top.degs$baseMean)
  colnames(top.degs.baseMean) <- "BaseMean"
  
  #Extract just the first column and the counts columns
  top.degs.subset2 <- top.degs[,c(1,8:ncol(top.degs))]
  rownames(top.degs.subset2) <- top.degs.subset2$Row.names
  top.degs.subset2$Row.names <- NULL
  
  #Transpose, center, and scale the normalized counts
  scaled <- t(apply(top.degs.subset2, 1, scale))
  colnames(scaled) <- colnames(top.degs.subset2)
  
  #Map colors to values
  l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
                              max(top.degs.log2fc)),
                            c("white", "red"))
  
  mean.colors <- colorRamp2(c(quantile(top.degs.baseMean)[1],
                              quantile(top.degs.baseMean)[4]),
                            c("white", "red"))
  
  #Isolate gene symbols for labeling
  labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(top.degs.subset2))
  top.degs.subset2$Symbols <- labels
  
  # Create a named vector for the colors
  color_vector <- c("tan", "black")
  names(color_vector) <- c(args[3], args[5])
  
  #Annotation column
  hmAnno <- HeatmapAnnotation(group = design$Group,
                              name = "", 
                              show_annotation_name = FALSE,
                              col = list(group = color_vector))
  
  
  
  #Set heatmap splitting pattern
  hmSplit <- rep(1:2, c(args[4], args[6]))
  
  #Heatmap for scaled data
  hmScaled <- Heatmap(scaled,
                      column_labels = colnames(scaled), 
                      name = "Z-score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      top_annotation = hmAnno,
                      column_split = hmSplit,
                      column_title = paste(args[2],"heatmap", sep = " "))  
  
  #Heatmap for log2FC values
  hml2FC <- Heatmap(top.degs.log2fc,
                    row_labels = rownames(top.degs.subset2),
                    cluster_rows = FALSE,
                    name = "log2FC",
                    col = l2FC.colors,
                    cell_fun = function(j,i, x, y, w, h, col) {
                      grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                gp = gpar(fontsize = 3, 
                                          col = "black"))})
  #Heatmap for average expression
  hmMean <- Heatmap(top.degs.baseMean,
                    row_labels = rownames(top.degs.subset2),
                    row_names_gp = gpar(fontsize = 3),
                    cluster_rows = FALSE,
                    name = "Avg Expression",
                    col = mean.colors,
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(round(top.degs.baseMean[i, j],1), x, y,
                                gp = gpar(fontsize = 3,
                                          col = "black"))})
  #Draw the final heatmap
  HM <- hmScaled + hml2FC + hmMean
  pdf(file = paste(args[2], "heatmap.pdf", sep = "."))
  
  #Draw the complex heatmap
  draw(HM)
  #Close the PDF device
  dev.off()
  
}

print("Operation completed!")













