#The purpose of this script is to conduct all Gene Set Enrichment analyses for Dr. Nakanishi's AJ/B6 Project

#Input files required: list of ALL DEGs generated from DESeq2 with AJ as the reference (denominator) level

#Load Libraries
#GSEA analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(enrichplot))

#Data manipulation
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))


#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/")

#Read in the input file
DEGs <- read.csv("AJKO_vs_B6KO_All_DEGs_NoUnknowns.csv")
res.ordered <- DEGs

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
gse.readable <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gseDF <- as.data.frame(gse.readable)
gseDF <- gseDF[order(gseDF$NES, decreasing = TRUE), ]

#Set labels
custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in B6", suppressed = "Enriched in AJ")
)

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

######  AJKO GSEA CNET
#Prepare significant genes
AJ_gsea <- gse
AJ_core.genes <- str_split(as.data.frame(AJ_gsea)[,"core_enrichment"], "/")
my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.5])

#For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
filtered.core.genes <- sapply(lapply(AJ_core.genes, function(x) x[x %in% my.selected.genes]),
                              paste, collapse = "/")
AJ_gsea@result$core_enrichment <- filtered.core.genes
AJ_gse.filtered <- setReadable(AJ_gsea, "org.Mm.eg.db", "ENTREZID")
AJreadable <- as.data.frame(AJ_gse.filtered)
AJreadable <- AJreadable[order(AJreadable$NES, decreasing = TRUE), ]

AJ_GOCNET <- cnetplot(AJ_gse.filtered, node_label = "gene", foldChange = unique_genes, cex_label_gene = 0.6, 
                   colorEdge = TRUE, 
                   circular = TRUE, showCategory = 10)
AJ_GOCNET.title <- AJ_GOCNET + ggtitle("AJKO CNET")

#FOR B6#
#Take the subset of the gse results that have + NES score
B6_gsea <- gse
B6_gse.df <- as.data.frame(B6_gsea)
B6gse.pos <- B6_gse.df[B6_gse.df$NES > 0,]
B6_gsea@result <- B6gse.pos

#Filter for core genes
B6core.genes <- str_split(as.data.frame(B6_gsea)[,"core_enrichment"], "/")
my.selected.genesB6 <- names(unique_genes[abs(unique_genes) > 2.0])
filtered.core.genesB6 <- sapply(lapply(B6core.genes, function(x) x[x %in% my.selected.genesB6]),
                                paste, collapse = "/")

#Set the GSEA results to be the filtered B6 results
B6_gsea@result$core_enrichment <- filtered.core.genesB6
B6_gse.filteredB6 <- setReadable(B6_gsea, "org.Mm.eg.db", "ENTREZID")
B6readable <- as.data.frame(B6_gse.filteredB6)

#Plot the B6 enrichment results
B6_GOCNET <- cnetplot(B6_gse.filteredB6, node_label = "gene", foldChange = unique_genes, cex_label_gene = 0.6, 
                     colorEdge = TRUE, 
                     circular = TRUE, showCategory = 10)
B6_GOCNET.title <- B6_GOCNET + ggtitle("B6KO CNET")

m_t2gC7 <- msigdbr(species = "Mus musculus", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

C7_gsea <- GSEA(unique_genes, TERM2GENE = m_t2gC7)
C7_gsea.readable <- as.data.frame(setReadable(C7_gsea, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in B6:KO", suppressed = "Enriched in A/J:KO")
)

# Create the dotplot with custom facet labels
C7.dotplot <- dotplot(C7_gsea, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 7.5, 
                        label_format = 90,
                        title = "Immunological Signatures Enrichment") + 
  facet_grid(.~.sign, labeller = custom_labels)
ggsave("AJKO_ImmunologicalDotplot.png", C7.dotplot)


#Filter for core genes
C7_core.genes <- str_split(as.data.frame(C7_gsea)[,"core_enrichment"], "/")
my.selected.genesC7 <- names(unique_genes[abs(unique_genes) > 2.0])
filtered.core.genesC7 <- sapply(lapply(C7_core.genes, function(x) x[x %in% my.selected.genesC7]),
                                paste, collapse = "/")

#Set the GSEA results to be the filtered B6 results
C7_gsea@result$core_enrichment <- filtered.core.genesC7
C7_gsea.filtered <- setReadable(C7_gsea, "org.Mm.eg.db", "ENTREZID")
test.df <- as.data.frame(C7_gsea.filtered)

#Plot the B6 enrichment results
GOCNETC7 <- cnetplot(C7_gsea.filtered, node_label = "gene", foldChange = unique_genes ,cex_label_gene = 1.0, 
                     colorEdge = TRUE, layout = "kk",
                     circular = FALSE, showCategory = 10)
GOCNETC7.title <- GOCNETC7 + ggtitle("AJKO Immunological Signatures Enrichment CNET")
ggsave("ImmunologicalGSEA_AJKO_CNET.pdf", GOCNETC7.title, width = 16, height = 12)

GOCNETC7


library("DOSE")
gse2 <- simplify(gse)
cnetplot(gse, foldChange=unique_genes)

B6_heatplot <- heatplot(B6_gse.filteredB6, foldChange = unique_genes)
AJ_heatplot <- heatplot(AJ_gse.filtered, foldChange = unique_genes)

cowplot::plot_grid(AJ_heatplot, B6_heatplot)

B6 <- heatplot(B6_gse.filteredB6, symbol = "rect", foldChange = unique_genes, pvalue = NULL, label_format = 50)
B6 <- B6 + ggtitle("B6")
AJ <- heatplot(AJ_gse.filtered, symbol = "rect", foldChange = unique_genes, pvalue = NULL, label_format = 50)
AJ <- AJ + ggtitle("AJ")
  
merged_plot <- cowplot::plot_grid(AJ, B6)
merged_plot <- merged_plot +
  theme(axis_text.x = element_text(size = 6))
ggsave("AJKO_B6KO_GSEA_heatplot.pdf", merged_plot, width = 20, height = 12)


#Merge
filteredGSE <- rbind(AJreadable, B6readable)
gse.readable@result$core_enrichment <- filteredGSE


#FOR B6#
#Take the subset of the gse results that have + NES score
B6_gsea <- gse
B6_gse.df <- as.data.frame(B6_gsea)
B6gse.pos <- B6_gse.df[B6_gse.df$NES > 0,]
B6_gsea@result <- B6gse.pos

#Filter for core genes
B6core.genes <- str_split(as.data.frame(B6_gsea)[,"core_enrichment"], "/")
my.selected.genesB6 <- names(unique_genes[abs(unique_genes) > 2.0])
filtered.core.genesB6 <- sapply(lapply(B6core.genes, function(x) x[x %in% my.selected.genesB6]),
                                paste, collapse = "/")

#Set the GSEA results to be the filtered B6 results
B6_gsea@result$core_enrichment <- filtered.core.genesB6
B6_gse.filteredB6 <- setReadable(B6_gsea, "org.Mm.eg.db", "ENTREZID")
B6readable <- as.data.frame(B6_gse.filteredB6)


AJ_gsea <- gse
AJ_gse.df <- as.data.frame(AJ_gsea)
AJgse.pos <- B6_gse.df[AJ_gse.df$NES < 0,]
AJ_gsea@result <- AJgse.pos

AJ_gsea <- gse
AJ_core.genes <- str_split(as.data.frame(AJ_gsea)[,"core_enrichment"], "/")
my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.5])

#For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
filtered.core.genes <- sapply(lapply(AJ_core.genes, function(x) x[x %in% my.selected.genes]),
                              paste, collapse = "/")
AJ_gsea@result$core_enrichment <- filtered.core.genes
AJ_gse.filtered <- setReadable(AJ_gsea, "org.Mm.eg.db", "ENTREZID")
AJreadable <- as.data.frame(AJ_gse.filtered)


B6 <- heatplot(B6_gse.filteredB6, symbol = "rect", foldChange = unique_genes, pvalue = NULL, label_format = 50)
B6 <- B6 + ggtitle("Enriched in B6")
AJ <- heatplot(AJ_gse.filtered, symbol = "rect", foldChange = unique_genes, pvalue = NULL, label_format = 50)
AJ <- AJ + ggtitle("Enriched in AJ")

#Reformat x axis labels
AJtest <- AJ + theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6))
B6test <- B6 + theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6))

mergedPlot <- cowplot::plot_grid(AJtest, B6test)
ggsave("AJKO_B6KO_GSEA_heatplot.pdf", mergedPlot, width = 20, height = 12)


kegg <- gseKEGG(unique_genes,
                organism = "mmu",
                keyType = "kegg",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 1e-300,
                verbose = TRUE,
                by = "fgsea")
kegg.readable <- as.data.frame(setReadable(kegg, "org.Mm.eg.db", "ENTREZID"))

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in B6:KO", suppressed = "Enriched in A/J:KO")
)

kegg@result$Description <- gsub(" - Mus musculus (house mouse$", "", kegg@result$Description)

#Remove the "-Mus musculus (house mouse)" part from every KEGG description
desc <- kegg@result$Description
desc <- sub(" - .*", "", desc)
kegg@result$Description <- desc



# Create the dotplot with custom facet labels
KEGG.dotplot <- dotplot(kegg, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 7.5, 
                        label_format = 90,
                        title = "A/J:KO Enriched KEGG Terms") +
  facet_grid(.~.sign, labeller = custom_labels)
KEGG.dotplot
ggsave("KO_keggDotplot.png", KEGG.dotplot)



###





#Read in the original KEGG results so I can change the labels
oldKegg <- read.csv("/users/mikemartinez/Desktop/AJB6_KnockOut/Results/AJB6__fullAnalysis/AJKO_vs_B6KO/GeneSet_EnrichmentAnalysis/KEGG/AJKO_vs_B6KO_KEGG_enrichment_results.csv", header = TRUE, sep = ",")
kegg@result <- oldKegg

oldC7 <- read.csv("/users/mikemartinez/Desktop/AJB6_KnockOut/Results/AJB6__fullAnalysis/AJKO_vs_B6KO/GeneSet_EnrichmentAnalysis/MSigDB_ImmunologicSignatures/AJKO_vs_B6KO_C7_gsea_results.csv", header = TRUE, sep = ",")
C7_gsea@result <- oldC7

C7.dotplot <- dotplot(C7_gsea, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 7.5, 
                        label_format = 90,
                        title = "MSigDB Immunological Signatures") +
  facet_grid(.~.sign, labeller = custom_labels)
C7.dotplot
ggsave("AJKO_B6KO_ImmuneSigs_Dotplot.png", C7.dotplot)

















