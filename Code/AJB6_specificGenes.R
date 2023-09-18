#The purpose of this script is to deterine the AJ and B6 specific genes in the AJ/B6 Project
#The DESeq2 analysis was already conducted on ALL AJ and B6 strains (WT and KO) and is saved under the file
  #AJB6_AllvsAll_Top_results.csv where DEGs were filtered padj < 0.05 and |log2FC| > 2
#This script aims to do the following: find genes where almost all the reads for AJ are 0
  #These genes will be designated as B6 specific.
#The same principles will apply to B6 genes with low read counts (AJ-specific)
###################################################################################################################

#Load libraries
#GSEA analysis
library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationDbi")

#Plotting
library("ggplot2")
library("ggrepel")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")


#Heatmap
library("pheatmap")


#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/")

#Read in the input file (DEG results from All vs All)
res <- read.csv("RNASeq/AJB6_AllvsAll_top_results.csv", header = TRUE, sep = ",")
  
#Clean up the data frame a little bit
res <- res[,c(2:4,8:ncol(res))]
  #Columns 7:ncol(res) are the columns that contain counts data
  #7:30 are AJ (there are 24 samples)
  #31:53 are B6 (there are 23 samples)

#Find the genes that have low read counts in AJ
minimumCountpergene <- 0
MinSampleWithminimumgeneCounts <- 12

#Low read counts in AJ counts
B6_specific <- res[rowSums(data.frame(res[,7:18] = minimumCountpergene)) > MinSampleWithminimumgeneCounts
                   & rowSums(data.frame(res[,31:41] > minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#Low read counts in B6 counts
AJ_specific <- res[rowSums(data.frame(res[,31:53] < minimumCountpergene)) > MinSampleWithminimumgeneCounts
                   & rowSums(data.frame(res[,7:30] > minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#Filter out Rikens and other unknowns
B6_specific <- B6_specific[!grepl("[A-Za-z0-9]+Rik", B6_specific$Symbols),]
AJ_specific <- AJ_specific[!grepl("[A-Za-z0-9]+Rik", AJ_specific$Symbols),]
AJ_specific <- AJ_specific[!grepl("AA465934", AJ_specific$Symbols),]
AJ_specific <- AJ_specific[!grepl("AK157302", AJ_specific$Symbols),]

#Write csv
write.csv(B6_specific, file = "B6_specificGenes.csv")
write.csv(AJ_specific, file = "AJ_specificGenes.csv")


AJ_spec_counts <- AJ_specific[,c(1:4,7:ncol(AJ_specific))]
AJ_spec_counts <- AJ_spec_counts[order(AJ_spec_counts$log2FoldChange, decreasing = TRUE),]
B6_spec_counts <- B6_specific[,c(1:4,7:ncol(B6_specific))]
B6_spec_counts <- B6_spec_counts[order(B6_spec_counts$log2FoldChange, decreasing = TRUE),]

AJB6 <- rbind(AJ_spec_counts, B6_spec_counts)

AJ_counts <- AJ_spec_counts[,5:ncol(AJ_spec_counts)]
rownames(AJ_counts) <- AJ_spec_counts$Symbols
AJ_mat <- as.matrix(AJ_counts)

B6_counts <- B6_spec_counts[,5:ncol(B6_spec_counts)]
rownames(B6_counts) <- B6_spec_counts$Symbols
B6_mat <- as.matrix(B6_counts)


#Split the data frames into just AJ and B6
AJ <- res[,c(1:30)]
B6 <- res[,c(1:6,31:ncol(res))]


#Find the genes that have low read counts in AJ
minimumCountpergene <- 1
MinSampleWithminimumgeneCounts <- 11

#Low read counts in AJ counts
B6_test <- AJ[rowSums(data.frame(AJ[,7:18] < minimumCountpergene)) > MinSampleWithminimumgeneCounts
              & rowSums(data.frame(AJ[,19:30] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]
B6_test.filt <- B6_test[!grepl("[A-Za-z0-9]+Rik", B6_test$Symbols),]
B6_specific_symbols <- B6_test.filt$Symbols
B6_specific <- B6[B6$Symbols %in% B6_specific_symbols,]
B6_specific_genes <- merge(B6_test.filt, B6_specific[,7:ncol(B6_specific)], by = 0, all = TRUE)

#Find the genes that have low read counts in B6
minimumCountpergene <- 1
MinSampleWithminimumgeneCounts <- 11
AJ_test <- B6[rowSums(data.frame(B6[,7:17] < minimumCountpergene)) > MinSampleWithminimumgeneCounts,]# &
AJ_test.filt <- AJ_test[!grepl("[A-Za-z0-9]+Rik", AJ_test$Symbols),]
AJ_specific_symbols <- AJ_test.filt$Symbols
AJ_specific <- AJ[AJ$Symbols %in% AJ_specific_symbols,]
AJ_specific_genes <- merge(AJ_test.filt, AJ_specific[,7:ncol(AJ_specific)], by = 0, all = TRUE)

#Write csv
write.csv(B6_specific_genes, file = "B6_specificGenes_new.csv")
write.csv(AJ_specific_genes, file = "AJ_specificGenes_new.csv")


###



