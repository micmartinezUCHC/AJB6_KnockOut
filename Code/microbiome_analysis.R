#The purpose of this script is to analyze the AJ/B6 Microbiome data from Shoreline

#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/Microbiome/")

#Load libraries
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ggh4x")
library("tidyverse")
library("RColorBrewer")
###################
###################

######################
#Read in the RAW microbiome counts
  #test <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/phylumData.csv", header = TRUE, sep = ",")
  #names <- test$taxon

  #sums <- rowSums(test[,2:ncol(test)])
  #relabund2 <- test[,2:ncol(test)]/sums

  #rownames(relabund2) <- names
  #relabund2 <- rownames_to_column(relabund2, var = "taxa")
######################

#Pivot longer the data frame
relabund2 <- relabund2 %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "count") 
relabund2$Level <- "Phylum"
  #write.csv(relabund2, file = "test_long.csv")

#Load in the pivoted data
test.long <- read.csv("phylum_longData.csv", header = TRUE, sep = ",")

#Get the mean frequency of each taxa to order the phylum legend
mean_freqs <- test.long %>%
  group_by(Sample, Phylum) %>%
  summarise(phyl_count = sum(count),
            .groups = 'drop') %>%
  group_by(Sample) %>%
  summarise(phyl_freq = phyl_count / sum(phyl_count),
            Phylum = Phylum,
            .groups = 'drop') %>%
  group_by(Phylum) %>%
  summarise(mean = mean(phyl_freq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)

phylum_long <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Phylum_Level/Phylum.long.csv")
#Calculate alpha diversity metrics
richness <- function(x){
  sum(x > 0)
}

shannon <- function(x){
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}


shared <- phylum_long %>%
  group_by(Sample) %>%
  summarize(sobs = richness(count),
            Shannon = shannon(count),
            Simpson = simpson(count))
write.csv(shared, file = "alphaDiversity_Metrics_PhylumLevel.csv")

#Read in the edited alphaDiversity_metric_FamilyLevel.csv
alphaDiv <- read.csv("alphaDiversity_Metrics_PhylumLevel.csv", header = TRUE, sep = ",")

alphaDiv <- alphaDiv %>%
  group_by(Phenotype)

alphaDiv <- alphaDiv %>%
  mutate(Age = factor(Age, levels = Age_order))
alphaDiv <- alphaDiv %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

#See if the data is parametric
ShanHist <- hist(alphaDiv$Shannon)
SimpHist <- hist(alphaDiv$Simpson)
shapiro.test(alphaDiv$Shannon) #If P-value is > 0.05, we assume normality
shapiro.test(alphaDiv$Simpson)

#Plot distribution of Shannon Diversity
Shannons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + , nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  ylim(2.0,3.1) +
  labs(title = "Phylum level Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "bottom")
Shannons_boxplot
ggsave("Shannon_FamilyLevel_diversity.pdf", Shannons_boxplot, width = 12, height = 8)

#Plot distribution of Simpson Diversity
Simpsons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  labs(title = "Family level Simpson Diversity") +
  ylim(0,0.31) +
  theme_bw() +
  theme(legend.position = "bottom")
Simpsons_boxplot
ggsave("Simpsons_FamilyLevel_diversity.pdf", Simpsons_boxplot, width = 12, height = 8)


#Convert phyla to factors with a specific order based on the order obtained above (mean frequency)
Phyla_order <- c("Bacteroidetes",
                 "Proteobacteria",
                 "Actinobacteria",
                 "Firmicutes",
                 "Tenericutes",
                 "Verrucomicrobia",
                 "Unknown_phylum",
                 "Deferribacteres",
                 "Bacteria_unclassified",
                 "Cyanobacteria",
                 "Candidatus_Saccharibacteria")

#Re-order
test.long <- test.long %>%
  mutate(Phylum = factor(Phylum, levels = Phyla_order))

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")

#Re-order
test.long <- test.long %>%
  mutate(Age = factor(Age, levels = Age_order))

#Convert Sample to factor and set desired order
Sample_order <- c("AK1", "AK2", "AK3", "AK4", "AK5", "AK6", "AK7", "AK8", "AK9", "AK10", "AK11", "AK12",
                  "AW1", "AW2", "AW3", "AW4", "AW5", "AW6", "AW7", "AW8", "AW9", "AW10", "AW11", "AW12",
                  "BK1", "BK2", "BK3", "BK4", "BK5", "BK6", "BK7", "BK8", "BK9", "BK10", "BK11", "BK12",
                  "BW1", "BW2", "BW3", "BW4", "BW5", "BW6", "BW7", "BW8", "BW9", "BW10", "BW11", "BW12")

#Re-order sample name
test.long <- test.long %>%
  mutate(Sample = factor(Sample, levels = Sample_order))

#Convert Phenotype order to factor and set desired order
Phenotype_order <- c("WT", "KO")

#Re-order phenotype order
test.long <- test.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))


#Plot
ggplot(test.long, (aes(x = Sample, y = count))) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Phylum-level Relative Abundnace",
       subtitle = "A/J and C57BL/6 mice with and without Ptges KO")




#Read in the akkermansia data
akk <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Microbiome/Akkermansia.csv", sep = ",")

namesAkk <- akk$taxon

#sumsAkk <- rowSums(akk[,2:ncol(akk)])
#relAkk <- akk[,2:ncol(akk)]/sumsAkk

#rownames(akk) <- namesAkk
#akk <- rownames_to_column(akk, var = "taxon")

akk.long <- akk %>%
  pivot_longer(-taxon, names_to = "Sample", values_to = "count") 
write.csv(akk.long, file = "AkkermansiaTEST.long.csv")

#Read akk long
akk.long <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Akkermansia_muciniphila/Akkermansia.long.csv", header = TRUE, sep = ",")
  #This is the relative abundance one. For read counts one, see Akkermansia.long.csv

#Reorder factors
akk.long <- akk.long %>%
  mutate(Age = factor(Age, levels = Age_order))
akk.long <- akk.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

#Replace the underscore in species name with a space
akk.long$Species <- gsub("_", " ", akk.long$Species)

#Test for normality/parametric data
shapiro.test(akk.long$count)
  #If p-val is significant, we cannot assume normality is met. Need a non-parametric test like Mann-Whitney U
  


#Plot box and whisker faceted by age/strain
AkkBox <- ggplot(akk.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Akkermansia muciniphila"))) +
  theme_bw()
AkkBox <- AkkBox + theme(legend.position = "none")
AkkBox

#Plot box and whisker faceted by strain/age
ggplot(akk.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Akkermansia muciniphila."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#Split the data by 8 weeks and 20 weeks for individual graphs
akk8 <- akk.long[1:44,]
akk20 <- akk.long[45:nrow(akk.long),]

#8 Weeks
ggplot(akk8, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Akkermansia muciniphila."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#20 Weeks
ggplot(akk20, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Akkermansia muciniphila"), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()


# E. coli
coli <- read.csv("/users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/E.coli_Counts.csv", header = TRUE, sep = ",")
namesColi <- coli$taxon

coli.long <- coli %>%
  pivot_longer(-taxon, names_to = "Sample", values_to = "count") 
write.csv(coli.long, file = "Ecoli.long.csv")

#Read in Ecoli.long
coli.long <- read.csv("/users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Ecoli/Ecoli.long.csv", header = TRUE, sep = ",")

#Reorder factors
coli.long <- coli.long %>%
  mutate(Age = factor(Age, levels = Age_order))
coli.long <- coli.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

EcoliBox <- ggplot(coli.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Escherichia coli"))) +
  theme_bw()
EcoliBox <- EcoliBox + theme(legend.position = "none")


###########
###########
###########
###########
###########

#Prevotella species
prevo <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Microbiome/Prevotella.csv", sep = ",")

prevo.long <- prevo %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "count") 
write.csv(prevo.long, file = "Prevotella.long.csv")

#Read in Prevotella long format
prev.long <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Prevotella_sp/Prevotella.long.csv", header = TRUE, sep = ",")

#Reorder factors
prev.long <- prev.long %>%
  mutate(Age = factor(Age, levels = Age_order))
prev.long <- prev.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

#Plot box and whisker faceted by age/strain
PrevoBox <- ggplot(prev.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste(italic("Prevotella sp."), "Unclassified"))) +
  theme_bw()
PrevoBox <- PrevoBox + theme(legend.position = "none")



#Individual plots
#Split the data by 8 weeks and 20 weeks for individual graphs
prev8 <- prev.long[1:44,]
prev20 <- prev.long[45:nrow(akk.long),]

#8 Weeks
ggplot(prev8, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Prevotella sp."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#20 Weeks
ggplot(prev20, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Prevotella sp."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#########
#########
#########
#########
#########

#Read in the long file that contains both
both.long <- read.csv("bothLong.csv", header = TRUE, sep = ",")

#Reorder factors
both.long <- both.long %>%
  mutate(Age = factor(Age, levels = Age_order))
both.long <- both.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

#Plot box and whisker faceted by age/strain
ggplot(both.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = "Presence of key bacterial species",
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#Plot box and whisker faceted by strain/age
ggplot(both.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = "Presence of key bacterial species",
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#########
#########
#########
#########
#########

#Read in data for Bacteroides sp.
bact <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Microbiome/Bacteroides.csv", sep = ",")

bact.long <- bact %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "count") 
write.csv(bact.long, file = "Bacteroides.long.csv")

#Read in Prevotella long format
bact.long <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Bacteroides_sp/Bacteroides.long.csv", header = TRUE, sep = ",")

#Reorder factors
bact.long <- bact.long %>%
  mutate(Age = factor(Age, levels = Age_order))
bact.long <- bact.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))


#Plot box and whisker faceted by age/strain
BactBox <- ggplot(bact.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste(italic("Bacteroides sp."), "Unclassified"))) +
  theme_bw()
BactBox <- BactBox + theme(legend.position = "none")

#Plot box and whisker faceted by strain/age
ggplot(bact.long, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Bacteroides sp."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#Individual plots
#Split the data by 8 weeks and 20 weeks for individual graphs
bact8 <- bact.long[1:44,]
bact20 <- bact.long[45:nrow(akk.long),]

#8 Weeks
ggplot(bact8, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Bacteroides sp."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()

#20 Weeks
ggplot(bact20, aes(x = Phenotype, y = count, fill = Species)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "blue") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Age + Strain, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(paste("Presence of ", italic("Bacteroides sp."), "Unclassified")),
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO"))) +
  theme_bw()



setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/")
combinedInfo <- cowplot::plot_grid(AkkBox, PrevoBox, BactBox, EcoliBox)
ggsave("Akkermansia_muciniphila_Boxplot.pdf", AkkBox,  width = 12, height = 8)
ggsave("Prevotella_species_unclassified_Boxplot.pdf", PrevoBox, width = 12, height = 8)
ggsave("Bacteroides_species_unclassified_Boxplot.pdf", BactBox, width = 12, height = 8)
ggsave("Escherichia_coli_Boxplot.pdf", EcoliBox, width = 12, height = 8)



####Looking into proteobacteria
#Read in the proteobacteria_orderData.csv file
Proteo <- read.csv("Proteobacteria_orderData.csv", header = TRUE, sep = ",")
namesProteo <- Proteo$Order
Proteo$Order <- NULL

#Get relative abundance
total <- rowSums(Proteo[,2:ncol(Proteo)])
relativeAbund <- Proteo[,2:ncol(Proteo)]/total

rownames(relativeAbund) <- namesProteo
relativeAbund <- rownames_to_column(relativeAbund, var = "Order")

#Pivot longer the data frame
relativeAbund <- relativeAbund %>%
  pivot_longer(-Order, names_to = "Sample_ID", values_to = "count") 
write.csv(relativeAbund, file = "order.long.csv")

#Read in the order.long dataframe after editing by hand
ord.long <- read.csv("order.long.csv", header = TRUE, sep = ",")

#Reorder factors
ord.long <- ord.long %>%
  mutate(Age = factor(Age, levels = Age_order))
ord.long <- ord.long %>%
  mutate(Phenotype = factor(Phenotype, levels = Phenotype_order))

#Get the mean frequency of each taxa to order the Order legend
mean_freqs <- ord.long %>%
  group_by(Sample, Order) %>%
  summarise(order_count = sum(count),
            .groups = 'drop') %>%
  group_by(Sample) %>%
  summarise(order_freq = order_count / sum(order_count),
            Order = Order,
            .groups = 'drop') %>%
  group_by(Order) %>%
  summarise(mean = mean(order_freq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)

#Get the mean frequency of each taxa to order the Order legend
mean_freqs <- ord.long %>%
  group_by(Sample, Class) %>%
  summarise(class_count = sum(count),
            .groups = 'drop') %>%
  group_by(Sample) %>%
  summarise(class_freq = class_count / sum(class_count),
            Order = Class,
            .groups = 'drop') %>%
  group_by(Order) %>%
  summarise(mean = mean(class_freq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)

#Convert Class to factors with a specific order based on the order obtained above (mean frequency)
Class_order <- c("Alphaproteobacteria",
                 "Betaproteobacteria",
                 "Gammaproteobacteria",
                 "Deltaproteobacteria",
                 "Epsilonproteobacteria")


#Re-order
ord.long <- ord.long %>%
  mutate(Class = factor(Class, levels = Class_order))


#Convert Orders to factors with a specific order based on the order obtained above (mean frequency)
Order_order <- c("Desulfovibrionales",
                 "Campylobacterales",
                 "Rhizobiales",
                 "Burkholderiales",
                 "Caulobacterales",
                 "Kiloniellales",
                 "Rhodospirillales",
                 "Neisseriales",
                 "Kordiimonadales",
                 "Pseudomonadales",
                 "Enterobacterales",
                 "Rhodobacterales",
                 "Other")

#Re-order
ord.long <- ord.long %>%
  mutate(Order = factor(Order, levels = Order_order))
ord.long$Order[is.na(ord.long$Order)] <- "Other"

#Convert Sample to factor and set desired order
Sample_order <- c("AK1", "AK2", "AK3", "AK4", "AK5", "AK6", "AK7", "AK8", "AK9", "AK10", "AK11", "AK12",
                  "AW1", "AW2", "AW3", "AW4", "AW5", "AW6", "AW7", "AW8", "AW9", "AW10", "AW11", "AW12",
                  "BK1", "BK2", "BK3", "BK4", "BK5", "BK6", "BK7", "BK8", "BK9", "BK10", "BK11", "BK12",
                  "BW1", "BW2", "BW3", "BW4", "BW5", "BW6", "BW7", "BW8", "BW9", "BW10", "BW11", "BW12")

#Re-order sample name
ord.long <- ord.long %>%
  mutate(Sample = factor(Sample, levels = Sample_order))

#Viridis
Ordern <- 13
Classn <- 5
Ordercolors <- c(microshades_palette("micro_blue", 5),
            microshades_palette("micro_purple",5),
            microshades_palette("micro_cvd_turquoise", 5))
Classcolors <- c(microshades_palette("micro_blue", 1),
                 microshades_palette("micro_purple",2),
                 microshades_palette("micro_cvd_turquoise", 2),
                 microshades_palette("micro_cvd_purple", 2))
                 

colors <- viridis(n = n, option = "C")

#For order level abundance
colors <- c("#4773aa", "#e4888b", "#7cb287", "#fbdb88", "#9355dc",
                     "#5fc7e3", "#fa9441", "#22723d", "#722257", "#e0699f",
                     "#ffdab3", "#a3b899", "#d4af37", "#5b6988", "#ddb9b3")
                     



#Plot for class-level abundance
ggplot(ord.long, (aes(x = Sample, y = count))) +
  geom_bar(aes(fill = Class), stat = "identity", position = "fill", width = 1) +
  scale_fill_manual(values = c("#4773aa", "#e4888b", "#7cb287", "#fbdb88", "#9355dc")) +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  facet_nested_wrap(~ Strain + Phenotype + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Proteobacteria: Class-level Relative Abundnace",
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO")))

#Plot for Order-level abundance
ggplot(ord.long, (aes(x = Sample, y = count))) +
  geom_bar(aes(fill = Order), stat = "identity", position = "fill", width = 1) +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  facet_nested_wrap(~ Strain + Phenotype + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Proteobacteria: Order-level Relative Abundnace",
       subtitle = expression(paste("A/J and C57BL/6 mice WT and", italic(" Ptges "), "KO")))
