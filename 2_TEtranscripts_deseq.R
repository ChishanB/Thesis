####################################################################################

#date: 25/6/24
#purpose: To import TEtranscripts raw data, conduct differential analysis, 
#construct piecharts showing what TE families are most represented
#email: chishanburch@gmail.com

#packages########

# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(DESeq2)

# Read the data
a <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")
#b <- read.delim("./1_data/150bp/raw_counts/TE_Class_only_counts.txt")
#c <- read.delim("./1_data/150bp/raw_counts/TE_Family_Class_only_counts.txt")
#d <- read.delim("./1_data/150bp/raw_counts/TE_only_and_all_counts.txt")
sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")


# Set treatment group levels
  as.data.frame(sample_metadata) %>%
  mutate(Treatment_group = factor(Treatment_group, levels = c("Control", "Acr")))


# Filter for TEs only, those with at least counts of 5 across 3 samples
# Also do one containing TEs and the rest f the genes
TE_subset <- subset(a, grepl('TE', Type, ignore.case = TRUE))
all_genes <- a

# Remove the columns from a which are not counts or name info
TE_subset <- TE_subset %>% select(-1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)
all_genes <- all_genes %>% select(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)


# Lets filter for TEs with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

rownames(TE_subset) <- NULL
rownames(all_genes) <- NULL

# Set the columns as the row names
TE_subset <- 
  TE_subset %>%
  data.frame() %>%
  column_to_rownames(var = "Name")

all_genes <- 
  all_genes %>%
  data.frame() %>%
  column_to_rownames(var = "ensembl.TE")

# Filtering step
TE_subset_filtered <- TE_subset %>%
  filter(rowSums(. >= min.counts) >= min.samples)

all_genes <- all_genes %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Remove rows containing NA values
if (any(is.na(TE_subset_filtered))) {
  TE_subset_filtered <- TE_subset_filtered[complete.cases(TE_subset_filtered), ]
}

if (any(is.na(all_genes))) {
  all_genes <- all_genes[complete.cases(all_genes), ]
}

# Create gene_metadata
gene_metadata <- 

# Perform DESeq analysis WITHOUT FORGETTING dds <- DESeq(dds) STEP

dds <- DESeqDataSetFromMatrix(countData = TE_subset_filtered,
                              colData = sample_metadata,
                              desig = ~ Treatment_group) #DO NOT FORGET THE ~ *********

dds2 <- DESeqDataSetFromMatrix(countData = all_genes,
                              colData = sample_metadata,
                              desig = ~ Treatment_group)


dds <- DESeq(dds)

dds2 <- DESeq(dds2)

#res <- results(dds) 

res <- results(dds, contrast = c("Treatment_group","Control","Acr")) 

levels(sample_metadata$Treatment_group)


res2 <- results(dds2, contrast = c("Treatment_group","Control","Acr")) 

levels(sample_metadata$Treatment_group)



res_df <- as.data.frame(res) # Coerce into a more easily accessible format.
res_df2 <- as.data.frame(res2)

#Save results

#write.csv(res_df, file = "./3_results/TE_only_DESeq.csv", row.names = TRUE)
#write.csv(res_df2, file = "./3_results/allgenes_DESeq.csv", row.names = TRUE)



###### These don't work, so we'll try using a different way
vsd <- vst(dds, blind = T) #TEs only
vst_matrix <- assay(vsd)


vsd2 <- vst(dds, blind = T) #All genes
vst_matrix2 <- assay(vsd)
######

# Perform variance standardisation with blind = T so the computer is not aware 
# of which treatment group the samples belong to (i.e., unsupervised analysis).

vsd <- varianceStabilizingTransformation(dds, blind = T) #it appears varianceStabilizingTransformation
# is a slightly different operation to vst
vst_matrix <- assay(vsd)
vsd2 <- varianceStabilizingTransformation(dds2, blind = T) #it appears varianceStabilizingTransformation
# is a slightly different operation to vst
vst_matrix2 <- assay(vsd)


# Save the variance stabilized DESeq outputs
#write.csv(vst_matrix, file = "./3_results/TE_only_vst_DESeq.csv", row.names = TRUE)
#write.csv(vst_matrix2, file = "./3_results/TE_only_vst_DESeq.csv", row.names = TRUE)

# DESeq2 PCA plot
colours <- c("Control" = "#745E96", "Acr" = "#36013F")


TE_DESeq_output <- 
  plotPCA(vsd, intgroup = "Treatment_group") +
  scale_color_manual(values = colours) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.8) + theme_minimal() +
  geom_label_repel(aes(label=colnames(vsd)), label.size = 0.5, box.padding = 0.25, label.padding = 0.5) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

allgenes_DESeq_output <- 
  plotPCA(vsd2, intgroup = "Treatment_group") +
  scale_color_manual(values = colours) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.8) + theme_minimal() +
  geom_label_repel(aes(label=colnames(vsd)), label.size = 0.5, box.padding = 0.25, label.padding = 0.5) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))


base::plot(allgenes_DESeq_output)
#Filter for padj <= 0.05

res_df_significant <- res_df %>%
  filter(padj <= 0.05) # This is for TEs

res_df_significant2 <- res_df2 %>%
  filter(padj <= 0.05) # This is for all genes meeting the thresholds


