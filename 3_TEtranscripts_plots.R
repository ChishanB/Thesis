
# Purpose: To improve readability of the TE-analysis pipeline, 
# clean up techniques from previous scripts and make the pipeline more readily reproducible. 
# (At the moment the Rmarkdown report is only legible to me)

##################################

# Load packages
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

# Load raw data
TEs_and_genes <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")
#Rename column because it's mouse data, not human. 
 colnames(TEs_and_genes)[colnames(TEs_and_genes) == "hgnc.TE"] <- "mouse.TE"
#write.csv(TEs_and_genes, file = "./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.csv", row.names = FALSE)

# Separate TEs and genes for easier use
TE_subset <- base::subset(TEs_and_genes, grepl('TE', Type, ignore.case = TRUE))
gene_subset <- base::subset(TEs_and_genes, grepl('gene', Type, ignore.case = TRUE))

# Only 1 TE was considered statistically significant after DESeq analysis and 
# padj = 0.05 (RLTR6-int), hence it was decided it wouldn't be productive to 
# compare relative presence of different TE classes and families in the data.

# Gene class information wasn't available, so it was decided it might be worth 
# obtaining a list of sex-linked genes that may have been statistically significant. 

# List of sex-linked genes if needed (unfiltered for padj)
#X_linked <- gene_subset[gene_subset$Chromosome == "X", ]
#Y_linked <- gene_subset[gene_subset$Chromosome == "Y", ]

# For volcano plots
all_DESeq_results <- read.csv("./3_results/allgenes_DESeq.csv")
TE_DESeq_results <- read.csv("./3_results/TE_only_DESeq.csv")



# Add a column labelling whether TEs are upregulated or downregulated
all_DESeq_results <- all_DESeq_results %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange < -0.585, "Sign_Down",
                           if_else(padj < 0.05 & log2FoldChange > 0.585, "Sign_Up", "NS"))) 

TE_DESeq_results <- TE_DESeq_results %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange < -0.585, "Sign_Down",
                           if_else(padj < 0.05 & log2FoldChange > 0.585, "Sign_Up", "NS"))) 

# Create a -log10 column for y-axis of volcano plots
all_DESeq_results <- all_DESeq_results %>% mutate(log_padj = -log10(padj))
TE_DESeq_results <- TE_DESeq_results %>% mutate(log_padj = -log10(padj))


# Declare colours
colours_1 <- c("Sign_Up" = "#AA3377", "Sign_Down" = "#228833", "NS" = "#BBBBBB") 
opacity <- c("Sign_Up" = 1, "Sign_Down" = 1, "NS" = 0.5)

# Declare padj threshold
p_adj_threshold <- 0.05
options(ggrepel.max.overlaps = Inf)

# Plot
TE_DESeq_results %>%
  filter(is.na(sign_DE) != TRUE) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sign_DE, alpha = sign_DE, fill = sign_DE)) +
  geom_point(size = 5, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + theme_bw() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 


all_DESeq_results %>%
  filter(is.na(sign_DE) != TRUE) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sign_DE, alpha = sign_DE, fill = sign_DE)) +
  geom_point(size = 5, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + theme_bw()  +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 



# need to do an inner join. all_DESeq_results contains the results but not the gene classes. we will inner 
# join this with the og (raw) data table so we can take the gene classes from it and join them with the results table
# while omitting other data from the raw data table

#rename id column in results table to match raw data
names(all_DESeq_results)[names(all_DESeq_results) == "X"] <- "ensembl.TE"


joined_data <- inner_join(all_DESeq_results, TEs_and_genes, by = "ensembl.TE")

# remove columns we dont need again
joined_data <- joined_data %>% select(-11, -12, -13, -14, -17, -18, -19, -20, -22, -23, -24, -25, -26, -27, -28, -29)

# remove rows which are NA in columns resulting from DE analysis to isolate our filtered data
# chose to omit rows with "NA" padj, because this indicates padj could not compute a meaniingful padj for those entries. 
joined_data_clean <- joined_data %>% filter(!is.na(padj))

# save that
#write.csv(joined_data_clean, file = "./3_results/clean_joined_data.csv", row.names = FALSE)

# isolate the gene classes and then make a pie chart of them!


### Pie charts time #####################################################################

#load clean joined data created in 3_TEtranscripts_plots.R

#joined_data_clean <- read.csv("./3_results/clean_joined_data.csv")

# Obtain no. upregulated and downregulated
nrow(joined_data_clean[joined_data_clean$log2FoldChange > 0.585, ])
#51 upregulated
nrow(joined_data_clean[joined_data_clean$log2FoldChange < -0.585, ])
#208 downregulated

UPREG <- joined_data_clean %>% filter(joined_data_clean$log2FoldChange > 0.585)
DOWNREG <- joined_data_clean %>% filter(joined_data_clean$log2FoldChange < -0.585)


# Prepare data for pie charts (remove the ones which are NA for gene class)
UPREG <- UPREG %>% filter(!is.na(gene.class))
UPREG <- table(UPREG$gene.class)

DOWNREG <- DOWNREG %>% filter(!is.na(gene.class))
DOWNREG <- table(DOWNREG$gene.class)


upreg_types <- data.frame(Type = names(UPREG), Count = as.numeric(UPREG))
downreg_types <- data.frame(Type = names(DOWNREG), Count = as.numeric(DOWNREG))



# Plot upregulated 
## Declare custom colours
custom_colors <- c(protein_coding = "#77AADD", lincRNA = "#EE8866", antisense = "#FFAABB", miRNA = "#EEDD88", 
                   processed_pseudogene = "#99DDFF", processed_transcript = "#44BB99", 
                   transcribed_unprocessed_pseudogene = "#DDDDDD", processed_pseudogene = "#AAAA00", IG_C_gene = "#BBCC33", 
                   unprocessed_pseudogene = "#545454")





# Add a column for legend labels
upreg_types$legend_labels <- paste(upreg_types$Type, "(", upreg_types$Count, ")", sep = "")

up_plot <- ggplot(upreg_types, aes(x = " ", y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 1) +
  coord_polar(theta = "y") +
  labs(title = "Upregulated gene classes", fill = "Type") +
  theme_void() + 
  scale_fill_manual(values = custom_colors, labels = upreg_types$legend_labels) + 
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) 


# Add a column for legend labels
downreg_types$legend_labels <- paste(downreg_types$Type, "(", downreg_types$Count, ")", sep = "")

down_plot <- ggplot(downreg_types, aes(x = " ", y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 1) +
  coord_polar(theta = "y") +
  labs(title = "Downregulated gene classes", fill = "Type") +
  theme_void() +
  scale_fill_manual(values = custom_colors, labels = downreg_types$legend_labels) +
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#ggsave(plot = up_plot, filename = "./3_results/up_geneclasses.png", width = 6, height = 7, dpi = 300)
#ggsave(plot = down_plot, filename = "./3_results/down_geneclasses.png", width = 6, height = 7, dpi = 300)
  

# heatmap time

# redo the deseq
# vst normalise


# Load packages
library(pheatmap)
library(DESeq2)


# Load raw data
TEs_and_genes <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")
sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")
# loaf joined_data_clean if not done already

# Prepare raw data and metadata
sample_metadata <- as.data.frame(sample_metadata) %>% column_to_rownames(., var = "SampleID")


# Remove the columns from a which are not counts or name info
a <- TEs_and_genes %>% select(-1, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)


# Check for duplicates in the column we want to make the row names
#duplicated_values <- a$hgnc.TE[duplicated(a$hgnc.TE)]
#print(duplicated_values)

# Make row names unique (e.g., pakap triplicates are now pakap, pakap.1, pakap.2)
a <- a %>%
  mutate(rowname = make.unique(as.character(hgnc.TE))) %>%
  column_to_rownames(var = "rowname")

# Lets filter for TEs with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

# Filtering step
a_filtered <- a %>%
  filter(rowSums(. >= min.counts) >= min.samples)
# Make DESeqDataSetFromMatrix(countData = a_filtered, colData = sample_metadata,  : ncol(countData) == nrow(colData) true
# i.e. make the number of columns the same as the other data frame
a_filtered <- a_filtered %>% select(-1)

# Re create DESeq results
dds <- DESeqDataSetFromMatrix(countData = a_filtered, # Needs to be raw counts with rownames set as RNA ids
                              colData = sample_metadata, # Data frame of meta data
                              desig = ~ Treatment_group) 
dds <- DESeq(dds)


# vst normalise the DEseq object
vsd <- vst(dds, blind = T)

# Convert it to a data frame
str(assay(vsd))
Deseq_vst <- assay(vsd)


# Subset for names of DE genes (we know this data is clean)
gene_names_vector <- joined_data_clean$hgnc.TE

vst_subset <- Deseq_vst[rownames(Deseq_vst) %in% gene_names_vector, ]


# Create Mean centred and z-score (to create two different plots)
##Mean_centred <- vst_subset- rowMeans(vst_subset)

#Z_scored <- t(scale(t(vst_subset)))
####Z_scored <- t(scale(t(-vst_subset))) # 20/11/23 changed to negative to reflect correct direction



# Use pheatmap to plot, using default hierarchical clustering and using the following variables:

## WHAT AM I MISSING HERE??


sample_annotation <- sample_metadata 


pheatmap(vst_subset, scale = "none", clustering_distance_rows = "euclidean", 
                       clustering_distance_cols = "euclidean", annotation_col = sample_annotation)



#UP TO THIS PARTT

# Create function to save heatmap
save_pheatmap_png <- function(x, filename, width=980, height=1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Save heatmap
#save_pheatmap_png(pheatmap_1, "./2_figures/pheatmap_1.png")

# Select the types from gene_metadata, and assign ids as rownames
gene_metadata_df <- gene_metadata %>%
  filter(sRNA.ID_type %in% rownames(vst_subset)) %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type") %>%
  dplyr::select(Type)

sample_annotation <- sample_metadata 
miRNA_annotation <- gene_metadata_df


# Mean-centred pheatmap
pheatmap_mean_centred <- pheatmap(Mean_centred, scale = "none", clustering_distance_rows = "euclidean", 
                                  clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
                                  annotation_row = miRNA_annotation, annotation_colors = custom_colours)
# Save
save_pheatmap_png(pheatmap_mean_centred, "./2_figures/pheatmap_mean_centred.png")

# Z-scored pheatmap
pheatmap_zscored <- pheatmap(Z_scored, scale = "none", clustering_distance_rows = "euclidean", 
                             clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
                             annotation_row = miRNA_annotation, annotation_colors = custom_colours,
                             annotation_names_row = FALSE, angle_col = 90, annotation_names_col = FALSE)

pheatmap_zscored2 <- pheatmap(Z_scored, scale = "none", clustering_distance_rows = "euclidean", 
                              clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
                              annotation_row = miRNA_annotation, annotation_colors = custom_colours,
                              annotation_names_row = FALSE, show_rownames = FALSE, angle_col = 45,
                              annotation_names_col = FALSE)


# Save
save_pheatmap_png(pheatmap_zscored, "./2_figures/pheatmap_zscored.png")
save_pheatmap_png(pheatmap_zscored2, "./2_figures/pheatmap_zscored2.png")

# TPM values pheatmap (ensure the TPM filtered counts csv was used in Deseq etc. 
#Other 3 plots use raw filtered counts csv)

vst_subset <- log2(vst_subset)

pheatmap_TPMcounts <- pheatmap(vst_subset, scale = "none", clustering_distance_rows = "euclidean", 
                               clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
                               annotation_row = miRNA_annotation, annotation_colors = custom_colours)
# Save
#save_pheatmap_png(pheatmap_TPMcounts, "./2_figures/pheatmap_TPMcounts_log2.png")
