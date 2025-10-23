# Purpose: Editing to use genes + TE Deseq results instead of just the genes results

# Date: 23.10.25


# Author: Chishan Burch


#############################################
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(data.table)
library(clusterProfiler)

## Data preparation steps 1-2
# get significantly DE (padj < 0.05) genes
DE_results <- fread("./3_results/april_DESeq_genesandTEs.csv")


# Remove TEs
DE_results <- DE_results[grepl("ENSMUS", DE_results$ensembl.TE), ]

# Apply basic QC
DE_results <- DE_results[!is.na(DE_results$padj) & DE_results$padj < 0.05, ]

# Remove extra column 
DE_results <- DE_results[, !"V1"]


# get validated miRNAs from previous analyses

miRNA_targets_up <- fread(file = "./3_results/multimiR_mature_mirnas_up_validated.csv")
miRNA_targets_down <- fread(file = "./3_results/multimiR_mature_mirnas_down_validated.csv")

# remove NAs
if (any(is.na(miRNA_targets_up))) {
  miRNA_targets_up <- miRNA_targets_up[complete.cases(miRNA_targets_up), ]
}

if (any(is.na(miRNA_targets_down))) {
  miRNA_targets_down <- miRNA_targets_down[complete.cases(miRNA_targets_down), ]
}

# combine upregulated and downregulated target dataframes
miRNA_targets <- rbind(miRNA_targets_up,miRNA_targets_down)

# rename column to match other dataframe so I can use a join

DE_results <- as.data.frame(DE_results)  # convert to data frame

DE_results <- DE_results %>%
  dplyr::rename(target_ensembl = ensembl.TE)
head(DE_results)

# this shows that out of 185 DE genes in the dataset, just 100 of them were validated miRNA targets
DE_validated <- semi_join(miRNA_targets, DE_results, by = "target_ensembl") %>%
  distinct(target_ensembl, .keep_all = TRUE)

# Save this 
#write.csv(DE_validated, file = "./3_results/oct_multimiR_DE_targets_of_DE_mature_miRNA.csv")


## KEGG and GO-ORA (step 3)

# Set up background list
genes <- fread(file = "./1_data/Transcriptomics/TE_and_gene_counts_Updated_and_Raw.txt")

# Isolating the gene names instead of the ensembl ids.
background_list <- genes$ensembl.TE

# BACKGROUND LIST 
head(background_list)

# Gotta convert background list to ensembl
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)  
organism <- org.Mm.eg.db    
#library(clusterProfiler)

# Convert gene names to KEGG IDs
KEGG_background_list <- bitr(
  background_list,
  fromType = "ENSEMBL",  # Source ID type (e.g., HGNC symbol)
  toType = "ENTREZID",      # Target ID type
  OrgDb = organism  
)

KEGG_background_list <- KEGG_background_list$ENTREZID

# set up hit list
unique_multimiR_ids <- unique(combined_df$target_entrez)

hit_list <- unique(DE_validated$target_entrez)


# KEGG analysis
KEGG_enrichment_results <- enrichKEGG(
  gene = hit_list,
  universe = KEGG_background_list,
  organism = "mmu",  # Mouse
  keyType = "kegg",  # Key type in KEGG
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Save KEGG results

KEGG_results <- as.data.frame(KEGG_enrichment_results@result)

#write.csv(KEGG_results, file = "./3_results/DE_validated_KEGG_september.csv")


#KEGG_barplot <- barplot(KEGG_enrichment_results, showCategory = 20)
KEGG_dotplot <- dotplot(KEGG_enrichment_results, color = "pvalue")

# Save KEGG barplot
#png("KEGG_barplot.png", width = 1200, height = 900, res = 150)
#barplot(KEGG_enrichment_results, showCategory = 20)
#dev.off()

# Save KEGG dotplot
png("./2_figures/KEGG_validated_DE_targets_dotplot.png", width = 1200, height = 900, res = 150)
dotplot(KEGG_enrichment_results, color = "pvalue")
dev.off()



#KEGG_results <- as.data.frame(KEGG_enrichment_results@result)

#write.csv(KEGG_results, file = "./3_results/sncRNA_July_KEGG.csv")

## GO-ORA

# Set up hit list and background list for GO-ORA, which takes ensembl ids
#genes <- fread(file = "./1_data/Transcriptomics/TE_and_gene_counts_Updated_and_Raw.txt")
#background_list <- genes$ensembl.TE

GO_background_list <- genes$ensembl.TE

# set up hit list
GO_hit_list <- unique(DE_validated$target_ensembl)



go_enrich <- enrichGO(gene = GO_hit_list,
                      universe = GO_background_list,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


# Extract and save results table

go_results <- go_enrich@result
go_results$GeneRatio <- gsub("/","//",go_results$GeneRatio)
go_results$BgRatio <- gsub("/","//",go_results$BgRatio)

# Save table
#write.csv(go_results, file = "./3_results/DE_validated_GO_ORA_september.csv")



## Dot plot
enrichplot::dotplot(go_enrich)

#ggsave(
#  filename = "./2_figures/GO_ORA_validated_DE_targets_dotplot.png",
#  plot = last_plot(),   # or your plot object
#  width = 10,           # width in inches
#  height = 8,           # height in inches
#  dpi = 300             # resolution (dots per inch)
#)


######## Adding in a volcano plot
library(ggplot2)

DE_results <- fread("./3_results/april_DESeq_genes.csv")

DE_results <- DE_results %>%
  rename(target_ensembl = V1)

DE_results <- DE_results %>%
  dplyr::mutate(sign_DE = if_else(log2FoldChange > 0, 
                                  "Upregulated", 
                                  "Downregulated"))

DE_results <- DE_results %>% 
  dplyr::mutate(log_padj = -log10(padj))


x <- semi_join(DE_results, DE_validated, by = "target_ensembl") %>%
  distinct(target_ensembl, .keep_all = TRUE)


## need to ad log_padj column

# Ensures the points which have padj < 0.05 are shown in the front for aesthetic reasons
x <- x %>%
  arrange(desc(padj))

# Volcano Plot
colours_1 <- c("Upregulated" = "#994455", "Downregulated" = "#0077BB") 
opacity <- c("Sign_Up" = 1, "Sign_Down" = 1)

plot_1 <- x %>%
  filter(!is.na(sign_DE)) %>%
  ggplot(aes(x = log2FoldChange, y = log_padj, col = sign_DE, alpha = sign_DE, fill = sign_DE)) +
  geom_point(size = 5, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + 
  labs(
    title = 'Volcano Plot of Differential Expression (DE targets of DE miRNAs)',
    x = "Log2 Fold Change",           # Change this to your desired x-axis label
    y = "-Log10 FDR",    # Change this to your desired y-axis label
    color = "Effect in treated group",           # Optional: legend title for color
    fill = "Effect in treated group",            # Optional: legend title for fill
    alpha = "Effect in treated group"            # Optional: legend title for alpha
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

print(plot_1)
#ggsave("./2_figures/DE_targets_of_DE_miRNAs.png", plot = plot_1, width = 20, height = 15, units = "cm", dpi = 300)
