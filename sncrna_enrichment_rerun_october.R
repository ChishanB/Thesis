# Purpose: Redoing the enrichment analysis for aim 2 (mirna target prediction) because previously, 
# I used multimiR query results which hadn't been filtered to exclude negative results from tarbase. 
# I thought querying for 'validated' interactions was enough. But it wanst! Yay!

# Directly ripping pieces of script from sncrna_june to try and keep pipeline as consisent as possible. 
# Need to do it with aim 3 too! YAy!

# Author: Chishan Burch no. 1 TWICE fan 
# Date: 30.10.25

###################################################################

# Load in libraries
library(org.Mm.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(openxlsx)
library(readxl)
library(edgeR)
library(clusterProfiler)

# Load mouse database
library(org.Mm.eg.db)  
organism <- org.Mm.eg.db    


# The file I'm loading has mature validated interactions between DE miRNAs and
# genes in the SV whether or not they are DE

miRNA <- read.csv("./3_results/oct_multimiR_all_targets_of_DE_mature_miRNA.csv")

# Includes all genes passing basic QC
genes <- read.csv("./3_results/julyrerun_DESeq_genesandTEs.csv")

# Remove TEs
genes <- genes %>%
  filter(grepl("ensmus", X, ignore.case = TRUE))

# Filter for DE genes so that I can do enrichment on DE targets only
DE_genes <- genes %>%
  filter(!is.na(padj) & padj <= 0.05)

# Set up background list for aim 2 (All SV genes which are targets of DE miRNAs)
background_list <- genes$X

# Convert background list to Entrez ids for KEGG. GO can use ensembl

# Convert gene names to KEGG IDs
KEGG_background_list <- bitr(
  background_list,
  fromType = "ENSEMBL",  # Source ID type (e.g., HGNC symbol)
  toType = "ENTREZID",      # Target ID type
  OrgDb = organism  
)

KEGG_background_list <- KEGG_background_list$ENTREZID

# To get miRNA hit list:

#Unique multimir targets
unique_target_ids <- miRNA$target_entrez

# Check if there is some overlap between the hit list and background list
length(intersect(unique_target_ids, KEGG_background_list)) 
# 496 overlaps. Meaning, 496 genes in the SV (DE or not) are targeted by DE miRNAs

KEGG_enrichment_results <- enrichKEGG(
  gene = unique_target_ids,
  universe = KEGG_background_list,
  organism = "mmu",  # Mouse
  keyType = "kegg",  # Key type in KEGG
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(KEGG_enrichment_results)

KEGG_dotplot <- dotplot(KEGG_enrichment_results, color = "p.adjust")

# Save KEGG dotplot
#png("2_figures/oct_KEGG_dotplot_aim2.png", width = 1200, height = 900, res = 150)
#dotplot(KEGG_enrichment_results, color = "p.adjust")
#dev.off()

KEGG_results <- as.data.frame(KEGG_enrichment_results@result)

#write.csv(KEGG_results, file = "./3_results/sncRNA_Oct_KEGG.csv")


# Now for GO-ORA

# To get miRNA hit list:

#Unique multimir targets
unique_target_ids <- miRNA$target_ensembl

go_enrich <- enrichGO(gene = unique_target_ids,
                      universe = background_list,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


go_results <- go_enrich@result
go_results$GeneRatio <- gsub("/","//",go_results$GeneRatio)
go_results$BgRatio <- gsub("/","//",go_results$BgRatio)


#write.xlsx(go_results, file = "./3_results/sncRNA_Oct_GOenrich.csv")

## Dot plot
enrichplot::dotplot(go_enrich)

#ggsave(file = "./2_figures/oct_GOenrich_dotplot_aim2.png")

# Save GO dotplot
#png("2_figures/oct_GOenrich_dotplot_aim2.png", width = 1200, height = 900, res = 150)
#dotplot(go_enrich, color = "p.adjust")
#dev.off()

# I'm gonna do the same for aim 3 in a separate script
