# Purpose: A total of 16 461 genes detected in the SV were known (validated)
# targets of miRNAs which were DE in response to acrylamide exposure. 
# A total of 185 genes were DE in response to acrylamide exposure. 
# 100 of these DE genes were targets of DE miRNAs.

# You can see how I got these results in the script KEGG_GO_of_DE_targets_only.R

# This script is just for a venn diagram

# Date: 22.9.25

#################################

library(ggvenn)

# Create dummy gene sets
# Large set: 16,461 miRNA targets
setA <- paste0("GeneA", 1:16461)

# DE genes: 185 total but i had to write it as 285 to make the correct number display. 
setB <- paste0("GeneB", 1:185)

# Force 100 overlap
setB[1:100] <- setA[1:100]   # first 100 DE genes overlap with miRNA targets

# Put sets into a list
gene_sets <- list(
  "DE miRNA targets" = setA,
  "DE genes" = setB
)

# Draw venn diagram
venn_plot <- ggvenn(
  gene_sets,
  fill_color = c("skyblue", "salmon"),
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 5
)

ggsave(
  filename = "./2_figures/venn_DE_miRNA_targets_DE_genes.png",
  plot = venn_plot,
  width = 8,    
  height = 6,   
  dpi = 300
)

