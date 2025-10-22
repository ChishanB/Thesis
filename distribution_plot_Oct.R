# Date: 28/7/25
# Author: Chishan Burch
# Purpose: 

# Literally the same as distribution_plot_July but I'm including hairpin miRNAs now

#
# (Instructions from John): Single omics analysis – link mature miRNA to mRNA targets 

#Generate a dataframe which shows us the common miRNA targets. What this means is that if we had 10 mature miRNA that were DE in the dataset, how many of these target the same gene. What I envision is that we would have a list of genes from multimiR data, and for each of those genes, you would have all the miRNA from the DE miRNA that target it – ie Gene 1 might be targeted by miR-223, miR-155 and miR-146a, while gene 2 might only be targeted by miR-223
#Plot this as a distribution plot showing the number of genes on y axis, and number of miRNA that target each gene on x. So you might have 300 genes that are targeted by only one miRNA, but 30 genes targeted by 2 and so on
#The final step will be to link the DE miRNA with appropriately regulated DE genes. So remember, if a gene is a real target of a down-regulated miRNA, you would expect it to be up-regulated and what we are hoping you will do is plot on a volcano plot all gene targets of the down-regulated miRNA. So the volcano plot will only show those gene targets and it will see how many of those are up-regulated. You will then repeat this with the up-regulated miRNA and down-regulated targets. 


library(openxlsx)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(openxlsx)
library(readxl)
library(ggrepel)
library(edgeR)
library(data.table)


miRNAs <- fread(file = "./3_results/july_DESeq_sncRNAs.csv")


miRNAs <- miRNAs %>%
  filter(sign_DE != "NS")

# filter for mature mirnas and hairpins and then rbind()
mature <- miRNAs %>%
  filter(grepl("mature", V1))

mature$V1 <- gsub(" mature$", "", mature$V1)


hairpins <- miRNAs %>%
  filter(grepl("hairpin", V1))

hairpins$V1 <- gsub(" hairpin$", "", hairpins$V1)

miRNAs <- rbind(mature,hairpins)

# miRNAs contains all differentially expressed miRNAs. 
# We want to know how many target the same genes.

mature_up <- fread(file = "./3_results/multimiR_mature_mirnas_up_validated.csv")
mature_down <- fread(file = "./3_results/multimiR_mature_mirnas_down_validated.csv")
hairpins_up <- fread(file = "./3_results/multimiR_hairpin_mirnas_up_validated.csv")
hairpins_down <- fread(file = "./3_results/multimiR_hairpin_mirnas_down_validated.csv")


miRNA_targets <- rbind(mature_up,mature_down,hairpins_up,hairpins_down)

# Seeing what the plots come out like if i only use hairpins
#miRNA_targets <- rbind(hairpins_up,hairpins_down)


 intersect(miRNAs$V1, miRNA_targets$mature_mirna_id)

 ### This is including only mature miRNAs
 #[1] "mmu-miR-128-3p"  "mmu-miR-145a-5p" "mmu-miR-199a-3p" "mmu-miR-199b-3p" 
 #"mmu-miR-19b-3p"  "mmu-miR-211-5p"  "mmu-miR-221-3p" 
#[8] "mmu-miR-410-3p"  "mmu-miR-532-5p" 
 
# This is when we include both mature and hairpins
 #>  intersect(miRNAs$V1, miRNA_targets$mature_mirna_id)
 #[1] "mmu-miR-199a-3p" "mmu-miR-199b-3p" "mmu-miR-532-5p"  "mmu-miR-128-3p" 
 #[5] "mmu-miR-211-5p"  "mmu-miR-145a-5p" "mmu-miR-19b-3p"  "mmu-miR-221-3p" 
 #[9] "mmu-miR-410-3p"  "mmu-miR-1a-3p"   "mmu-miR-3475-3p" "mmu-let-7b-3p"  
 #[13] "mmu-miR-23a-3p"  "mmu-miR-1b-5p"   "mmu-miR-504-5p" 
 
 
# So 8 of the DE mature miRNAs had targets in multimiR and only 1 hairpin had a 
 # DE target in multimiR
 

common_ids <- intersect(miRNAs$V1, miRNA_targets$mature_mirna_id)

filtered_targets <- miRNA_targets %>%
  filter(mature_mirna_id %in% common_ids)

nrow(filtered_targets)
miRNA_target_numbers <- table(filtered_targets$mature_mirna_id)

#write.csv(miRNA_target_numbers, file = "./3_results/oct_DE_matureandhairpin_miRNA_hits_freq.csv")


# Summarize: for each gene, list all DE miRNAs that target it
gene_miRNA_df <- filtered_targets %>%
  group_by(target_symbol) %>%
  summarise(
    miRNAs = paste(unique(mature_mirna_id), collapse = ", "),
    n_miRNAs = n_distinct(mature_mirna_id)
  ) %>%
  arrange(desc(n_miRNAs))  # optional: sort by number of miRNAs


#write.csv(gene_miRNA_df, file = "./3_results/miRNA_targets_hairpinandmature_october.csv")


barplot <- ggplot(gene_miRNA_df, aes(x = n_miRNAs)) +
  geom_bar(fill = "#cc5500") +
  labs(
    title = "Number of DE miRNAs Targeting Each Gene detected in the mouse SV",
    x = "Number of DE miRNAs (hairpins and mature) targeting the same gene",
    y = "Number of genes"
  ) +
  theme_minimal()

plot(barplot)
#ggsave(
#  filename = "./2_figures/miRNA_distribution_plot_hairpinsandmature.png", # change to whatever you like
 #plot     = barplot,
#  width    = 8,    # in inches
#height   = 6,    # in inches
# dpi      = 300   # resolution
#)


# For each group of DE miRNAs (upregulated or downregulated), look at their gene 
# targets, then check whether those genes are appropriately regulated in your DE mRNA dataset.

DE_mRNAs <- fread("./3_results/april_DESeq_genes.csv")


miRNAs <- miRNAs %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange < 0, "Downregulated",
                           if_else(padj < 0.05 & log2FoldChange > 0, "Upregulated", "NS")))

DE_mRNAs <- DE_mRNAs %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange < 0, "Downregulated",
                           if_else(padj < 0.05 & log2FoldChange > 0, "Upregulated", "NS")))
# Downregulated miRNAs (should release repression, allowing targets to go UP)
down_miRNAs <- miRNAs %>%
  filter(sign_DE == "Downregulated") %>%
  pull(V1)

# Upregulated miRNAs (should enhance repression, so targets go DOWN)
up_miRNAs <- miRNAs %>%
  filter(sign_DE == "Upregulated") %>%
  pull(V1)


# Get target genes of each miRNA group
down_targets <- filtered_targets %>%
  filter(mature_mirna_id %in% down_miRNAs) %>%
  distinct(target_ensembl)

up_targets <- filtered_targets %>%
  filter(mature_mirna_id %in% up_miRNAs) %>%
  distinct(target_ensembl)

# Genes targeted by downregulated miRNAs
down_miRNA_gene_targets <- DE_mRNAs %>%
  filter(V1 %in% down_targets$target_ensembl)

# Genes targeted by upregulated miRNAs
up_miRNA_gene_targets <- DE_mRNAs %>%
  filter(V1 %in% up_targets$target_ensembl)

colours_1 <- c("Sign_Up" = "#994455", "Sign_Down" = "#0077BB", "NS" = "#BBBBBB") 

opacity <- c("Sign_Up" = 1, "Sign_Down" = 1, "NS" = 1)


p1 <- ggplot(down_miRNA_gene_targets, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = log2FoldChange > 0), alpha = 1) +
  scale_color_manual(values = c("TRUE" = "#994455", "FALSE" = "grey")) +
  labs(title = "Targets of Downregulated mature and hairpin miRNAs (Expected ↑)",
       x = "log2 Fold Change (mRNA)",
       y = "-log10 adjusted p-value") +
  theme_minimal()

p2 <- ggplot(up_miRNA_gene_targets, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = log2FoldChange < 0), alpha = 1) +
  scale_color_manual(values = c("TRUE" = "#0077BB", "FALSE" = "grey")) +
  labs(title = "Targets of Upregulated mature and hairpin miRNAs (Expected ↓)",
       x = "log2 Fold Change (mRNA)",
       y = "-log10 adjusted p-value") +
  theme_minimal()

#ggsave("./3_results/upregulated_miRNA_targets_hairpinsandmature.png",
#       plot = p1,
#       width = 6, height = 5, dpi = 300)
#ggsave("./3_results/downregulated_miRNA_hairpinsandmature.png",
#       plot = p2,
#       width = 6, height = 5, dpi = 300)


###########


genesandTEs <- fread(file = "./3_results/april_DESeq_genesandTEs.csv")


genesandTEs <- genesandTEs %>%
  filter(padj < 0.05)


genesandTEs <- genesandTEs %>%
  filter(grepl("ENSMUS", ensembl.TE))


miRNA_targets_up <- fread(file = "./3_results/multimiR_mature_mirnas_up_validated.csv")
miRNA_targets_down <- fread(file = "./3_results/multimiR_mature_mirnas_down_validated.csv")

miRNA_targets <- rbind(miRNA_targets_up,miRNA_targets_down)

intersect(genesandTEs$hgnc.TE, miRNA_targets$target_symbol)
#[1] "Col1a1"    "Capns1"    "Ssr4"      "Kmt2a"     "Arid1a"    "Piezo1"    "Psmb6"     "Chd3"      "Dusp6"     "Timp3"     "Ccng1"     "Rack1"     "Sec14l1"  
#[14] "Nme2"      "Pfkp"      "Hivep1"    "Tkt"       "Bcl6"      "Adcy6"     "Chd1"      "Svil"      "Ehbp1l1"   "Sorbs1"    "Gaa"       "Itih5"     "Igfbp5"   
#[27] "Pam"       "Atp2b4"    "Col5a1"    "Abca2"     "Fbn1"      "Sord"      "Tpx2"      "Sec62"     "Uox"       "Rps20"     "Coro2a"    "Plin2"     "Rps6"     
#[40] "Hspg2"     "Map7d1"    "Epb41"     "Mmp17"     "Ncor2"     "Polr1d"    "Col1a2"    "Bhlhe40"   "Cand2"     "Pianp"     "Msn"       "Flna"      "Rab11fip1"
#[53] "Col4a1"    "Col4a2"    "Cyb5b"     "Nxpe2"     "Map4"      "Arhgef17"  "Prrg3"     "Jcad"      "Ppp1r18"   "Rps27l"    "Egr1"      "Rps21"     "Prune2"   
#[66] "Srrm2"     "Lrp1"      "Spen"      "Fmod"      "Pgm5"      "Selenok"   "Snai1"     "Rps2"      "Megf8"     "S1pr1"     "Rpl18a"    "Wnk1"      "Camk2n1"  
#[79] "Kmt2d"     "F2r"       "Vamp8"     "Klf13"     "Tmco1"     "Sec61b"    "Smagp"     "Capn12"    "Tead1"     "Tns1"      "Palld"     "Rpl5"      "Rpl18"    
#[92] "Kif13b"    "Tpt1"      "H2-K1"     "Fbln2"     "Slfn9"     "Atxn1l"    "Fryl"      "C4b"       "Ppp1r12b"  "Heg1"      "Plekhn1"   "Zfp970"    "Lrrc32"   
#[105] "Ccl21a"   


common_ids <- intersect(genesandTEs$hgnc.TE, miRNA_targets$target_symbol)

filtered_targets <- miRNA_targets %>%
  filter(target_symbol %in% common_ids)

nrow(filtered_targets)
miRNA_target_numbers <- table(filtered_targets$mature_mirna_id)

#write.csv(miRNA_target_numbers, file = "./3_results/oct_DE_matureandhairpin_miRNA_hits_DE_genes_freq.csv")

miRNA_target_numbers <- as.data.frame(miRNA_target_numbers)
head(miRNA_target_numbers)



gene_miRNA_df <- filtered_targets %>%
  group_by(target_symbol) %>%
  summarise(
    miRNAs = paste(unique(mature_mirna_id), collapse = ", "),
    n_miRNAs = n_distinct(mature_mirna_id)
  ) %>%
  arrange(desc(n_miRNAs))  # optional: sort by number of miRNAs




barplot <- ggplot(gene_miRNA_df, aes(x = n_miRNAs)) +
  geom_bar(fill = "#994455") +
  labs(
    title = "Number of DE miRNAs Targeting Each DE Gene in the mouse SV",
    x = "Number of DE miRNAs targeting the same DE gene",
    y = "Number of genes"
  ) +
  theme_minimal()

plot(barplot)
#ggsave(
#  filename = "./2_figures/miRNA_distribution_plot_DEgenes.png", 
#  plot     = barplot,
#  width    = 8,    # in inches
#  height   = 6,    # in inches
#  dpi      = 300   # resolution
#)


