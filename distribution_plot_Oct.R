# Date: 23/10/25
# Author: Chishan Burch
# Purpose: 


# Restarting cause i think I did it all wrong.

# I realised my distribution plots were wrong while I was writing
# my discussion and counting up how many targets a particular miRNA had, 
# and saw when comparing the corresponding multimir results to my plots 
#that the numbers were inconsistent. It turns out I was counting in gneral 
#how many unique genes there were and not necessarily unique genes per miRNA. 
#So I regenerated them and changed the colours to be different to the ones 
#I used to indicate up and downregulation. 


# (Instructions from John): Single omics analysis – link mature miRNA to mRNA targets 


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

#miRNAs <- rbind(mature,hairpins)

# miRNAs contains all differentially expressed miRNAs. 



# We want to know how many target the same genes.

mature_up <- fread(file = "./3_results/multimiR_mature_mirnas_up_validated.csv")
mature_down <- fread(file = "./3_results/multimiR_mature_mirnas_down_validated.csv")
#hairpins_up <- fread(file = "./3_results/multimiR_hairpin_mirnas_up_validated.csv")
#hairpins_down <- fread(file = "./3_results/multimiR_hairpin_mirnas_down_validated.csv")

# Decide whether u want hairpins included
#miRNA_targets <- rbind(mature_up,mature_down,hairpins_up,hairpins_down)
miRNA_targets <- rbind(mature_up,mature_down)
miRNA_targets <- miRNA_targets[ ,!"V1"]

DE_mature_miRNA_multimiR_results <- fread(file = "./3_results/oct_multimiR_DE_targets_of_DE_mature_miRNA.csv")
DE_mature_miRNA_multimiR_results <- DE_mature_miRNA_multimiR_results[ ,!"V1"]

print(miRNA_targets) # for general distribution plot
print(DE_mature_miRNA_multimiR_results) # for DE miRNA x DE mRNA distribution plot

############## General distribution plot

# Remove duplicates
miRNA_targets_unique <- miRNA_targets %>%
  distinct(mature_mirna_id, target_symbol)

# Count unique targets per miRNA
miRNA_target_counts <- miRNA_targets_unique %>%
  group_by(mature_mirna_id) %>%
  summarise(n_unique_targets = n()) %>%
  arrange(n_unique_targets)  # optional: sort by target count

# Horizontal bar plot
p1 <- ggplot(miRNA_target_counts, aes(x = reorder(mature_mirna_id, n_unique_targets), y = n_unique_targets)) +
  geom_bar(stat = "identity", fill = "#cc5500") +
  coord_flip() +  # flips axes
  labs(
    title = "Number of Unique Targets per miRNA",
    x = "miRNA",
    y = "Number of unique target genes"
  ) +
  theme_minimal()

#ggsave(
#  filename = "./2_figures/miRNA_distribution_plot_matureDEmiRNAs_x_allmRNAs.png", 
#plot     = p1,
#  width    = 8,    # in inches
#height   = 6,    # in inches
# dpi      = 300   # resolution
#)


###################### DE miRNAs linking to DE mRNAs

# Remove duplicates
DE_mature_miRNA_multimiR_results_unique <- DE_mature_miRNA_multimiR_results %>%
  distinct(mature_mirna_id, target_symbol)

# Count unique targets per miRNA
DE_mature_miRNA_multimiR_results_unique_counts <- DE_mature_miRNA_multimiR_results_unique %>%
  group_by(mature_mirna_id) %>%
  summarise(n_unique_targets = n()) %>%
  arrange(n_unique_targets)  # optional: sort by target count

# Horizontal bar plot
p2 <- ggplot(DE_mature_miRNA_multimiR_results_unique_counts, aes(x = reorder(mature_mirna_id, n_unique_targets), y = n_unique_targets)) +
  geom_bar(stat = "identity", fill = "#4B0082") +
  coord_flip() +  # flips axes
  labs(
    title = "Number of Unique Targets per miRNA",
    x = "miRNA",
    y = "Number of unique target genes"
  ) +
  theme_minimal()

#ggsave(
#  filename = "./2_figures/miRNA_distribution_plot_matureDEmiRNAs_x_DEmRNAs.png", 
#  width    = 8,    # in inches
#height   = 6,    # in inches
 #dpi      = 300   # resolution
#)
