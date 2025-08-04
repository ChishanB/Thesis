# Date: 13.1.25
# Author: Chishan
# Contact: chishanburch@gmail.com
# Purpose: Construct a volcano plot from miRMaster data using our DE pipeline 

# Data sourced from:
#  https://drive.google.com/drive/folders/1e1354MEyKLuyauwAsU5VEiHUWBqFU2T-
# (google drive) / Data - Chishan Burch / sncRNA / sncRNA_sample_metadata.csv
# (google drive) / Data - Chishan Burch / sncRNA / miRMaster_run_6.8.24


###############################################################################
library(data.table)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(openxlsx)

miRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_miRNAs_raw.csv")
piRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_piRNAs_raw.csv")
rRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_rRNAs_raw.csv")
snoRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_snoRNAs_raw.csv")
tRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_tRNAs_raw.csv")
snRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_snRNAs_raw.csv")
miscRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_miscRNAs_raw.csv")
scaRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_scaRNAs_raw.csv")
circRNA <- read.csv("./1_data/miRMaster_raw_formatted/miRmaster_circRNAs_raw.csv")

# Stack those data frames
combined_data <- rbind(miRNA, piRNA, rRNA, snoRNA, tRNA, snRNA, miscRNA, scaRNA, circRNA)
##might need to add another column tagging RNA type if I ever want to colour scatter plot based on rna subtype

# Lets filter for rnas with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

# Filtering step

filtered_data <- combined_data %>%
  filter(rowSums(combined_data >= min.counts) >= min.samples)


#when combined, decide whether its more appropriate to change NAs to 0 or remove all rows containing 0
#----- rework above

miRNA <- fread("./1_data/miRmaster_miRNAs_raw.csv")
sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

# Set the columns as the row names with tibble
miRNA <- 
  miRNA %>%
  data.frame() %>%
  column_to_rownames(var = "miRNA")

# Define thresholds 
min.counts <- 5
min.samples <- 3

# Initial filter
miRNA_filtered <- miRNA %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Remove rows containing NA values
if (any(is.na(miRNA))) {
  miRNA <- miRNA[complete.cases(miRNA), ]
}

# DE analysis
dds <- DESeqDataSetFromMatrix(countData = miRNA,
                              colData = sample_metadata,
                              desig = ~ Treatment_group)

dds <- DESeq(dds)


res <- results(dds, contrast = c("Treatment_group","Acr","Control")) 

levels(sample_metadata$Treatment_group)

res_df <- as.data.frame(res) # Coerce into a more easily accessible format.

# Remove rows containing NA values
if (any(is.na(res_df))) {
  res_df <- res_df[complete.cases(res_df), ]
}

#VST matrix
vsd <- varianceStabilizingTransformation(dds, blind = T) 
vst_matrix <- assay(vsd)

#
# Add a column labelling whether miRNAs are upregulated or downregulated
res_df <- res_df %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange < -0.585, "Sign_Down",
                           if_else(padj < 0.05 & log2FoldChange > 0.585, "Sign_Up", "NS"))) 


# Create a -log10 column for y-axis of volcano plots
res_df <- res_df %>% mutate(log_padj = -log10(padj))



# Declare colours
colours_1 <- c("Sign_Up" = "#AA3377", "Sign_Down" = "#228833", "NS" = "#BBBBBB") 
opacity <- c("Sign_Up" = 1, "Sign_Down" = 1, "NS" = 0.5)

# Declare padj threshold
p_adj_threshold <- 0.05
options(ggrepel.max.overlaps = Inf)

# Plot
plot <- res_df %>%
  filter(is.na(sign_DE) != TRUE) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sign_DE, alpha = sign_DE, fill = sign_DE)) +
  geom_point(size = 5, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + theme_bw() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 

