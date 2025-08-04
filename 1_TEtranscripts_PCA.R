####################################################################################

#date: 25/6/24
#purpose: To import TEtranscripts raw output, construct PCA plots,
#email: chishanburch@gmail.com

#packages########

# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

# Read the data
a <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")
b <- read.delim("./1_data/150bp/raw_counts/TE_Class_only_counts.txt")
#c <- read.delim("./1_data/150bp/raw_counts/TE_Family_Class_only_counts.txt")
#d <- read.delim("./1_data/150bp/raw_counts/TE_only_and_all_counts.txt")
sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

# Rename 'sample.IDs' to 'SampleID'
sample_metadata <- sample_metadata %>% 
  rename(SampleID = sample.IDs)


# Check column names and structure of sample_metadata
print(names(sample_metadata))
print(head(sample_metadata))

# Assuming the column is named correctly, assign group labels
# Change 'SampleID' to the correct column name if it differs

if ("SampleID" %in% colnames(sample_metadata)) {
  sample_metadata$Treatment_group <- ifelse(grepl('p|P', sample_metadata$SampleID), 'Control', 'Acr')
} else {
  stop("Column 'SampleID' not found in sample_metadata.")
}

# Set colors for the groups
level_colors <- c("Acr" = "#CC79A7", "Control" = "#009E73")
level_colors_darker <- c("Acr" = "#663C53", "Control" = "#004f39")

# Preprocess data
b <- b %>% select(-1)

# Perform PCA
test_PCA <- prcomp(t(b), center = FALSE)

# Create PCA data frame
pca_data <- data.frame(test_PCA$x)
pca_data$SampleID <- rownames(pca_data)

# Ensure SampleID is a common column for merging
if ("SampleID" %in% colnames(pca_data) & "SampleID" %in% colnames(sample_metadata)) {
  pca_data <- merge(pca_data, sample_metadata, by = "SampleID")
} else {
  stop("Column 'SampleID' not found in one of the data frames.")
}

####
test_PCA <- prcomp(pca_data[, c("PC1", "PC2")], center = TRUE, scale. = TRUE)

# Extract percentage of variance explained
variance_explained <- 100 * test_PCA$sdev^2 / sum(test_PCA$sdev^2)
####

# Create PCA plot using ggplot2
test_PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment_group, label = SampleID)) +
  geom_point(size = 3) +
  geom_label_repel(size = 3, box.padding = 0.25, label.padding = 0.5) +
  scale_color_manual(values = level_colors_darker) +
  labs(color = "Treatment group") +  labs(color = "Treatment group", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"), y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  guides(color = guide_legend(title = "Treatment group")) +
  theme_minimal() + stat_ellipse(aes(group = Treatment_group), alpha = 0.2, linetype = 1) 
 
# Print the plot
print(test_PCA_plot) #SOMETIMES SAMPLES A2 DISAPPEARS WHEN I RUN THIS AND I DONT KNOW WHY

ggsave(plot = test_PCA_plot, filename = "./2_figures/TEtranscripts_PCA.png", width = 8, height = 6, dpi = 300)          


#############################
#b doesnt have the names of the TEs, but a does. Can I repeat this process with a?
# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

# Rename 'sample.IDs' to 'SampleID'
sample_metadata <- sample_metadata %>% 
  rename(SampleID = sample.IDs)


# Check column names and structure of sample_metadata
print(names(sample_metadata))
print(head(sample_metadata))

# Assuming the column is named correctly, assign group labels
# Change 'SampleID' to the correct column name if it differs

if ("SampleID" %in% colnames(sample_metadata)) {
  sample_metadata$Treatment_group <- ifelse(grepl('p|P', sample_metadata$SampleID), 'Control', 'Acr')
} else {
  stop("Column 'SampleID' not found in sample_metadata.")
}


# Read the data
a <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")

# Remove the columns from a which are not counts or name info
a <- a %>% select(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)


# Lets filter for TEs with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

# Set the columns as the row names
a <- 
  a %>%
  data.frame() %>%
  column_to_rownames(var = "ensembl.TE")

# Filtering step
a_filtered <- a %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Start plotting
# Perform PCA
test_PCA <- prcomp(t(a_filtered), center = FALSE)

# Create PCA data frame
pca_data <- data.frame(test_PCA$x)
pca_data$SampleID <- rownames(pca_data)

# Ensure SampleID is a common column for merging
if ("SampleID" %in% colnames(pca_data) & "SampleID" %in% colnames(sample_metadata)) {
  pca_data <- merge(pca_data, sample_metadata, by = "SampleID")
} else {
  stop("Column 'SampleID' not found in one of the data frames.")
}

# Set colors for the groups
level_colors <- c("Acr" = "#CC79A7", "Control" = "#009E73")
level_colors_darker <- c("Acr" = "#663C53", "Control" = "#004f39")

# Plot 
a_filtered <- prcomp(pca_data[, c("PC1", "PC2")], center = TRUE, scale. = TRUE)

# Extract percentage of variance explained
variance_explained <- 100 * a_filtered$sdev^2 / sum(a_filtered$sdev^2)


# Create PCA plot using ggplot2
TE_PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment_group, label = SampleID)) +
  geom_point(size = 3) +
  geom_label_repel(size = 3, box.padding = 0.25, label.padding = 0.5) +
  scale_color_manual(values = level_colors_darker) +
  labs(color = "Treatment group") +  labs(color = "Treatment group", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"), y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  guides(color = guide_legend(title = "Treatment group")) +
  theme_minimal() + stat_ellipse(aes(group = Treatment_group), alpha = 0.2, linetype = 1) 

# Print the plot
print(TE_PCA_plot) 

#ggsave(plot = TE_PCA_plot, filename = "./2_figures/TE_PCA_plot.png", width = 8, height = 6, dpi = 300)          

####TEs alone#######################################################################################
# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

# Rename 'sample.IDs' to 'SampleID'
sample_metadata <- sample_metadata %>% 
  rename(SampleID = sample.IDs)


# Check column names and structure of sample_metadata
print(names(sample_metadata))
print(head(sample_metadata))

# Assuming the column is named correctly, assign group labels
# Change 'SampleID' to the correct column name if it differs

if ("SampleID" %in% colnames(sample_metadata)) {
  sample_metadata$Treatment_group <- ifelse(grepl('p|P', sample_metadata$SampleID), 'Control', 'Acr')
} else {
  stop("Column 'SampleID' not found in sample_metadata.")
}


### can i just look at TEs alone???
e <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")

TE_subset <- subset(e, grepl('TE', Type, ignore.case = TRUE))

# Remove the columns from a which are not counts or name info
#TE_subset <- TE_subset %>% select(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)

# Remove the columns from a which are not counts or name info
TE_subset <- TE_subset %>% select(-1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)


# Lets filter for TEs with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

rownames(TE_subset) <- NULL

# Set the columns as the row names
TE_subset <- 
 TE_subset %>%
  data.frame() %>%
 column_to_rownames(var = "Name")

#rownames(TE_subset) <- TE_subset$Names

# Filtering step
TE_subset_filtered <- TE_subset %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Remove rows containing NA values
if (any(is.na(TE_subset_filtered))) {
  TE_subset_filtered <- TE_subset_filtered[complete.cases(TE_subset_filtered), ]
}

# Start plotting
# Perform PCA
TE_PCA <- prcomp(t(TE_subset_filtered), center = FALSE)

# Create PCA data frame
pca_data <- data.frame(TE_PCA$x)
pca_data$SampleID <- rownames(pca_data)

# Ensure SampleID is a common column for merging
if ("SampleID" %in% colnames(pca_data) & "SampleID" %in% colnames(sample_metadata)) {
  pca_data <- merge(pca_data, sample_metadata, by = "SampleID")
} else {
  stop("Column 'SampleID' not found in one of the data frames.")
}

# Set colors for the groups
level_colors <- c("Acr" = "#CC79A7", "Control" = "#009E73")
level_colors_darker <- c("Acr" = "#663C53", "Control" = "#004f39")

# Plot 
TE_subset_filtered <- prcomp(pca_data[, c("PC1", "PC2")], center = TRUE, scale. = TRUE)

# Extract percentage of variance explained
variance_explained <- 100 * TE_subset_filtered$sdev^2 / sum(TE_subset_filtered$sdev^2)


# Create PCA plot using ggplot2
TE_PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment_group, label = SampleID)) +
  geom_point(size = 3) +
  geom_label_repel(size = 3, box.padding = 0.25, label.padding = 0.5) +
  scale_color_manual(values = level_colors_darker) +
  labs(color = "Treatment group") +  labs(color = "Treatment group", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"), y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  guides(color = guide_legend(title = "Treatment group")) +
  theme_minimal() + stat_ellipse(aes(group = Treatment_group), alpha = 0.2, linetype = 1) 

# Print the plot
print(TE_PCA_plot) 

#ggsave(plot = TE_PCA_plot, filename = "./2_figures/TE_PCA_plot.png", width = 8, height = 6, dpi = 300)          

####genes alone#######################################################################################
# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

# Rename 'sample.IDs' to 'SampleID'
sample_metadata <- sample_metadata %>% 
  rename(SampleID = sample.IDs)


# Check column names and structure of sample_metadata
print(names(sample_metadata))
print(head(sample_metadata))

# Assuming the column is named correctly, assign group labels
# Change 'SampleID' to the correct column name if it differs

if ("SampleID" %in% colnames(sample_metadata)) {
  sample_metadata$Treatment_group <- ifelse(grepl('p|P', sample_metadata$SampleID), 'Control', 'Acr')
} else {
  stop("Column 'SampleID' not found in sample_metadata.")
}


### can i just look at non-TE genes alone???
f <- read.delim("./1_data/150bp/raw_counts/TE_and_gene_counts_Updated_and_Raw.txt")

gene_subset <- subset(f, grepl('gene', Type, ignore.case = TRUE))

# Remove the columns from a which are not counts or name info
#TE_subset <- TE_subset %>% select(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)

# Remove the columns from a which are not counts or name info
gene_subset <- gene_subset %>% select(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)


# Lets filter for TEs with at least 5 counts across 3 samples.
# Set parameters
min.counts <- 5
min.samples <- 3

#rownames(gene_subset) <- NULL

# Set the columns as the row names
gene_subset <- 
  gene_subset %>%
  data.frame() %>%
  column_to_rownames(var = "ensembl.TE")

#rownames(TE_subset) <- TE_subset$Names

# Filtering step
gene_subset_filtered <- gene_subset %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Remove rows containing NA values
if (any(is.na(gene_subset_filtered))) {
  gene_subset_filtered <- gene_subset_filtered[complete.cases(gene_subset_filtered), ]
}

# Start plotting
# Perform PCA
gene_PCA <- prcomp(t(gene_subset_filtered), center = FALSE)

# Create PCA data frame
pca_data <- data.frame(gene_PCA$x)
pca_data$SampleID <- rownames(pca_data)

# Ensure SampleID is a common column for merging
if ("SampleID" %in% colnames(pca_data) & "SampleID" %in% colnames(sample_metadata)) {
  pca_data <- merge(pca_data, sample_metadata, by = "SampleID")
} else {
  stop("Column 'SampleID' not found in one of the data frames.")
}

# Set colors for the groups
level_colors <- c("Acr" = "#CC79A7", "Control" = "#009E73")
level_colors_darker <- c("Acr" = "#663C53", "Control" = "#004f39")

# Plot 
gene_subset_filtered <- prcomp(pca_data[, c("PC1", "PC2")], center = TRUE, scale. = TRUE)

# Extract percentage of variance explained
variance_explained <- 100 * gene_subset_filtered$sdev^2 / sum(gene_subset_filtered$sdev^2)


# Create PCA plot using ggplot2
gene_PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment_group, label = SampleID)) +
  geom_point(size = 3) +
  geom_label_repel(size = 3, box.padding = 0.25, label.padding = 0.5) +
  scale_color_manual(values = level_colors_darker) +
  labs(color = "Treatment group") +  labs(color = "Treatment group", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"), y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  guides(color = guide_legend(title = "Treatment group")) +
  theme_minimal() + stat_ellipse(aes(group = Treatment_group), alpha = 0.2, linetype = 1) 

# Print the plot
print(gene_PCA_plot) 

#ggsave(plot = gene_PCA_plot, filename = "./2_figures/nonTE_gene_PCA_plot.png", width = 8, height = 6, dpi = 300)          

