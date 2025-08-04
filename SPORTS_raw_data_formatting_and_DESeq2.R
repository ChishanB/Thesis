library(data.table)
library(dplyr)
library(tibble)
library(DESeq2)
#library(ggplot2)
library(openxlsx)
#library(pheatmap)
#library(ggrepel)
library(edgeR)

library(data.table)



#data <- fread("./1_data/SPORTS_data/Schjenken_A2_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_A3_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_A4_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_A5_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_P1_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_P3_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_P4_summary.txt", header = FALSE)
#data <- fread("./1_data/SPORTS_data/Schjenken_P5_summary.txt", header = FALSE)

# Remove rows containing 'unmatch' in a specific column
data <- data[!grepl("Unmatch", data$V1), ]

data <- data[-(1:4), -1]

data <- data %>%
  filter(V2 != "" & !is.na(V2))

data <- data %>%
  filter(V2 != "-")

data <- as.data.frame(data)

colnames(data) <- NULL

colnames(data) <- c("sncRNA", "Count")

#write.csv(data, file = "./1_data/SPORTS_data/P5_raw.csv", row.names = FALSE)


# Now we need to put all of the sample info together potentially with a join
A2_data <- read.csv("./1_data/SPORTS_data/A2_raw.csv")
A3_data <- read.csv("./1_data/SPORTS_data/A3_raw.csv")
A4_data <- read.csv("./1_data/SPORTS_data/A4_raw.csv")
A5_data <- read.csv("./1_data/SPORTS_data/A5_raw.csv")
P1_data <- read.csv("./1_data/SPORTS_data/P1_raw.csv")
P3_data <- read.csv("./1_data/SPORTS_data/P3_raw.csv")
P4_data <- read.csv("./1_data/SPORTS_data/P4_raw.csv")
P5_data <- read.csv("./1_data/SPORTS_data/P5_raw.csv")

# Rename "Count" to the sample name
colnames(A2_data)[colnames(A2_data) == "Count"] <- "A2"
colnames(A3_data)[colnames(A3_data) == "Count"] <- "A3"
colnames(A4_data)[colnames(A4_data) == "Count"] <- "A4"
colnames(A5_data)[colnames(A5_data) == "Count"] <- "A5"

colnames(P1_data)[colnames(P1_data) == "Count"] <- "P1"
colnames(P3_data)[colnames(P3_data) == "Count"] <- "P3"
colnames(P4_data)[colnames(P4_data) == "Count"] <- "P4"
colnames(P5_data)[colnames(P5_data) == "Count"] <- "P5"


# Repeat the process for other samples (data_a3, data_p1, etc.) and rename the 'Count' column accordingly

# Merge all sample data frames by sncRNA
combined_data <- A2_data %>%
  full_join(A3_data, by = "sncRNA") %>%
  full_join(A4_data, by = "sncRNA") %>%
  full_join(A5_data, by = "sncRNA") %>%
  full_join(P1_data, by = "sncRNA") %>%
  full_join(P3_data, by = "sncRNA") %>%
  full_join(P4_data, by = "sncRNA") %>%
  full_join(P5_data, by = "sncRNA")

# View combined data
head(combined_data)

#write.csv(combined_data, file = "./1_data/SPORTS_data/raw_sports_counts.csv", row.names = FALSE)

combined_data <- column_to_rownames(combined_data, var = "sncRNA")

sample_metadata <- fread("./1_data/sncRNA_sample_metadata.csv") %>% 
  as.data.frame() %>%
  mutate(Treatment_group = factor(Treatment_group, levels = c("Control", "Acr")))


combined_data <- na.omit(combined_data)

min.counts <- 5
min.samples <- 3

# Initial filter
combined_data <- combined_data %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# DESeq2 needs integers
combined_data <- round(combined_data)

# After formatting, do the DESEq2 and save that as a spreadsheet. Send this to Sean with the miRMaster one

dds <- DESeqDataSetFromMatrix(countData = combined_data, # Needs to be raw counts with rownames as RNAs
                              colData = sample_metadata, #data frame of meta data (rownames are samples)
                              desig = ~ Treatment_group) # Column in metadata for treatment


dds <- DESeq(dds)

res <- results(dds) 

res_df <- as.data.frame(res) # Coerce into a more easily accessible format.

#write.csv(res_df, file = "./1_data/SPORTS_data/SPORTS_DESeq2_output.csv", row.names = TRUE)

# Filter dds results based on a significance threshold of p-adjusted < 0.05.
res_df_significant <- res_df %>%
  filter(padj <= 0.05)


res_df_significant2 <- res_df %>%
  filter(pvalue <= 0.05)

# Save dds results meeting significance threshold. 
 #write.csv(res_df_significant, file = "./1_data/SPORTS_data/SPORTS_DESeq2_output_under0.05.csv", row.names = TRUE)









