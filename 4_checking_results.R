# Purpose: To search the results generated in 3_Tetranscripts_plots.R for DE genes 
# encoding PIWI proteins, and generate a table showing whether they were DE as 
# expected

#############################################

# Load packages
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

# Read in the DESeq2 output (filtered for 5 across 3, includes TE and non-TE genes)

data <- read.csv("./3_results/clean_joined_data.csv")

ago_piwi_etc <- data %>%
  filter(grepl("piwi|ago|hsp90|miwi|Fkbp6|Zfp809|KRAB", mouse.TE, ignore.case = TRUE))

# Remove columns that arent needed
ago_piwi_etc <- ago_piwi_etc[ ,-c(1, 2, 4, 5, 6, 9)]
print(ago_piwi_etc)


#plugging in names of known repressors
ERV_repressors <- data %>%
  filter(grepl("Zfp809|Zfp819|Zfp932|Zfp708|Gm15446|Zfp42|YY1|Oct4|Myc|Zfp281|Rap1|Tin2", mouse.TE, ignore.case = TRUE))

# Remove columns that arent needed
ERV_repressors <- ERV_repressors[ ,-c(1, 2, 4, 5, 6, 9)]
print(ERV_repressors)

ERV_repressors_filtered <- ERV_repressors %>% filter(padj < 0.05)
print(ERV_repressors_filtered)
