# Purpose: Assign subtype labels to deseq2 data using grep() so I can make a pie chart in excel
# 


library(dplyr)
library(openxlsx)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)


deseq_result <- fread(file = "./3_results/july_DESeq_sncRNAs.csv")
deseq_result <- arrange(deseq_result, padj)
head(deseq_result)


x <- deseq_result %>%
  mutate(Annotation = word(V1, -1)) %>%  
  mutate(
    Subtype = case_when(
      grepl("hairpin", Annotation) ~ "hairpin",
      grepl("mature", Annotation)  ~ "mature",
      grepl("piRNA", Annotation)   ~ "piRNA",
      grepl("snoRNA", Annotation)  ~ "snoRNA",
      grepl("snRNA", Annotation)   ~ "snRNA",
      grepl("rRNA", Annotation)    ~ "rRNA",
      grepl("Rfam other", Annotation) ~ "Rfam other than sncRNA",
      TRUE ~ Annotation  
    )
  )

write.xlsx(x, file = "./3_results/july_DESeq_sncRNAs_labelled.xlsx")

