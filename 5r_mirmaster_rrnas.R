
# Purpose: To construct a reusable pipeline to convert raw miRMaster output 
# (excel spreadsheet, xlsx file) into a format suitable for DE analysis using 
# the package DESeq2. Data features 8 biological replicates labelled as group 
# A (test) or P (control).

# Example constructed using rRNA data. Edit as needed. 

# Author: Chishan Burch
# Contact: chishanburch@gmail.com
# Date: 13/01/2025


##################################

# Install these if you don't have them already
#install.packages("purrr")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("openxlsx")
#install.packages("tibble")

# Load packages
library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(openxlsx)
library(purrr)

# Take the miRMaster output from https://drive.google.com/drive/folders/1b7CZSGZ4DOdIXZTOaLohfrFm2gSaNs0L
# In this example we have just downloaded raw rRNA data from miRMAaster (rRNA_quantification_raw.xlsx)
# and stored it in C:/Users/jabry/OneDrive/UniRstudio/thesis_work/1_data

rRNA_data <- read.xlsx("./1_data/rRNA_quantification_raw.xlsx")

# We will create individual objects for each replicate, and name them A2, A3 etc. 

A2 <- rRNA_data %>%
  select(Reference, A2)%>%
  rename(Count = A2)
A3 <- rRNA_data %>%
  select(Reference, A3)%>%
  rename(Count = A3)
A4 <- rRNA_data %>%
  select(Reference, A4)%>%
  rename(Count = A4)
A5 <- rRNA_data %>%
  select(Reference, A5)%>%
  rename(Count = A5)

P1 <- rRNA_data %>%
  select(Reference, P1)%>%
  rename(Count = P1)
P3 <- rRNA_data %>%
  select(Reference, P3)%>% 
  rename(Count = P3)
P4 <- rRNA_data %>%
  select(Reference, P4)%>%
  rename(Count = P4)
P5 <- rRNA_data %>% 
  select(Reference, P5)%>%
  rename(Count = P5)

# Start by removing the unneeded columns and aggregating data. We only want to
# keep the id and count columns because that's what DESeq2 needs.


A2_agg <- A2 %>%
  select(Reference, Count) %>% 
  group_by(Reference) %>%
  summarize(Count = sum(Count))%>%
  rename(A2 = Count)

A3_agg <- A3 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count)) %>%
  rename(A3 = Count)

A4_agg <- A4 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count)) %>%
  rename(A4 = Count)

A5_agg <- A5 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count)) %>%
  rename(A5 = Count)

P1_agg <- P1 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count)) %>%
  rename(P1 = Count)

P3_agg <- P3 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count))  %>%
  rename(P3 = Count)

P4_agg <- P4 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count))  %>%
  rename(P4 = Count)

P5_agg <- P5 %>%
  select(Reference, Count) %>%
  group_by(Reference) %>%
  summarize(Count = sum(Count))  %>%
  rename(P5 = Count)

# List of all aggregated data frames
data_frames <- list(A2_agg, A3_agg, A4_agg, A5_agg, P1_agg, P3_agg, P4_agg, P5_agg)

# Full join all data frames
df_combined <- reduce(data_frames, full_join, by = "Reference")

# Print the result
print(df_combined)

# save that
#write.csv(df_combined, file = "./1_data/miRMaster_raw_formatted/miRmaster_rRNAs_raw.csv", row.names = FALSE)

# Reuse and repeat this for other sncRNA types
# Get rid of NAs in next step when we combine them all and do DESeq2.
