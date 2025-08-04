
# Purpose: To construct a reusable pipeline to convert raw miRMaster output 
# (excel spreadsheet, xlsx file) into a format suitable for DE analysis using 
# the package DESeq2. Designed to work with 8 samples. 

# Example constructed using piRNA data. Edit as needed. 

# Author: Chishan Burch
# Contact: chishanburch@gmail.com


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

# First, go into each excel file and remove any merged rows. 
# Read in miRMaster outputs for miRNA (xlsx files)
A2 <- read.xlsx("./1_data/miRMaster_data/A2/A2_pirnas.xlsx")
A3 <- read.xlsx("./1_data/miRMaster_data/A3/A3_pirnas.xlsx")
A4 <- read.xlsx("./1_data/miRMaster_data/A4/A4_pirnas.xlsx")
A5 <- read.xlsx("./1_data/miRMaster_data/A5/A5_pirnas.xlsx")
P1 <- read.xlsx("./1_data/miRMaster_data/P1/P1_pirnas.xlsx")
P3 <- read.xlsx("./1_data/miRMaster_data/P3/P3_pirnas.xlsx")
P4 <- read.xlsx("./1_data/miRMaster_data/P4/P4_pirnas.xlsx")
P5 <- read.xlsx("./1_data/miRMaster_data/P5/P5_pirnas.xlsx")


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

# Full join all data frames on 'miRNA'
df_combined <- reduce(data_frames, full_join, by = "Reference")

# Print the result
print(df_combined)

# save that
#write.csv(df_combined, file = "./1_data/miRMaster_raw_formatted/miRmaster_piRNAs_raw.csv", row.names = FALSE)

# Reuse and repeat this for other sncRNA types
