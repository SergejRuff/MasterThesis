rm(list = ls())

library(Virusparies)
library(tidyverse)

file_paths <- list.files(path = "/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/data/RNAvirus_Mammals_newJan2023/mammals/",
                         pattern = "virushunter.tsv", 
                         full.names = TRUE, 
                         recursive = TRUE)

# Use the existing ImportVirusTable function to read each file and store them in a list
virus_tables <- lapply(file_paths, ImportVirusTable)

# Combine the individual data frames into one
combined_virus_table <- do.call(rbind, virus_tables)

# Check the first few rows (optional)
head(combined_virus_table)

combined_virus_table <- VhgSubsetHittable(combined_virus_table,ViralRefSeq_E_criteria = 1e-5)