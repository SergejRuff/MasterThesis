rm(list = ls())

library(Virusparies)
library(tidyverse)
library(gt)

path <- "data/RNAvirus_Mammals_newJan2023/mammals/orthomyxo_20july_RNAvirus_nofil_1/virusgatherer-cap3.tsv"
out <- "output/mammals/Orthomyxo/plots/rnavirus/Gatherer/"

file <- ImportVirusTable(path)

file <- VhgPreprocessTaxa(file,taxa_rank = "Family")

file_filtered <- VhgSubsetHittable(file,ViralRefSeq_E_criteria = 1e-5)

result <- file_filtered %>%
  # Extract the part of the string after the first '|'
  mutate(ViralRefSeq_subject = str_extract(ViralRefSeq_subject, "(?<=\\|).*")) %>%
  # Prefix each entry with `//` and replace pipe delimiter with newline characters
  mutate(ViralRefSeq_subject = paste0("//", str_replace_all(ViralRefSeq_subject, "\\|", "\n// "))) %>%
  # Group by taxonomy and concatenate subjects with newlines
  group_by(ViralRefSeq_taxonomy) %>%
  summarise(ViralRefSeq_subjects = paste(unique(ViralRefSeq_subject), collapse = "\n"))%>%
  separate_rows(ViralRefSeq_subjects, sep = "//")%>%
  ungroup() %>% 
  filter(ViralRefSeq_subjects != "")


gt_table <- VhgTabularRasa(result,title = "Viral Subject",names_ = c("Reference taxonomy","Viral subject"))

gt_table

ExportVirusGt(gt_table,filename = "subject.png",path = out,export_gt_obj = TRUE,create.dir = TRUE)
ExportVirusGt(gt_table,filename = "subject.docx",path = out,export_gt_obj = TRUE,create.dir = TRUE)