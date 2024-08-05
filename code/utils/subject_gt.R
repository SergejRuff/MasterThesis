rm(list = ls())

library(Virusparies)
library(tidyverse)
library(gt)

path <- "data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv"
out <- "output/TaubertDatacombined/plots/smalldna/Gatherer/"

file <- ImportVirusTable(path)

file <- VhgPreprocessTaxa(file,taxa_rank = "Family")

file_filtered <- VhgSubsetHittable(file,ViralRefSeq_E_criteria = 1e-5)

# Process and count matching subjects
result <- file_filtered %>%
  mutate(ViralRefSeq_subject = str_extract(ViralRefSeq_subject, "(?<=\\|).*")) %>%
  mutate(ViralRefSeq_subject = str_replace_all(ViralRefSeq_subject, "\\|", "\n")) %>%
  group_by(ViralRefSeq_taxonomy, ViralRefSeq_subject) %>%
  summarise(subject_count = n(), .groups = 'drop')


gt_table <- VhgTabularRasa(result,title = "Viral Subject",names_ = c("Reference taxonomy","Viral subject","counts"))

gt_table

ExportVirusGt(gt_table,filename = "subject.png",path = out,export_gt_obj = TRUE,create.dir = TRUE)
ExportVirusGt(gt_table,filename = "subject.docx",path = out,export_gt_obj = TRUE,create.dir = TRUE)



