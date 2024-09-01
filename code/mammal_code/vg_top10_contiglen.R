rm(list = ls())

library(Virusparies)
library(tidyverse)

vg <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/combined_ga.tsv")

vg <- VhgAddPhylum(vg,"ViralRefSeq_taxonomy")

valid_phyla_rna <- c("Ambiviricota", "Duplornaviricota", "Kitrinoviricota",
                     "Lenarviricota", "Negarnaviricota", "Pisuviricota")

vg_top10_contiglen <-vg %>%
  filter(Phylum %in% valid_phyla_rna) %>%          
  group_by(Phylum) %>%
  filter(ViralRefSeq_E < 1e-5) %>%
  arrange(Phylum, desc(contig_len)) %>%
  slice_head(n = 10) %>%
  ungroup()
  

VhgBoxplot(vg,x_column = ViralRefSeq_taxonomy,y_column = contig_len,group_unwanted_phyla = "rna",reorder_criteria = "phylum_median")

