rm(list = ls())

library(Virusparies)
library(tidyverse)

vg <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/combined_ga_with_source.tsv")

vg <- VhgAddPhylum(vg,"ViralRefSeq_taxonomy")

valid_phyla_rna <- c("Ambiviricota", "Duplornaviricota", "Kitrinoviricota",
                     "Lenarviricota", "Negarnaviricota", "Pisuviricota")

vg_top10_contiglen <- vg %>%
  filter(Phylum %in% valid_phyla_rna) %>%
  filter(ViralRefSeq_E < 1e-5) %>%
  arrange(Phylum, desc(contig_len)) %>%
  group_by(Phylum) %>%
  distinct(ViralRefSeq_taxonomy, .keep_all = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()
  

VhgBoxplot(vg_top10_contiglen,x_column = ViralRefSeq_taxonomy,y_column = contig_len,group_unwanted_phyla = "rna",reorder_criteria = "phylum_median")


modified_contig_ids <- sub("_cap3_Contig-[0-9]+$", "", vg_top10_contiglen$contig_id)

# View the modified contig IDs
print(modified_contig_ids)

writeLines(modified_contig_ids, "contig_ids_top10.txt")

ExportVirusDataFrame(vg_top10_contiglen,"longest_con_gvtable.tsv",dir_path = "output/mammals/")
