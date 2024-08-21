rm(list = ls())

library(Virusparies)
library(gt)
# set sys env for brave executable for gtsave 
# Sys.setenv(
#   CHROMOTE_CHROME = "C:/Program Files/BraveSoftware/Brave-Browser/Application/brave.exe"
# )

rna <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virusgatherer-cap3.tsv")
smalldna <- ImportVirusTable("data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv")
largedna <- ImportVirusTable("data/hittables_taubert/largedna/combined_virusgatherer-cap3.tsv")

rna <- VhgSubsetHittable(rna,ViralRefSeq_E_criteria = 1e-5)
smalldna <- VhgSubsetHittable(smalldna,ViralRefSeq_E_criteria = 1e-5)
largedna <- VhgSubsetHittable(largedna,ViralRefSeq_E_criteria = 1e-5)

rna <- VhgPreprocessTaxa(rna,taxa_rank = "Family")
smalldna <- VhgPreprocessTaxa(smalldna,taxa_rank = "Family")
largedna <- VhgPreprocessTaxa(largedna,taxa_rank = "Family")


rna_subject <- VhgGetSubject(rna,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE)
smalldna_subject <- VhgGetSubject(smalldna,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE)
largedna_subject <- VhgGetSubject(largedna,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE)

rna_gt <- VhgTabularRasa(rna_subject,title = "Taubert - Found subjects for RNA viruses",names_ = c("Viral reference taxonomy",
                                                                               "Viral reference subject",
                                                                               "Count"),title_colour = "black")
smalldna_gt <-VhgTabularRasa(smalldna_subject,title = "Taubert - Found subjects for small DNA viruses",names_ = c("Viral reference taxonomy",
                                                                                    "Viral reference subject",
                                                                                    "Count"),title_colour = "black")
largedna_gt <-VhgTabularRasa(largedna_subject,title = "Taubert - Found subjects for large DNA viruses",names_ = c("Viral reference taxonomy",
                                                                                    "Viral reference subject",
                                                                                    "Count"),title_colour = "black")

rna_gt<- rna_gt%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

smalldna_gt <- smalldna_gt%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

largedna_gt <- largedna_gt%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

gtsave(smalldna_gt,filename = "smalldna_gt.png",path = "output/TaubertDatacombined/plots/smalldna/Gatherer/")
gtsave(largedna_gt,filename = "largedna_gt.png",path = "output/TaubertDatacombined/plots/largedna/Gatherer/")
gtsave(rna_gt,filename = "rna_gt.png",path = "output/TaubertDatacombined/plots/rnavirus/Gatherer/")



combined_subjects <- rbind(rna_subject,smalldna_subject,largedna_subject)


cgt <-  VhgTabularRasa(combined_subjects,title = "Taubert - Found subjects",names_ = c("Viral reference taxonomy",
                                                                                                 "Viral reference subject",
                                                                                                 "Count"),title_colour = "black")


cgt<- cgt%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

cgt <- cgt |>
  # Add row groups
  tab_row_group(
    label = "RNA Viruses",
    rows = 1:5
  ) |>
  tab_row_group(
    label = "Small DNA Viruses",
    rows = 6
  ) |>
  tab_row_group(
    label = "Large DNA Viruses",
    rows = 7:11
  ) |>
  # Apply custom style for row group labels
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_row_groups(groups = c("RNA Viruses", "Small DNA Viruses", "Large DNA Viruses"))
  )

gtsave(cgt,filename = "cgt_gt.png",path = "output/TaubertDatacombined/")