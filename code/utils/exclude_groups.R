rm(list = ls())
library(Virusparies)
library(readxl)
library(stringr)
library(tidyverse)


vh_file <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/combined_ga.tsv")

ICTV_data <- read_excel("/media/sergej/My Book/Rskripte_f-r_masterarbeit-main/data//ICTV_Master_Species_List_2023_MSL39.v2.xlsx",sheet = 2)


vh_file <- VhgAddPhylum(vh_file,"ViralRefSeq_taxonomy")


valid_phyla <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
  "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")




vh_file <- vh_file %>%
  mutate(ViralRefSeq_taxonomy = if_else(
    grepl(paste0(valid_phyla, collapse = "|"), Phylum, ignore.case = TRUE),
    ViralRefSeq_taxonomy,
    "Non-RNA-virus"
  ))

# Print the updated data frame
print(vh_file)

# Filter out entries that are not in the valid phyla list, excluding "unclassified"
filtered_phyta_colors <- phyla_colors[names(phyla_colors) %in% c(valid_phyla, "unclassified")]

# Add the "Non-RNA-virus" entry with the color black
filtered_phyta_colors["Non-RNA-virus"] <- "#000000"

# Print the result
print(filtered_phyta_colors)

remove_non_group <- function(file,groupby,chosen_group,label_vector){
  
  valid_phyla_rna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                    "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")
  
  
  valid_phyla_smalldna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                        "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")
  
  valid_phyla_largedna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                        "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")
  
  chosen_list <- switch(chosen_group,
                        "rna" = valid_phyla_rna,
                        "smalldna" = valid_phyla_smalldna,
                        "largedna" = valid_phyla_largedna,
                        stop("Invalid chosen_group value. Use 'rna', 'smalldna', or 'largedna'."))
  
  change_label <- switch(chosen_group,
                         "rna" = "Non-RNA-virus",
                         "smalldna" = "Non-Small-DNA-Virus",
                         "largedna" = "Non-Large-DNA-Virus")
  
  file <- VhgAddPhylum(file,"ViralRefSeq_taxonomy")
  
  file <- file %>%
    mutate(ViralRefSeq_taxonomy = if_else(
      grepl(paste0(valid_phyla, collapse = "|"), Phylum, ignore.case = TRUE),
      .data[[groupby]],
      change_label
    ))
  
  # Filter out entries that are not in the valid phyla list, excluding "unclassified"
  label_vector <- label_vector[names(label_vector) %in% c(valid_phyla, "unclassified")]
  
  # Add the "Non-RNA-virus" entry with the color black
  label_vector[change_label] <- "#000000"
  
  return(list(file =file,label=label_vector))
  
  
}
