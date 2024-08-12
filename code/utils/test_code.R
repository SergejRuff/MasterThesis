


#' pivort ICVT Data
#'
#' @param taxa_rank taxa_rank
#'
#' @return pivotted ICTV_data
#'
#' @importFrom tidyr pivot_longer
#' @noRd
format_ICTV <- function(taxa_rank){
  
  
  
  return(ICTV_data %>%
           select(.data$Phylum:.data$Subgenus) %>%
           pivot_longer(.data$Subphylum:.data$Subgenus, names_to = "level", values_to = "name") %>%
           filter(str_detect(.data$name, paste0("\\w+", taxa_rank), negate = TRUE)) %>%
           na.omit() %>%
           distinct())
  
  
  
}


#' define taxa string
#'
#' @param taxa_rank a character indicating the taxa rank
#'
#' @return a character with the taxa suffix
#'
#' @noRd
taxonomy_rank_hierarchy <- function(taxa_rank) {
  
  # Define the list with values
  taxa_list <- list(
    "viricotina",
    "viricetes",
    "viricetidae",
    "virales",
    "virineae",
    "viridae",
    "virinae",
    "virus"
  )
  
  # Define the names for each element
  taxa_names <- c(
    "Subphylum",
    "Class",
    "Subclass",
    "Order",
    "Suborder",
    "Family",
    "Subfamily",
    "Genus"
  )
  
  # Assign names to the list elements using setNames()
  taxa_list <- setNames(taxa_list, taxa_names)
  
  valid_ranks <- names(taxa_list)
  
  if (!(taxa_rank %in% valid_ranks)) {
    stop("Error: Invalid taxa rank provided. Please provide one of: Subphylum, Class, Subclass, Order, Suborder, Family, Subfamily, Genus")
  }
  
  return(taxa_list[[taxa_rank]])
}



#' @title VhgPreprocessTaxa: preprocess ViralRefSeq_taxonomy elements
#'
#' @details
#' Process the `ViralRefSeq_taxonomy` column.
#'
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param taxa_rank (optional): Specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#'
#'
#' @details
#' Besides `best_query`, the user can utilize the `ViralRefSeq_taxonomy` column as `x_column` or `groupby` in plots.
#' This column needs preprocessing because it is too long and has too many unique elements for effective grouping.
#' The element containing the taxa suffix specified by the `taxa_rank` argument is used. NA values are replaced by "unclassified".
#'
#' This function is used internally by every function that can use the `ViralRefSeq_taxonomy` column as input.
#' The user can also apply this function independently to process the taxonomy column and filter for the selected taxa rank in their data.
#'
#'
#' @return file with preprocessed ViralRefSeq_taxonomy elements
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' file_filtered <- VhgPreprocessTaxa(file,"Family")
#'
#' print("ViralRefSeq_taxonomy before processing:\n")
#' print(head(file$ViralRefSeq_taxonomy,5))
#'
#' print("ViralRefSeq_taxonomy after processing:\n")
#' print(head(file_filtered$ViralRefSeq_taxonomy,5))
#'
#'
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom tidyr unnest pivot_longer
#' @importFrom stringr  str_extract str_remove_all str_detect
#' @export
VhgPreprocessTaxa <- function(file,taxa_rank) {
  
  
  taxa_rank <- taxonomy_rank_hierarchy(taxa_rank)
  
  
  
  ictv_formatted <- format_ICTV(taxa_rank)
  
  
  taxon_filter <- paste(unique(ictv_formatted$name), collafilpse = "|")
  
  
  file <- file %>%
    mutate(
      ViralRefSeq_taxonomy = str_remove_all(ViralRefSeq_taxonomy, "taxid:\\d+\\||\\w+\\s\\w+\\|"),
      name = str_extract(ViralRefSeq_taxonomy, taxon_filter),
      ViralRefSeq_taxonomy = str_extract(ViralRefSeq_taxonomy, paste0("\\w+", taxa_rank))
    ) %>%
    left_join(ictv_formatted, join_by("name" == "name")) %>%
    mutate(
      ViralRefSeq_taxonomy = case_when(
        is.na(ViralRefSeq_taxonomy) & is.na(.data$Phylum) ~ "unclassified",
        is.na(ViralRefSeq_taxonomy) ~ paste("unclassified", .data$Phylum),
        .default = ViralRefSeq_taxonomy
      )
    ) %>% select(-c(name:level))%>%
    mutate(
      ViralRefSeq_taxonomy = if_else(ViralRefSeq_taxonomy == "unclassified unclassified", "unclassified", ViralRefSeq_taxonomy),
      ViralRefSeq_taxonomy = if_else(ViralRefSeq_taxonomy == "unclassified NA", "unclassified", ViralRefSeq_taxonomy)
    )
  
  
  
  return(file)
}



################################################






## test 

# base_str_extract <- function(txt, pattern, perl = TRUE, ...) {
#   x <- regexpr(pattern, txt, perl = perl, ...)
#   x[which(x == -1)] <- NA
#   x[which(x != -1)] <- regmatches(txt, x)
#   as.character(x)
# }

##################################################################################
##################################################################################
##################################################################################
##################################################################################

