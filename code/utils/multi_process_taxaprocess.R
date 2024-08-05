# Load necessary packages
library(dplyr)
library(stringr)
library(parallel)

# Define the function to process each chunk
process_chunk <- function(chunk, ictv_formatted, taxa_rank) {
  taxon_filter <- paste(unique(ictv_formatted$name), collapse = "|")
  
  chunk_processed <- chunk %>%
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
    ) %>% select(-c(name:level)) %>%
    mutate(
      ViralRefSeq_taxonomy = if_else(ViralRefSeq_taxonomy == "unclassified unclassified", "unclassified", ViralRefSeq_taxonomy),
      ViralRefSeq_taxonomy = if_else(ViralRefSeq_taxonomy == "unclassified NA", "unclassified", ViralRefSeq_taxonomy)
    )
  
  return(chunk_processed)
}

# Define the main function
VhgPreprocessTaxa2 <- function(file, taxa_rank, num_cores = 1) {
  taxa_rank <- taxonomy_rank_hierarchy(taxa_rank)
  ictv_formatted <- format_ICTV(taxa_rank)
  
  # Validate num_cores
  if (num_cores < 1) {
    stop("Number of cores must be at least 1")
  }
  
  # Use default number of cores if specified is greater than available
  num_cores <- min(num_cores, detectCores() - 1)
  
  # Split the data into chunks
  num_chunks <- num_cores
  chunks <- split(file, rep(1:num_chunks, length.out = nrow(file)))
  
  # Process chunks in parallel
  processed_chunks <- mclapply(chunks, process_chunk, ictv_formatted = ictv_formatted, taxa_rank = taxa_rank, mc.cores = num_cores)
  
  # Combine processed chunks
  file_processed <- bind_rows(processed_chunks)
  
  return(file_processed)
}
