# Load necessary library
library(dplyr)
library(Virusparies)

# Define the path to the root directory
root_dir <- "data/hittables_taubert/smalldna/taubert"

# Recursively find all virushunter.tsv files
tsv_files <- list.files(path = root_dir, pattern = "virusgatherer-cap3.tsv", recursive = TRUE, full.names = TRUE)
# tsv_files <- list.files(path = root_dir, pattern = "virushunter.tsv", recursive = TRUE, full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()


# Loop through each file and read it into a data frame, then add to the list
for (file in tsv_files) {
  df <- ImportVirusTable(file)
  
  # Extract the folder name from the file path
  folder_name <- basename(dirname(dirname(file)))
  
  # Add the folder name as a new column
  df <- df %>%
    mutate(folder = folder_name)
  
  data_list <- append(data_list, list(df))
}

# Combine all data frames into one using row binding
combined_df <- do.call(rbind,data_list)

# Save the combined data frame to a new file
write.table(combined_df, file = "data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(combined_df, file = "data/hittables_taubert/smalldna/combined_virushunter.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Print a message indicating completion
cat("Combined file saved as combined_virushunter.tsv")
