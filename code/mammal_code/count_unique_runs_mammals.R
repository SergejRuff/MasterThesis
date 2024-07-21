# Define the directory containing the .txt files
directory_path <- "D:/git_clone/RscriptsMasterthesis/data/RNAvirus_Mammals_newJan2023/sra_ids"

# List all .txt files in the directory
file_names <- list.files(path = directory_path, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty data frame to hold all data
all_data <- data.frame()

# Read each file and combine them
for (file in file_names) {
  file_data <- read_delim(file, delim = "\t", col_names = FALSE)  # Assuming tab-separated values, adjust if necessary
  all_data <- bind_rows(all_data, file_data)
}

# Assuming SRA accession numbers are in the first column
# Rename the column for clarity
colnames(all_data) <- "SRA_Run"

# Get unique SRA runs
unique_sra_runs <- all_data %>% distinct(SRA_Run)

# Count the number of unique SRA runs
num_unique_sra_runs <- nrow(unique_sra_runs)

# Print the number of unique SRA runs
print(num_unique_sra_runs)