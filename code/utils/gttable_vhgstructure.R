library(gt)
library(Virusparies)
library(dplyr)

test <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Bunya/virushunter.tsv")
test2 <- ImportVirusTable("data/Florian_data/RNA/virushunter.tsv")

test <- head(test,5)
test2 <- head(test2,5)


test$run_id <- test2$run_id

# vg_file <- test[,c(16,1:15)] for hunter
vg_file <- test[,c(14,1:13)] # for gatherer

# Function to get example values
get_examples <- function(column) {
  if (is.factor(column) || is.character(column)) {
    return(paste0(head(column, 3), collapse = ", "))
  } else if (is.numeric(column) || is.integer(column)) {
    # Use scientific notation format if necessary
    formatted <- formatC(head(column, 3), format = "e", digits = 2)
    return(paste0(formatted, collapse = ", "))
  } else {
    return(paste0(head(column, 3), collapse = ", "))
  }
}

# Create a data frame with column names, data types, and example values
column_info <- data.frame(
  Column_Name = names(vg_file),
  Data_Type = sapply(vg_file, function(col) {
    class(col) %>%
      gsub("character", "chr", .) %>%
      gsub("numeric", "dbl", .) %>%
      gsub("integer", "int", .) %>%
      gsub("factor", "fct", .) %>%
      gsub("Date", "date", .)
  }),
  Examples = sapply(vg_file, get_examples),
  stringsAsFactors = FALSE
)

# Create the gt table
gt_table <- column_info %>%
  gt() %>%
  tab_header(
    title = "Column Overview"
  ) %>%
  cols_label(
    Column_Name = "Column Name",
    Data_Type = "Data Type",
    Examples = "Example Values"
  ) %>%
  fmt_markdown(
    columns = c(Examples)
  ) %>%
  cols_align(
    align = "left",
    columns = c(Column_Name, Data_Type, Examples)
  ) %>%
  fmt_number(
    columns = c(Examples),
    decimals = 2
  )

# Print the gt table
print(gt_table)

path <- "misc/vhg_structure/"
ExportVirusGt(gtable=gt_table,filename="vg_structure.docx",path = path)