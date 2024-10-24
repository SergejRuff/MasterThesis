rm(list=ls())

library(Virusparies)
library(ggplot2)
library(tidyverse)


m1 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virusgatherer-cap3.tsv")
m2 <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virusgatherer-cap3.tsv")

m1 <- m1[c(4,31,76),]

m1 <- head(m1, 3)
m2 <- head(m2, 3)

m1 <- m2 %>%
  select(run_id) %>%  # Select the run_id column from m2
  bind_cols(m1) 

# Step 1: Extract the first 3 rows
first_three <- head(m1, 3)

# Create summary_data as before
summary_data <- data.frame(
  Column = names(first_three),
  Type = sapply(first_three, function(x) {
    dtype <- class(x)[1]
    case_when(
      dtype == "numeric" ~ "<num>",
      dtype == "integer" ~ "<int>",
      dtype == "character" ~ "<chr>",
      dtype == "factor" ~ "<fctr>",
      dtype == "Date" ~ "<date>",
      dtype == "POSIXct" | dtype == "POSIXt" ~ "<dttm>",
      TRUE ~ paste0("<", dtype, ">")
    )
  }),
  Example_1 = sapply(first_three, function(x) as.character(x[1])),
  Example_2 = sapply(first_three, function(x) as.character(x[2])),
  Example_3 = sapply(first_three, function(x) as.character(x[3]))
)

# Combine the examples into a single string for each column
summary_data <- summary_data %>%
  mutate(Examples = paste(Example_1, Example_2, Example_3, sep = ", ")) %>%
  select(Column, Type, Examples)

# Convert the summary to a gt table and apply styling
summary_data <- summary_data %>%
  gt() %>%
  tab_header(title = "Gatherer Hittable Structure") %>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  ) %>%
  # Highlight the first row
  tab_style(
    style = list(
      cell_fill(color = "#FFD700")  # Background color

    ),
    locations = cells_body(
      rows = 1
    )
  ) %>%
  # Highlight rows 2 to 4 with a different color
  tab_style(
    style = list(
      cell_fill(color = "#CEA2FD")  # Background color for rows 2 to 4

    ),
    locations = cells_body(
      rows = 2:6
    )
  )


gtsave(filename = "gatherer_structure.png",data = summary_data,path = "data/")