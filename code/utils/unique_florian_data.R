rm(list=ls())

library(readxl)
library(Virusparies)
library(tidyverse)
library(gt)
library(emo)

florian <- read_excel("data/Florian_data/taubert_daten_überblick.xlsx",col_names = FALSE)

test <- gsub("_R[12].fastq.gz","",florian$...1)

test <- gsub("_R[12]_001.fastq.gz","",test)

test <- gsub("_[12].fq.gz","",test)

test <- gsub(".fastq.gz","",test)
test <- gsub("_L00[12]","",test)


test <- unique(test)

print(test)

print(length(test))


library(tidyverse)

fl <- florian %>% arrange(...1)

rna <- ImportVirusTable("data/Florian_data/RNA/virushunter.tsv")
smalldna <- ImportVirusTable("data/Florian_data/smalldna//virushunter.tsv")
largedna <- ImportVirusTable("data/Florian_data/largedna//virushunter.tsv")

rna <- VhgSubsetHittable(rna,ViralRefSeq_E_criteria = 1e-5)
smalldna <- VhgSubsetHittable(smalldna,ViralRefSeq_E_criteria = 1e-5)
largedna <- VhgSubsetHittable(largedna,ViralRefSeq_E_criteria = 1e-5)

rna$run_id <- gsub("_R[12]", "", rna$run_id)
rna$run_id <- gsub("_R[12]_001", "", rna$run_id)
rna$run_id <- gsub("_[12]", "", rna$run_id)
#rna$run_id <- gsub(".fastq.gz", "", rna$run_id)
#rna$run_id <- gsub("_R[12]_001", "", rna$run_id)
#rna$run_id <- gsub("_R[12]", "", rna$run_id)
#rna$run_id <- gsub("_L00[12]_001", "", rna$run_id)

smalldna$run_id <- gsub("_R[12]", "", smalldna$run_id)
smalldna$run_id <- gsub("_R[12]_001", "", smalldna$run_id)
smalldna$run_id <- gsub("_[12]", "", smalldna$run_id)
#smalldna$run_id <- gsub("_L00[12]_001", "", smalldna$run_id)


largedna$run_id <- gsub("_R[12]", "", largedna$run_id)
largedna$run_id <- gsub("_R[12]_001", "", largedna$run_id)
largedna$run_id <- gsub("_[12]", "", largedna$run_id)
#largedna$run_id <- gsub("_L00[12]_001", "", largedna$run_id)

largedna <- VhgPreprocessTaxa(largedna,"Family")
rna <- VhgPreprocessTaxa(rna,"Family")
smalldna <- VhgPreprocessTaxa(smalldna,"Family")

# Summarize the number of unique ViralRefSeq_taxonomy for each run_id
summary_result_rna <- aggregate(ViralRefSeq_taxonomy ~ run_id, data = rna, function(x) length(unique(x)))

# Print the summary result
print(summary_result_rna)

# Summarize the number of unique ViralRefSeq_taxonomy for each run_id
summary_result_large <- aggregate(ViralRefSeq_taxonomy ~ run_id, data = largedna, function(x) length(unique(x)))

# Print the summary result
print(summary_result_large)

# Summarize the number of unique ViralRefSeq_taxonomy for each run_id
summary_result_small <- aggregate(ViralRefSeq_taxonomy ~ run_id, data = smalldna, function(x) length(unique(x)))

# Print the summary result
print(summary_result_small)

singleend <- c("Tx00205.S5388.1.fastq.gz","Tx00232.S5388.1.fastq.gz")

df <- data.frame(test)

# Merge all data frames together
combined_df <- Reduce(function(x, y) merge(x, y, by.x = "test", by.y = "run_id", all.x = TRUE, all.y = TRUE), 
                      list(df, summary_result_small, summary_result_large, summary_result_rna))

# Replace all NA values with 0
combined_df[is.na(combined_df)] <- 0

combined_df <- combined_df[-16,]

# Print the combined data frame
print(combined_df)


combined_df <- combined_df %>%
  mutate(pairend = c(rep(emo::ji("white heavy check mark"), 12),
                     rep(emo::ji("cross mark"), 2),
                     rep(emo::ji("white heavy check mark"), 3)))
                  

# Reorder the columns to make 'pairend' the second column
combined_df <- combined_df %>%
  select(test, pairend, everything())

combined_df <- combined_df %>% 
  gt()%>% 
  tab_spanner(
    label = "Unique HMM profiles",
    columns = 3:5,
    id = 'spannerA'
  ) %>% 
  cols_label(
    test = "Sample",
    pairend = "Pair-end?",
    ViralRefSeq_taxonomy.x = "Small DNA viruses",
    ViralRefSeq_taxonomy.y = "Large DNA viruses",
    ViralRefSeq_taxonomy = "RNA viruses"
  ) %>% 
  tab_header(
    title = md("***Florian Data***"),
    subtitle = md("Total number of files: **32**")

  ) %>% 
  opt_align_table_header(align = "left") %>% 
  tab_options(
    # These were the ones we applied in the first chapter
    data_row.padding = px(2),
    summary_row.padding = px(3), # A bit more padding for summaries
    row_group.padding = px(4)    # And even more for our groups
  ) %>% opt_stylize(style = 4) %>% 
  tab_options(
    heading.title.font.size = px(20)
  )


rna_p <- VhgRunsBarplot(rna,theme_choice = "linedraw_dotted",plot_text = 1,
               title = "Florian data: RNA viruses",title_size = 14,axis_title_size = 12,xtext_size = 12)

sm <- VhgRunsBarplot(smalldna,theme_choice = "linedraw_dotted",plot_text = 1,
               title = "Florian data: Small DNA viruses",title_size = 14,axis_title_size = 12,xtext_size = 12)

la <- VhgRunsBarplot(largedna,theme_choice = "linedraw_dotted",plot_text = 1,
                     title = "Florian data: Large DNA viruses",title_size = 14,axis_title_size = 12,xtext_size = 12)


path <- "output/Florian/"


ExportVirusPlot(plot = rna_p$plot,file_name = "flo_rna.png",path=path,
                width = 7,height = 7,units = "in")

ExportVirusPlot(plot = sm$plot,file_name = "flo_sm.png",path=path,
                width = 7,height = 7,units = "in")

ExportVirusPlot(plot = la$plot,file_name = "flo_la.png",path=path,
                width = 7,height = 7,units = "in")
