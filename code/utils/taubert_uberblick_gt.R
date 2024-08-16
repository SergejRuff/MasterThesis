rm(list = ls())

# devtools::install_github("hadley/emo")

library(readxl)
library(gt)
library(emo)
library(tidyverse)

taubert <- read_excel("data/hittables_taubert/taubert_daten.xlsx")

taubert <- filter(taubert, row_number() <= n()-1)

taubert_gt <- taubert %>% 
  gt()%>% 
  cols_add(
    comp1 = rep(emo::ji("white heavy check mark"),nrow(taubert)),
    comp2 = rep(emo::ji("white heavy check mark"),nrow(taubert)),
    comp3 = rep(emo::ji("white heavy check mark"),nrow(taubert))
  ) %>% 
  tab_spanner(
    label = "Completed",
    columns = c(comp1,comp2,comp3),
    id = 'spannerA'
  ) %>% 
  cols_label(
    Folder = "Folder",
    Number_of_FASTQfiles = "Number of Samples",
    comp1 = "Small DNA viruses",
    comp2 = "Large DNA viruses",
    comp3 = "RNA viruses"
    
    
  ) %>% 
  tab_header(
    title = md("***Taubert Data***"),
    subtitle = md("Total Number of FASTQ Files: **323**")
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

print(taubert_gt)

gtsave(taubert_gt,"taubert_uberblick.png",path = "data/hittables_taubert/")
