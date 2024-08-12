rm(list = ls())

library(data.table)
library(readxl)
library(Virusparies)
library(profvis)
library(tidyverse)
library(microbenchmark)
#source("code/utils/multi_process_taxaprocess.R")
#source("code/utils/test_code.R")


ICTV_data <- read_excel("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/ICTV_data/ICTV_Master_Species_List_2023_MSL39.v2.xlsx",sheet = 2)


filepath <- "data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv"

file <- ImportVirusTable(filepath)

results <- microbenchmark(
  SingleThreaded = VhgPreprocessTaxa(file, taxa_rank = "Family"),
  Parallel4Cores = VhgPreprocessTaxa2(file, taxa_rank = "Family", num_cores = 4),
  times = 10
)

for (i in 1:7){
  
  
  cat("\n cores:", i, "\n")
  
  res <-bench::mark(
    Parallel4Cores = VhgPreprocessTaxa2(file, taxa_rank = "Family", num_cores = i),memory = FALSE)
  
  print(res)
}


bench::mark(
  SingleThreaded = VhgPreprocessTaxa(file, taxa_rank = "Family"))

print(results)

# test current package version.
profvis::profvis(VhgPreprocessTaxa(file, taxa_rank = "Family"))

profvis::profvis(VhgPreprocessTaxa2(file, taxa_rank = "Family", num_cores = 7))

