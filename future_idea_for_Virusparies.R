# Idea: Add VhgHist to Virusparies with ability to plot
# 1 histgram for all
# hist by group
# seperate into ridge plot


rm(list = ls())

library(Virusparies)
library(tidyverse)
library(ggridges)

vh_file <- ImportVirusTable("/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv")

vh_file <- VhgPreprocessTaxa(vh_file,taxa_rank = "Family")

vh_file <- VhgSubsetHittable(vh_file,group_column = "ViralRefSeq_taxonomy",c("Flaviviridae","Myriaviridae",
                                                                             "Togaviridae","Bromoviridae"))

vh_file %>% ggplot(mapping = aes(ViralRefSeq_ident,fill = ViralRefSeq_taxonomy))+
  geom_histogram(bins = 100,position = "identity",alpha=0.7)


vh_file %>% ggplot(mapping = aes(x = ViralRefSeq_ident,y=ViralRefSeq_taxonomy,fill = ViralRefSeq_taxonomy))+
  geom_density_ridges(scale = 5,alpha=0.2)