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

vh_file %>% ggplot()+
  geom_histogram(mapping = aes(ViralRefSeq_ident,fill = ViralRefSeq_taxonomy),bins = 100,position = "identity",alpha=0.7)


vh_file %>% ggplot()+
  geom_density_ridges(mapping = aes(x = ViralRefSeq_ident,y=ViralRefSeq_taxonomy,fill = ViralRefSeq_taxonomy),scale = 0.9,alpha=0.2)

VhgHist <- function(file,
                    plot_type=NULL,
                    ){
  
  p <- file %>% ggplot()
  
  if (is.null(plot_type) || plot_type == "simple_histogram") {
    # Default case: Simple histogram without color differentiation
    p <- p + 
      geom_histogram(mapping = aes(ViralRefSeq_ident), bins = 100, alpha = 0.7)
  } else if (plot_type == "taxa_colored_histogram") {
    # Taxa-colored histogram: Bars are filled based on ViralRefSeq_taxonomy
    p <- p + 
      geom_histogram(mapping = aes(ViralRefSeq_ident, fill = ViralRefSeq_taxonomy), 
                     bins = 100, position = "identity", alpha = 0.7)
  } else if (plot_type == "ridge_plot") {
    # Ridge plot: Visualizes the distribution of ViralRefSeq_ident across different taxa
    p <- p + 
      geom_density_ridges(mapping = aes(x = ViralRefSeq_ident, y = ViralRefSeq_taxonomy, 
                                        fill = ViralRefSeq_taxonomy), 
                          scale = 0.9, alpha = 0.2)
  } else {
    stop("Invalid plot_type. Choose either 'simple_histogram', 'taxa_colored_histogram', or 'ridge_plot'.")
  }
  

}