rm(list = ls())

library(Virusparies)
library(ggplot2)

taubert_largedna <- ImportVirusTable("data/hittables_taubert/largedna/combined_virusgatherer-cap3.tsv")

taubert_smalldna <- ImportVirusTable("data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv")
taubert_rna <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virusgatherer-cap3.tsv")

export_path <- "output/TaubertDatacombined/plots/"

taubert_largedna <- VhgPreprocessTaxa(taubert_largedna,taxa_rank = "Family")
taubert_smalldna <- VhgPreprocessTaxa(taubert_smalldna,taxa_rank = "Family")

taubert_rna <- VhgPreprocessTaxa(taubert_rna,taxa_rank = "Family")





evalbox_large <- VhgBoxplot(taubert_largedna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "Large DNA viruses",axis_title_size = 12,ytext_size = 12,
                      xtext_size = 11,title = "Taubert data (Gatherer)\nViral group E-values")

evalbox_small <- VhgBoxplot(taubert_smalldna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "Small DNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = "Taubert data (Gatherer)\nViral group E-values")

evalbox_rna <- VhgBoxplot(taubert_rna ,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "RNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = "Taubert data (Gatherer)\nViral group E-values")




evalbox_large$boxp <- evalbox_large$boxp +
  scale_y_continuous(limits = c(-1, 9))



evalbox_small$boxp <- evalbox_small$boxp +
  scale_y_continuous(limits = c(0, 8))

evalbox_rna$boxp <- evalbox_rna$boxp +
  scale_y_continuous(limits = c(-1, 30))



plot_list <- list(evalbox_large$boxp,evalbox_small$boxp,evalbox_rna$boxp)

# Adjust y-axis title margin
plot_list <- lapply(plot_list, function(p) {
  p + theme(axis.title.y = element_text(margin = margin(r = 15, unit = "pt")))
})


ExportVirusPlot("taubert_evalue_grid.png",plot = plot_list,
                width = 8,height = 11,units = "in",limitsize = FALSE,ncol = 1,nrow = 3,greedy = TRUE,align="h",
                labels = c("A","B","C"),path = export_path)