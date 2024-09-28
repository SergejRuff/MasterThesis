rm(list = ls())

library(Virusparies)
library(ggplot2)

taubert <- ImportVirusTable("data/hittables_taubert/largedna/combined_virusgatherer-cap3.tsv")
taubert_gatherer <- ImportVirusTable("data/hittables_taubert/")

export_path <- "output/TaubertDatacombined/plots/largedna/"

taubert <- VhgPreprocessTaxa(taubert,taxa_rank = "Family")


# evalbox <- VhgBoxplot(taubert,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "RNA viruses",axis_title_size = 12,ytext_size = 12,
#                       xtext_size = 11,title = "Taubert data (Gatherer)\nViral group E-values")

sra <- VhgRunsBarplot(taubert,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                      theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",remove_group_labels = FALSE,ytext_size = 10,
                      title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "Large DNA viruses")


identity <- VhgBoxplot(taubert,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "",
                       axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

contiglen <- VhgBoxplot(taubert,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                      contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = FALSE,ytext_size = 12,
                      title = NULL,xtext_size = 11,legend_text_size = 12,xlabel = "Large DNA viruses")


# evalbox$boxp <- evalbox$boxp +
#   scale_y_continuous(limits = c(-1, 30))




plot_list <- list(sra$plot,identity$boxp,contiglen$boxp)

# Adjust y-axis title margin
plot_list <- lapply(plot_list, function(p) {
  p + theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
})


ExportVirusPlot("taubert_largedna_grid.png",plot = plot_list,
                width = 12,height = 8,units = "in",limitsize = FALSE,ncol = 2,nrow = 2,greedy = TRUE,align="h",
                labels = c("A","B","C"),path = export_path)