rm(list = ls())

library(Virusparies)
library(ggplot2)

taubert_small <- ImportVirusTable("data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv")

taubert_large <- ImportVirusTable("data/hittables_taubert/largedna/combined_virusgatherer-cap3.tsv")

taubert_rna <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virusgatherer-cap3.tsv")

export_path <- "output/TaubertDatacombined/plots/"

taubert_small <- VhgPreprocessTaxa(taubert_small,taxa_rank = "Family")

taubert_large  <- VhgPreprocessTaxa(taubert_large ,taxa_rank = "Family")

taubert_rna  <- VhgPreprocessTaxa(taubert_rna ,taxa_rank = "Family")


#### SmallDNA

taubert_small_p1 <- VhgRunsBarplot(taubert_small,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                   theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",remove_group_labels = FALSE,ytext_size = 9,
                                   title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "Small DNA viruses",subtitle_size = 11,axis_title_size = 10)


taubert_small_p2 <- VhgBoxplot(taubert_small,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "none",xlabel = "",
                               axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

taubert_small_p3 <- VhgBoxplot(taubert_small,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                               contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = TRUE,ytext_size = 12,
                               title = NULL,xtext_size = 11,legend_text_size = 12,xlabel = "")



#### Large DNA

taubert_large_p1 <- VhgRunsBarplot(taubert_large,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                   theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",remove_group_labels = FALSE,ytext_size = 9,
                                   title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "Large DNA viruses",subtitle_size = 11,axis_title_size = 10)


taubert_large_p2 <- VhgBoxplot(taubert_large,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "none",xlabel = "",
                               axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

taubert_large_p3 <- VhgBoxplot(taubert_large,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                               contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = TRUE,
                               title = NULL,xtext_size = 11,legend_text_size = 12,xlabel ="")


##### RNAvirus

taubert_rna_p1 <- VhgRunsBarplot(taubert_rna,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                 theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",remove_group_labels = FALSE,ytext_size = 9,
                                 title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "RNA viruses",subtitle_size = 11,axis_title_size = 10)


taubert_rna_p2 <- VhgBoxplot(taubert_rna,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "none",xlabel = "",
                             axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

taubert_rna_p3 <- VhgBoxplot(taubert_rna,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                             contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = TRUE,ytext_size = 12,
                             title = NULL,xtext_size = 11,legend_text_size = 12,xlabel ="")



################################# Export

plot_list <- list(taubert_small_p1$plot,taubert_small_p2$boxp,taubert_small_p3$boxp,
                  taubert_large_p1$plot,taubert_large_p2$boxp,taubert_large_p3$boxp,
                  taubert_rna_p1$plot,taubert_rna_p2$boxp,taubert_rna_p3$boxp)

# Adjust y-axis title margin
plot_list <- lapply(plot_list, function(p) {
  p + theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
})


ExportVirusPlot("grid_9x9_taubert.png",plot = plot_list,
                width = 13,height = 9,units = "in",limitsize = FALSE,ncol = 3,nrow = 3,greedy = TRUE,align="h",
                labels = c("A","B","C","D","E","F","G","H","I"),path = export_path,rel_widths = c(1, 0.6, 1))