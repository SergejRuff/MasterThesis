rm(list = ls())

library(Virusparies)
library(ggplot2)
library(patchwork)
library(gt)

taubert_largedna <- ImportVirusTable("data/hittables_taubert/largedna/combined_virushunter.tsv")

taubert_smalldna <- ImportVirusTable("data/hittables_taubert/smalldna/combined_virushunter.tsv")
taubert_rna <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virushunter.tsv")

export_path <- "output/TaubertDatacombined/plots/"

taubert_largedna <- VhgPreprocessTaxa(taubert_largedna,taxa_rank = "Family")
taubert_smalldna <- VhgPreprocessTaxa(taubert_smalldna,taxa_rank = "Family")

taubert_rna <- VhgPreprocessTaxa(taubert_rna,taxa_rank = "Family")



evalbox_small <- VhgBoxplot(taubert_smalldna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "Small DNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = NULL)

evalbox_large <- VhgBoxplot(taubert_largedna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "Large DNA viruses",axis_title_size = 12,ytext_size = 12,
                      xtext_size = 11,title = NULL,subtitle = NULL)


evalbox_rna <- VhgBoxplot(taubert_rna ,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "RNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = NULL,subtitle = NULL)




evalbox_large$boxp <- evalbox_large$boxp +
  scale_y_continuous(limits = c(-0.5, 10))



evalbox_small$boxp <- evalbox_small$boxp +
  scale_y_continuous(limits = c(-0.5, 10))

evalbox_rna$boxp <- evalbox_rna$boxp +
  scale_y_continuous(limits = c(-0.5, 10))


eval_plot_co <- evalbox_small$boxp/evalbox_large$boxp/evalbox_rna$boxp+
  plot_annotation(tag_levels = "A")

ggsave("combined_eval_hunter_boxplots.png",plot = eval_plot_co,
       width = 10,height = 15,units = "in",limitsize = FALSE,path = export_path )

plot_list <- list(evalbox_large$boxp,evalbox_small$boxp,evalbox_rna$boxp)

# Adjust y-axis title margin
plot_list <- lapply(plot_list, function(p) {
  p + theme(axis.title.y = element_text(margin = margin(r = 15, unit = "pt")))
})


ExportVirusPlot("taubert_evalue_grid.png",plot = plot_list,
                width = 8,height = 11,units = "in",limitsize = FALSE,ncol = 1,nrow = 3,greedy = TRUE,align="h",
                labels = c("A","B","C"),path = export_path)



######################################################################################


rm(list = ls())

library(Virusparies)
library(ggplot2)

taubert_largedna <- ImportVirusTable("data/hittables_taubert/largedna/combined_virushunter.tsv")

taubert_smalldna <- ImportVirusTable("data/hittables_taubert/smalldna/combined_virushunter.tsv")
taubert_rna <- ImportVirusTable("data/hittables_taubert/rnavirus/combined_virushunter.tsv")

export_path <- "output/TaubertDatacombined/plots/"

taubert_largedna <- VhgPreprocessTaxa(taubert_largedna,taxa_rank = "Family")
taubert_smalldna <- VhgPreprocessTaxa(taubert_smalldna,taxa_rank = "Family")

taubert_rna <- VhgPreprocessTaxa(taubert_rna,taxa_rank = "Family")





evalbox_large <- VhgBoxplot(taubert_largedna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = "phylum_median",legend_position = "right",xlabel = "Large DNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = NULL)

evalbox_small <- VhgBoxplot(taubert_smalldna,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = "phylum_median",legend_position = "right",xlabel = "Small DNA viruses",axis_title_size = 12,ytext_size = 12,
                            xtext_size = 11,title = NULL,subtitle = NULL)

evalbox_rna <- VhgBoxplot(taubert_rna ,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",reorder_criteria = "phylum_median",legend_position = "right",xlabel = "RNA viruses",axis_title_size = 12,ytext_size = 12,
                          xtext_size = 11,title = NULL,subtitle = NULL)




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


ExportVirusPlot("taubert_evalue_hunter_grid.png",plot = plot_list,
                width = 9,height = 14,units = "in",limitsize = FALSE,ncol = 1,nrow = 3,greedy = TRUE,align="h",
                labels = c("A","B","C"),path = export_path)


testi <- CombineHittables(taubert_smalldna,taubert_largedna,taubert_rna)

tablo <- VhgRunsTable(testi,groupby = "ViralRefSeq_taxonomy",col_everyrow = TRUE,cell_colour = "white")



tablo <- tablo  |>
  tab_row_group(
    label = "Small DNA Viruses",
    rows = 1
  ) |>
  tab_row_group(
    label = "Large DNA Viruses",
    rows = c(4,5)
  ) |>
  tab_row_group(
    label = "RNA Viruses",
    rows = c(2,3,6)
  ) |>
  # Apply grey color to all row group headers
  tab_style(
    style = list(
      cell_fill(color = "lightgrey"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  )%>% 
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

gtsave(tablo,"summary_unique_runstaubert.png")