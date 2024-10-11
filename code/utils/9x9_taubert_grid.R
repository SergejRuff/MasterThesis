rm(list = ls())

library(Virusparies)
library(ggplot2)
library(patchwork)

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


taubert_large_p2 <- VhgBoxplot(taubert_large,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "",
                               axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

taubert_large_p3 <- VhgBoxplot(taubert_large,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                               contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = TRUE,
                               title = NULL,xtext_size = 11,legend_text_size = 12,xlabel ="")


##### RNAvirus

taubert_rna_p1 <- VhgRunsBarplot(taubert_rna,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                 theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",remove_group_labels = FALSE,ytext_size = 9,
                                 title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "RNA viruses",subtitle_size = 11,axis_title_size = 10)


taubert_rna_p2 <- VhgBoxplot(taubert_rna,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",theme_choice = "linedraw_dotted",reorder_criteria = NULL,legend_position = "right",xlabel = "",
                             axis_title_size = 12,ytext_size = 10,xtext_size = 11,title = NULL,remove_group_labels = TRUE)

taubert_rna_p3 <- VhgBoxplot(taubert_rna,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",theme_choice = "linedraw_dotted",reorder_criteria = NULL,
                             contiglen_log10_scale = TRUE,legend_position = "right",remove_group_labels = TRUE,ytext_size = 12,
                             title = NULL,xtext_size = 11,legend_text_size = 12,xlabel ="")



################################# Export

# plot_list <- list(taubert_small_p1$plot,taubert_small_p2$boxp,taubert_small_p3$boxp,
#                   taubert_large_p1$plot,taubert_large_p2$boxp,taubert_large_p3$boxp,
#                   taubert_rna_p1$plot,taubert_rna_p2$boxp,taubert_rna_p3$boxp)

taubert_large_p1$plot <- taubert_large_p1$plot + scale_y_continuous(limits = c(0, 10))
taubert_rna_p1$plot <- taubert_rna_p1$plot + scale_y_continuous(limits = c(0, 10))

plot_list <- list(taubert_small_p1$plot,taubert_small_p2$boxp,taubert_small_p3$boxp,
                  taubert_large_p1$plot,taubert_large_p2$boxp,taubert_large_p3$boxp
                  taubert_rna_p1$plot,taubert_rna_p2$boxp,taubert_rna_p3$boxp)

# Adjust y-axis title margin
plot_list <- lapply(plot_list, function(p) {
  p + theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
})

# 
# ExportVirusPlot("grid_9x9_taubert.png",plot = plot_list,
#                 width = 13,height = 9,units = "in",limitsize = FALSE,ncol = 3,nrow = 3,greedy = TRUE,align="h",
#                 labels = c("A","B","C","D","E","F","G","H","I"),path = export_path,rel_widths = c(1, 0.6, 1))

taubert_small_p1$plot <- taubert_small_p1$plot + scale_y_continuous(limits = c(0, 5))
taubert_large_p1$plot <- taubert_large_p1$plot + scale_y_continuous(limits = c(0, 5))
taubert_rna_p1$plot <- taubert_rna_p1$plot + scale_y_continuous(limits = c(0, 5))

taubert_small_p3$boxp <- taubert_small_p3$boxp + scale_y_continuous(limits = c(50, 270))
taubert_large_p3$boxp <- taubert_large_p3$boxp + scale_y_continuous(limits = c(50, 270))
taubert_rna_p3$boxp <- taubert_rna_p3$boxp + scale_y_continuous(limits = c(50, 270))

combined_gatherer_plots <- (taubert_small_p1$plot|taubert_small_p2$boxp|taubert_small_p3$boxp)/
  (taubert_large_p1$plot|taubert_large_p2$boxp|taubert_large_p3$boxp)/
  (taubert_rna_p1$plot|taubert_rna_p2$boxp|taubert_rna_p3$boxp)+
  plot_annotation(tag_levels = "A") 

ggsave("grid_9x9_taubert_gatherer.png", plot = combined_gatherer_plots, width = 13, height = 9, units = "in")


taubert_large_p15 <- VhSumHitsBarplot(taubert_large,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                      theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",ytext_size = 9,
                                      title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "Large DNA viruses",subtitle_size = 11,axis_title_size = 10,remove_group_labels = TRUE)



taubert_rna_p15 <- VhSumHitsBarplot(taubert_rna,groupby = "ViralRefSeq_taxonomy",reorder_criteria = NULL,
                                    theme_choice = "linedraw_dotted",plot_text = 1,legend_position = "none",ytext_size = 9,
                                    title = NULL,xtext_size = 11,legend_text_size = 12,plot_text_size = 5,xlabel = "RNA viruses",subtitle_size = 11,axis_title_size = 10,remove_group_labels = TRUE)


taubert_large_p15$plot <- taubert_large_p15$plot + scale_y_continuous(limits = c(0, 30))
taubert_rna_p15$plot <- taubert_rna_p15$plot + scale_y_continuous(limits = c(0, 30))

combined_plot <- (taubert_large_p1$plot |taubert_large_p15$plot | taubert_large_p2$boxp) / (taubert_rna_p1$plot |taubert_rna_p15$plot | taubert_rna_p2$boxp)+
  plot_annotation(tag_levels = "A") 

ggsave("grid_3x3_taubert_hunter.png", plot = combined_plot, width = 13, height = 9, units = "in")


pups <- unique(taubert_small$run_id)
pups <- as.data.frame(pups)
pups$pups <- gsub("_R[1,2]$","",pups$pups)
anszahl_an_samples <- gsub("_R[1,2]_001$","",pups$pups)
anszahl_an_samples <- unique(anszahl_an_samples)

testerino <- VhgSubsetHittable(peks,group_column = "ViralRefSeq_taxonomy",ViralRefSeq_E_criteria = 1e-5)

pups <- unique(testerino$run_id)
pups <- as.data.frame(pups)
pups$pups <- gsub("_R[1,2]$","",pups$pups)
anszahl_an_samples <- gsub("_R[1,2]_001$","",pups$pups)
anszahl_an_samples <- unique(anszahl_an_samples)


library(magick)
pipeline <- image_read("/media/sergej/My Book/Masterarbeit/VhVg_workflow_3.drawio(2).png") %>%
  image_ggplot()

combined_gatherer_plots <- (pipeline| (taubert_small_p1$plot|taubert_small_p2$boxp|taubert_small_p3$boxp)/
                (taubert_large_p1$plot|taubert_large_p2$boxp|taubert_large_p3$boxp)/
                (taubert_rna_p1$plot|taubert_rna_p2$boxp|taubert_rna_p3$boxp))+plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(face = 'bold',size = 25))

ggsave("plot_and-pipeline.tiff", plot = combined_gatherer_plots, width = 20, height = 15, units = "in")




########################################################################

combined_ident <- rbind(taubert_small_p2$summary_stats,taubert_large_p2$summary_stats,taubert_rna_p2$summary_stats)


tablo <- VhgTabularRasa(combined_ident,title = NULL,names_ = c("Viral reference taxonomy","Median","Q1","Q3",
                                                      "Mean","SD","Min","Max"))



tablo <- tablo  |>
  tab_row_group(
    label = "Small DNA Viruses",
    rows = 1
  ) |>
  tab_row_group(
    label = "Large DNA Viruses",
    rows = c(2,3)
  ) |>
  tab_row_group(
    label = "RNA Viruses",
    rows = c(4,5,6)
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
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4)
  ) %>%
  opt_stylize(style = 4) %>%
  tab_options(
    heading.title.font.size = px(20)
  )

gtsave(tablo,"taubert_ident_stats.png",path = "output/TaubertDatacombined/")



######################################################################################



combined_ident <- rbind(taubert_small_p3$summary_stats,taubert_large_p3$summary_stats,taubert_rna_p3$summary_stats)


tablo <- VhgTabularRasa(combined_ident,title = NULL,names_ = c("Viral reference taxonomy","Median","Q1","Q3",
                                                               "Mean","SD","Min","Max"))



tablo <- tablo  |>
  tab_row_group(
    label = "Small DNA Viruses",
    rows = 1
  ) |>
  tab_row_group(
    label = "Large DNA Viruses",
    rows = c(2,3)
  ) |>
  tab_row_group(
    label = "RNA Viruses",
    rows = c(4,5,6)
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
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4)
  ) %>%
  opt_stylize(style = 4) %>%
  tab_options(
    heading.title.font.size = px(20)
  )

gtsave(tablo,"taubert_con_stats.png",path = "output/TaubertDatacombined/")