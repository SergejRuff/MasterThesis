# --- clean glob. env. --- #
rm(list=ls())

# --- load packages --- #
library(Virusparies)



# --- Import and parameters --- #

path_to_hunter <- 
  "data/Florian_data/RNA/virusgatherer-cap3.tsv"
path_to_gatherer <- "data/Florian_data/RNA/virusgatherer-cap3.tsv"



sra_name <- "Taubert Data"
virustype <- "Small DNA (Hunter)"
virustype_gatherer <- "Small DNA (Gatherer)"
hunter_export_path <- "output/mammals/TaubertDatacombined/plots/smalldna/Hunter"
gatherer_export_path <- "output/mammals/TaubertDatacombined/plots/smalldna/Gatherer"

facet_column <- NULL



sra_name = sra_name
virustype = virustype
path = hunter_export_path
facet_column=facet_column


vg_file <- 
  ImportVirusTable(path_to_gatherer)


message("\n Performing VhgPreprocessTaxa for Hunter.\n")
vg_file <- VhgPreprocessTaxa(vg_file,"Family")

message("\n VhgRunsBarplot for Hunter.\n")

# Run Bar Chart - Number of viral groups detected across query sequences
srarun_bar <- VhgRunsBarplot(file = vg_file,groupby = "ViralRefSeq_taxonomy",
                             title = "A"
                             ,title_size = 20,facet_ncol = facet_column,ytext_size = 10,xtext_size = 10,plot_text_size = 3.5,axis_title_size = 12,
                             subtitle_size = 1,legend_text_size = 10,legend_title_size = 12,theme_choice = "linedraw_dotted",plot_text = 1,subtitle = " ",
                             reorder_criteria = NULL,legend_position = "none",xlabel = "RNA viruses")

message("\n VhSumHitsBarplot for Hunter.\n") 

# - boxplot for e values
box_evalue <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                           "B"
                         ,title_size = 20, y_column = "ViralRefSeq_E",facet_ncol = facet_column,ytext_size = 10,xtext_size = 10,axis_title_size = 12,
                         subtitle_size = 16,,legend_text_size = 10,legend_title_size = 12,,theme_choice = "linedraw_dotted",subtitle = " ",reorder_criteria = NULL,
                         legend_position = "right",xlabel = "")

message("\n VhgBoxplot Identity for Gatherer.\n")
# - boxplot for identity
iden_boxp <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                          "C"
                        ,title_size = 20, y_column = "ViralRefSeq_ident",facet_ncol = facet_column,ytext_size = 10,xtext_size = 10,axis_title_size = 12,
                        subtitle_size = 1,theme_choice = "linedraw_dotted",subtitle = " ",reorder_criteria = NULL,legend_position = "none",xlabel = "RNA viruses")

violin_contiglen <- VgConLenViolin(vg_file,title = 
                                    "D"
                                   ,title_size = 20,facet_ncol = facet_column,ytext_size = 10,xtext_size = 10,axis_title_size = 12,
                                   subtitle_size = 1,legend_text_size = 10,legend_title_size = 12,theme_choice = "linedraw_dotted",
                                   legend_position = "right",subtitle = " ",reorder_criteria = NULL,remove_group_labels = TRUE,xlabel = NULL)






plots <- list(srarun_bar$plot,box_evalue$boxp,iden_boxp$boxp,violin_contiglen$plot)

ExportVirusPlot("test.png",plot = plots,
                width = 12,height = 8,units = "in",limitsize = FALSE,ncol = 2,nrow = 2,greedy = FALSE,align="h")

