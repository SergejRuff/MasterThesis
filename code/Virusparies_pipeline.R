# Only to install the current version of Virusparies. Otherwise comment out
# remove.packages("Virusparies") # remove old version before installing new
# library(remotes)
# remotes::install_github("SergejRuff/Virusparies")


# --- clean glob. env. --- #
rm(list=ls())

# --- load packages --- #
library(Virusparies)


# --- Pipeline for Hunter --- #

vhPipeline <- function(vh_file,sra_name,virustype,path){
  
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vh_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 12 )
  
  # - sum of hits plot
  sumhitbar <- VhSumHitsBarplot(vh_file,groupby = "ViralRefSeq_taxonomy",title = 
                                  paste0(sra_name," - ",virustype,": Distribution of hits for each virus group"),title_size = 12)
  
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vh_file,groupby = "ViralRefSeq_taxonomy",title = 
                                            paste0(sra_name," - ",virustype,
                                            ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 16,legend_title = "Group")
  
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vh_file,groupby = "ViralRefSeq_taxonomy",title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 16,wrap_ncol = 3)
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vh_file,x_column  = "ViralRefSeq_taxonomy",title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 16, y_column = "ViralRefSeq_E",subtitle_size = 10)
  
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vh_file,x_column  = "ViralRefSeq_taxonomy",title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 16, y_column = "ViralRefSeq_ident")
  
  
  # --- Generate tables --- #
  
  # - table of viral groups detected - in which data set were they detected?
  run_table <- VhgRunsTable(vh_file,groupby  = "ViralRefSeq_taxonomy",
               title = paste0(sra_name," - ",virustype,
                              ": Table of viral groups detected across query sequences"),
              col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
              title_weight = "normal",title_size = 20)
  
  evalues <- VhgTabularRasa(box_evalue$summary_stats,
                            col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
                            title_weight = "normal",title_size = 20,title = 
                              paste0(sra_name," - ",virustype,": E-values for each group"),
                            names_ = c("group","median","Q1","Q3","mean","sd","max","min"))
  
  
  
  identity <- VhgTabularRasa(iden_boxp$summary_stats,
                            col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
                            title_weight = "normal",title_size = 20,title = 
                              paste0(sra_name," - ",virustype,": Identity [%] for each group"),
                            names_ = c("group","median","Q1","Q3","mean","sd","max","min"))
  
  
  
  
  # --- export --- #
  
  # Export plots
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "hunter_srarun.png",path=path,export_plotobj = TRUE,
                  width = 13,height = 10,units = "in")
  
  ExportVirusPlot(plot = sumhitbar$plot,file_name = "hunter_sumhitbar.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = identityplot$plot,file_name = "hunter_identityscatter.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = facetedscatterplot$plot,file_name = "hunter_facetedidentityscatter.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = box_evalue$boxp,file_name = "hunter_boxeval.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = iden_boxp$boxp,file_name = "hunter_identboxp.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  # export tables
  
  ExportVirusGt(gtable=run_table,filename="vh_runtable.docx",path = path)
  
  ExportVirusGt(gtable=evalues,filename="vh_evaluestable.docx",path = path)
  
  ExportVirusGt(gtable=identity,filename="vh_identitytable.docx",path = path)
  
  
 
  return(list(srarun_bar=srarun_bar,sumhitbar=sumhitbar,identityplot=identityplot,
              facetedscatterplot=facetedscatterplot,box_evalue=box_evalue,
              iden_boxp=iden_boxp,run_table=run_table))  
}


# --- Pipeline for Gatherer --- #

vgPipeline <- function(vg_file,sra_name,virustype,path){
  
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vg_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 12 )
  
  
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",conlen_bubble_plot = TRUE,title = 
                                            paste0(sra_name," - ",virustype,
                                                   ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 16,legend_title = "Group",legend_position = "right")
  
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",conlen_bubble_plot = TRUE,title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 16,wrap_ncol = 3)
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 16, y_column = "ViralRefSeq_E",subtitle_size = 10)
  
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 16, y_column = "ViralRefSeq_ident")
  
  violin_contiglen <- VgConLenViolin(vg_file,title = 
                                       paste0(sra_name," - ",virustype,
                                              ": Violinplot of contig length for each group")
                                     ,title_size = 16)
  
  
  # --- Generate tables --- #
  
  # - table of viral groups detected - in which data set were they detected?
  run_table <- VhgRunsTable(vg_file,groupby  = "ViralRefSeq_taxonomy",
               title = paste0(sra_name," - ",virustype,
                              ": Table of viral groups detected across query sequences"),
               ,
               col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
               title_weight = "normal",title_size = 20)
  
  
  evalues <- VhgTabularRasa(box_evalue$summary_stats,
                            col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
                            title_weight = "normal",title_size = 20,title = 
                              paste0(sra_name," - ",virustype,": E-values for each group"),
                            names_ = c("group","median","Q1","Q3","mean","sd","max","min"))
  
  
  
  identity <- VhgTabularRasa(iden_boxp$summary_stats,
                             col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
                             title_weight = "normal",title_size = 20,title = 
                               paste0(sra_name," - ",virustype,": Identity [%] for each group"),
                             names_ = c("group","median","Q1","Q3","mean","sd","max","min"))
  
  
  contig_len <- VhgTabularRasa(violin_contiglen$contiglen_stats,
                             col_everyrow = TRUE,column_colour = "white",cell_colour = "white",title_colour ="black",
                             title_weight = "normal",title_size = 20,title = 
                               paste0(sra_name," - ",virustype,": Contig Length [nt] for each group"),
                             names_ = c("group","median","Q1","Q3","mean","sd","max","min"))
  
  
  
  # --- export --- #
  
  # Export plots
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  

  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "gatherer_srarun.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = violin_contiglen$plot,file_name = "gatherer_violinplot.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = identityplot$plot,file_name = "gatherer_identityscattereplot.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = facetedscatterplot$plot,file_name = "gatherer_facetedidentityplot.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = box_evalue$boxp,file_name = "gatherer_box_eval.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  ExportVirusPlot(plot = iden_boxp$boxp,file_name = "gatherer_box_identity.png",path=path,export_plotobj = TRUE,
                  width = 12,height = 10,units = "in")
  
  # export tables
  
  ExportVirusGt(gtable=run_table,filename="vg_runtable.docx",path = path)
  
  ExportVirusGt(gtable=evalues,filename="vg_evaluestable.docx",path = path)
  
  ExportVirusGt(gtable=identity,filename="vg_identitytable.docx",path = path)
  
  ExportVirusGt(gtable=contig_len,filename="vg_contig_lentable.docx",path = path)
  
  

  
  
  return(list(srarun_bar=srarun_bar,violin_contiglen=violin_contiglen,identityplot=identityplot,
              facetedscatterplot=facetedscatterplot,box_evalue=box_evalue,
              iden_boxp=iden_boxp,run_table=run_table))  
  
}

# --- Import and parameters --- #

path_to_hunter <- 
  "data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virushunter.tsv"
path_to_gatherer <- ""



sra_name <- "Flavi Data"
virustype <- "RNA viruses (Hunter)"
virustype_gatherer <- "RNA viruses (Gatherer)"
hunter_export_path <- "output/Flavi/plots/Hunter"
gatherer_export_path <- "output/Flavi/plots/Gatherer"


# sra_name <- "Florian´s Data"
# virustype <- "Large DNA viruses (Hunter)"
# virustype_gatherer <- "Large DNA viruses (Gatherer)"
# hunter_export_path <- "output/Florian/largedna/plots/Hunter"
# gatherer_export_path <- "output/Florian/largedna/plots/Gatherer"



if(path_to_hunter != ""){
  
  vh_file <- 
    ImportVirusTable(path_to_hunter)
  
  vh_results <- vhPipeline(vh_file,sra_name = sra_name,virustype = virustype,path = hunter_export_path)
  
}


if(path_to_gatherer != ""){
  
  vg_file <- 
    ImportVirusTable(path_to_gatherer)
  
  vg_results <- vgPipeline(vg_file,sra_name = sra_name,virustype = virustype_gatherer,path = gatherer_export_path)
  
}
