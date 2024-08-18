# Only to install the current version of Virusparies. Otherwise comment out
# remove.packages("Virusparies") # remove old version before installing new
# library(remotes)
# remotes::install_github("SergejRuff/Virusparies")


# --- clean glob. env. --- #
rm(list=ls())

# --- load packages --- #
library(Virusparies)


# --- Pipeline for Hunter --- #

vhPipeline <- function(vh_file,sra_name,virustype,path,facet_column=NULL,groupby){
  
  message("\n Performing VhgPreprocessTaxa for Hunter.\n")
  vh_file <- VhgPreprocessTaxa(vh_file,"Family")
  
  message("\n VhgRunsBarplot for Hunter.\n")
  
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vh_file,groupby = groupby,
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 16,facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,plot_text_size = 3.5,axis_title_size = 16,
                               subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n VhSumHitsBarplot for Hunter.\n") 
  
  # - sum of hits plot
  sumhitbar <- VhSumHitsBarplot(vh_file,groupby = groupby,title = 
                                  paste0(sra_name," - ",virustype,": Distribution of hits for each virus group"),
                                title_size = 20,facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,plot_text_size = 3.5,axis_title_size = 16,
                                subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  
  
  message("\n VhgIdentityScatterPlot for Hunter.\n")
  
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vh_file,groupby = groupby,title = 
                                            paste0(sra_name," - ",virustype,
                                                   ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 20,legend_title = "Group",ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                                          subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n VhgIdenFacetedScatterPlot for Hunter.\n")
  
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vh_file,groupby = groupby,title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 20,wrap_ncol = 3,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                                                  subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n  VhgBoxplot E-Value for Hunter.\n")
  
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vh_file,x_column  = groupby,title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 20, y_column = "ViralRefSeq_E",facet_ncol = facet_column,ytext_size = 13,xtext_size = 15,axis_title_size = 16,
                           subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n  VhgBoxplot identity for Hunter.\n")
  
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vh_file,x_column  = groupby,title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 20, y_column = "ViralRefSeq_ident",facet_ncol = 
                            facet_column,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                          subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  
  # --- Generate tables --- #
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  message("\n  Generating Tables for Hunter.\n")
  
  if(!is.null(vh_file) ){
    
    
    message("\n  VhgRunsTable for Hunter.\n")
    
    # - table of viral groups detected - in which data set were they detected?
    run_table <- VhgRunsTable(vh_file,groupby  = "best_query",
                              title = paste0(sra_name," - ",virustype,
                                             ": Table of viral groups detected across query sequences"),
                              title_size = 20)
    
    ExportVirusGt(gtable=run_table,filename="vh_runtable.docx",path = path)
    
  }else{
    
    message("\n  Skipping VhgRunsTable for Hunter ( 0 obs.).\n")
    
  }
  
  
  
  if(!is.null(box_evalue$summary_stats)){
    
    
    message("\n  VhgTabularRasa for e-value sum stats for Hunter.\n")
    
    evalues <- VhgTabularRasa(box_evalue$summary_stats,
                              title = 
                                paste0(sra_name," - ",virustype,": E-values for each group"),
                              names_ = c("group","median","Q1","Q3","mean","sd","min","max"),
                              title_size = 20)
    
    ExportVirusGt(gtable=evalues,filename="vh_evaluestable.docx",path = path)
    
  }else{
    
    message("\n  Skipping VhgTabularRasa for e-value sum stats for Hunter ( 0 obs.).\n")
    
  }
  
  
  
  if(!is.null(iden_boxp$summary_stats)){
    
    
    message("\n  VhgTabularRasa for identity sum stats for Hunter.\n")
    
    identity <- VhgTabularRasa(iden_boxp$summary_stats,
                               title = 
                                 paste0(sra_name," - ",virustype,": Identity [%] for each group"),
                               names_ = c("group","median","Q1","Q3","mean","sd","min","max"),
                               title_size = 20)
    
    ExportVirusGt(gtable=identity,filename="vh_identitytable.docx",path = path)
    
  }else{
    
    message("\n  Skipping VhgTabularRasa for identity sum stats for Hunter ( 0 obs.).\n")
    
  }
  
  
  
  
  
  
  
  # --- export --- #
  
  # Export plots
  
  
  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "hunter_srarun.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = sumhitbar$plot,file_name = "hunter_sumhitbar.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = identityplot$plot,file_name = "hunter_identityscatter.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = facetedscatterplot$plot,file_name = "hunter_facetedidentityscatter.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = box_evalue$boxp,file_name = "hunter_boxeval.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = iden_boxp$boxp,file_name = "hunter_identboxp.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  
  
  
  
  
  
  
  return(list(srarun_bar=srarun_bar,sumhitbar=sumhitbar,identityplot=identityplot,
              facetedscatterplot=facetedscatterplot,box_evalue=box_evalue,
              iden_boxp=iden_boxp,run_table=run_table))  
}




######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################


# --- Pipeline for Gatherer --- #

vgPipeline <- function(vg_file,sra_name,virustype,path,facet_column=NULL){
  
  message("\n Performing VhgPreprocessTaxa for Gatherer.\n")
  vg_file <- VhgPreprocessTaxa(vg_file,"Family")
  
  message("\n VhgRunsBarplot for Gatherer.\n")
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vg_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 16,facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,plot_text_size = 3.5,axis_title_size = 16,
                               subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted" )
  
  message("\n VhgIdentityScatterPlot for Gatherer.\n")
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",conlen_bubble_plot = TRUE,title = 
                                            paste0(sra_name," - ",virustype,
                                                   ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 20,legend_title = "Group",legend_position = "right",ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                                          subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n VhgIdenFacetedScatterPlot for Gatherer.\n")
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",conlen_bubble_plot = TRUE,title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 20,wrap_ncol = 3,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                                                  subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n VhgBoxplot E-Value for Gatherer.\n")
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 20, y_column = "ViralRefSeq_E",facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                           subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  message("\n VhgBoxplot Identity for Gatherer.\n")
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 20, y_column = "ViralRefSeq_ident",facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                          subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  violin_contiglen <- VgConLenViolin(vg_file,title = 
                                       paste0(sra_name," - ",virustype,
                                              ": Violinplot of contig length for each group")
                                     ,title_size = 20,facet_ncol = facet_column,ytext_size = 15,xtext_size = 15,axis_title_size = 16,
                                     subtitle_size = 16,legend_text_size = 14,legend_title_size = 16,theme_choice = "linedraw_dotted")
  
  
  # --- Generate tables --- #
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  message("\n  Generating Tables for Hunter.\n")
  
  if(!is.null(vg_file) ){
    
    message("\n  VhgRunsTable for Gatherer.\n")
    
    # - table of viral groups detected - in which data set were they detected?
    run_table <- VhgRunsTable(vg_file,groupby  = "ViralRefSeq_taxonomy",
                              title = paste0(sra_name," - ",virustype,
                                             ": Table of viral groups detected across query sequences"),
                              title_size = 20)
    
    ExportVirusGt(gtable=run_table,filename="vg_runtable.docx",path = path)
    
  }else{
    
    message("\n  Skipping VhgRunsTable for Gatherer ( 0 obs.).\n")
    
  }
  
  if(!is.null(box_evalue$summary_stats)){
    
    message("\n  VhgTabularRasa for e-value sum stats for Hunter.\n")
    
    evalues <- VhgTabularRasa(box_evalue$summary_stats,
                             title = 
                                paste0(sra_name," - ",virustype,": E-values for each group"),
                              names_ = c("group","median","Q1","Q3","mean","sd","min","max"),
                             title_size = 20)
    
    ExportVirusGt(gtable=evalues,filename="vg_evaluestable.docx",path = path)
  }else{
    
    message("\n  Skipping VhgTabularRasa for e-value sum stats for Hunter ( 0 obs.).\n")
    
  }
  
  
  
  if(!is.null(iden_boxp$summary_stats)){
    
    message("\n  VhgTabularRasa for identity sum stats for Hunter.\n")
    
    identity <- VhgTabularRasa(iden_boxp$summary_stats,
                              title = 
                                 paste0(sra_name," - ",virustype,": Identity [%] for each group"),
                               names_ = c("group","median","Q1","Q3","mean","sd","min","max"),
                              title_size = 20)
    
    ExportVirusGt(gtable=identity,filename="vg_identitytable.docx",path = path)
    
    
    
  }else{
    
    message("\n  Skipping VhgTabularRasa for identity sum stats for Hunter ( 0 obs.).\n")
    
  }
  
  
  if(!is.null(violin_contiglen$contiglen_stats)){
    
    message("\n  VhgTabularRasa for contig length sum stats for Hunter.\n")
    
    contig_len <- VhgTabularRasa(violin_contiglen$contiglen_stats,
                                title = 
                                   paste0(sra_name," - ",virustype,": Contig Length [nt] for each group"),
                                 names_ = c("group","median","Q1","Q3","mean","sd","min","max"),
                                title_size = 20)
    
    ExportVirusGt(gtable=contig_len,filename="vg_contig_lentable.docx",path = path)
    
  }else{
    
    message("\n  Skipping VhgTabularRasa for contig length sum stats for Hunter ( 0 obs.).\n")
    
  }
  
  
  
  
  # --- export --- #
  
  # Export plots
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  
  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "gatherer_srarun.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = violin_contiglen$plot,file_name = "gatherer_violinplot.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = identityplot$plot,file_name = "gatherer_identityscattereplot.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = facetedscatterplot$plot,file_name = "gatherer_facetedidentityplot.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = box_evalue$boxp,file_name = "gatherer_box_eval.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  ExportVirusPlot(plot = iden_boxp$boxp,file_name = "gatherer_box_identity.png",path=path,export_plotobj = TRUE,
                  width = 7,height = 6,units = "in",scale = 2)
  
  
  
  
  
  
  
  return(list(srarun_bar=srarun_bar,violin_contiglen=violin_contiglen,identityplot=identityplot,
              facetedscatterplot=facetedscatterplot,box_evalue=box_evalue,
              iden_boxp=iden_boxp,run_table=run_table))  
  
}

# --- Import and parameters --- #

path_to_hunter <- 
  "data/hittables_taubert/smalldna/combined_virushunter.tsv"
path_to_gatherer <- "data/hittables_taubert/smalldna/combined_virusgatherer-cap3.tsv"



sra_name <- "Taubert Data"
virustype <- "Small DNA (Hunter)"
virustype_gatherer <- "Small DNA (Gatherer)"
hunter_export_path <- "output/mammals/TaubertDatacombined/plots/smalldna/Hunter"
gatherer_export_path <- "output/mammals/TaubertDatacombined/plots/smalldna/Gatherer"

facet_column <- NULL
group_ <- "ViralRefseq_taxonomy"


# sra_name <- "Florian´s Data"
# virustype <- "Large DNA viruses (Hunter)"
# virustype_gatherer <- "Large DNA viruses (Gatherer)"
# hunter_export_path <- "output/Florian/largedna/plots/Hunter"
# gatherer_export_path <- "output/Florian/largedna/plots/Gatherer"



if(path_to_hunter != ""){
  
  vh_file <- 
    ImportVirusTable(path_to_hunter)
  
  vh_results <- vhPipeline(vh_file,sra_name = sra_name,virustype = virustype,path = hunter_export_path,facet_column=facet_column,groupby=group_)
  
}


if(path_to_gatherer != ""){
  
  vg_file <- 
    ImportVirusTable(path_to_gatherer)
  
  vg_results <- vgPipeline(vg_file,sra_name = sra_name,virustype = virustype_gatherer,path = gatherer_export_path,facet_column=facet_column)
  
}