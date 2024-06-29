# Only to install the current version of Virusparies. Otherwise comment out
# remove.packages("Virusparies") # remove old version before installing new
# library(remotes)
# remotes::install_github("SergejRuff/Virusparies")


# --- clean glob. env. --- #
rm(list=ls())

# --- load packages --- #
library(Virusparies)


# --- Pipeline for Hunter --- #

vhPipeline <- function(vhg_file,sra_name,virustype,path){
  
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vhg_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 12 )
  
  # - sum of hits plot
  sumhitbar <- VhSumHitsBarplot(vhg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                  paste0(sra_name," - ",virustype,": Distribution of hits for each virus group"),title_size = 12)
  
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vhg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                            paste0(sra_name," - ",virustype,
                                            ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 12,legend_title = "Group")
  
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vhg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 12)
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vhg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 12, y_column = "ViralRefSeq_E",subtitle_size = 10)
  
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vhg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 12, y_column = "ViralRefSeq_ident")
  
  
  # --- Generate tables --- #
  
  # - table of viral groups detected - in which data set were they detected?
  VhgRunsTable(vhg_file,groupby  = "ViralRefSeq_taxonomy",
               title = paste0(sra_name," - ",virustype,
                              ": Table of viral groups detected across query sequences"),
               title_size = 16)
  
  
  # --- export --- #
  
  # Export plots
  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "test.png",path=path,export_plotobj = TRUE)
  
}


# --- Pipeline for Gatherer --- #

vgPipeline <- function(vg_file,sra_name,virustype,path){
  
  # Run Bar Chart - Number of viral groups detected across query sequences
  srarun_bar <- VhgRunsBarplot(file = vg_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste0(sra_name," - ",virustype,
                                              ": Distribution of viral groups detected across query sequences")
                               ,title_size = 12 )
  
  # - sum of hits plot
  sumhitbar <- VhSumHitsBarplot(vg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                  paste0(sra_name," - ",virustype,": Distribution of hits for each virus group"),title_size = 12)
  
  # - scatter plot e value vs identity
  identityplot <-  VhgIdentityScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                            paste0(sra_name," - ",virustype,
                                                   ": Scatterplot of viral reference e-values and identity")
                                          ,title_size = 12,legend_title = "Group")
  
  # - faceted scatter plot identity vs e values faceted by viral group
  facetedscatterplot <- VhgIdenFacetedScatterPlot(file = vg_file,groupby = "ViralRefSeq_taxonomy",title = 
                                                    paste0(sra_name," - ",virustype,
                                                           ": Faceted scatterplot of viral reference e-values and identity")
                                                  ,title_size = 12)
  # - boxplot for e values
  box_evalue <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                             paste0(sra_name," - ",virustype,
                                    ": Boxplot of viral reference e-values for each group")
                           ,title_size = 12, y_column = "ViralRefSeq_E",subtitle_size = 10)
  
  # - boxplot for identity
  iden_boxp <- VhgBoxplot(file = vg_file,x_column  = "ViralRefSeq_taxonomy",title = 
                            paste0(sra_name," - ",virustype,
                                   ": Boxplot of viral reference identity for each group")
                          ,title_size = 12, y_column = "ViralRefSeq_ident")
  
  
  # --- Generate tables --- #
  
  # - table of viral groups detected - in which data set were they detected?
  VhgRunsTable(vg_file,groupby  = "ViralRefSeq_taxonomy",
               title = paste0(sra_name," - ",virustype,
                              ": Table of viral groups detected across query sequences"),
               title_size = 16)
  
  
  # --- export --- #
  
  # Export plots
  
  ExportVirusPlot(plot = srarun_bar$plot,file_name = "test.png",path=path,export_plotobj = TRUE)
  
}

# --- Import and parameters --- #

path_to_hittable <- 
  "data/hittables_taubert/rnavirus/taubert/160149/hittables/virushunter.tsv"


vhg_file <- 
  ImportVirusTable(path_to_hittable)

sra_name <- "Taubert 150001"
virustype <- "RNA viruses"
path <- "output/test"

vhPipeline(vhg_file,sra_name = sra_name,virustype = Virustype,path = path)
