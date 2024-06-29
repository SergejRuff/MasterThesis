# --- clean glob. env. --- #
rm(list=ls())

# --- load packages --- #
library(Virusparies)

vhPipeline <- function(vhg_file,sra_name,virustype,path){
  
  srarun_bar <- VhgRunsBarplot(vh_file = vhg_file,groupby = "ViralRefSeq_taxonomy",
                               title = paste("Distribution of viral groups for",
                                             sra_name,"for",virustype) )
  
  ExportVirusPlot(plot = srarun_bar,file_name = "test.png",path=path)
  
}



path_to_hittable <- 
  "data/hittables_taubert/rnavirus/taubert/160149/hittables/virushunter.tsv"


vhg_file <- 
  ImportVirusTable(path_to_hittable)

sra_name <- "Taubert 150001"
Virustype <- "RNA virus"
path <- "output/test"

vhPipeline(vhg_file,sra_name = sra_name,virustype = Virustype,path = path)
