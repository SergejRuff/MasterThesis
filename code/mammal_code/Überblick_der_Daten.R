library(Virusparies)

file <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv")

file_filter <- VhgSubsetHittable(file,num_hits_min = 4,ViralRefSeq_E_criteria = 1e-5,
                                ViralRefSeq_ident_criteria = -90)




sra_runfil <- VhgRunsBarplot(file_filter,plot_text = 1,theme_choice = "linedraw_dotted",title = "Mammal data: Distribution of viral groups",
               title_size = 14,axis_title_size = 12,xtext_size = 12)


path <- "output/mammals/"


ExportVirusPlot(plot = sra_runfil$plot,file_name = "gesamtuberblick_mammals.png",path=path,
                width = 7,height = 7,units = "in")