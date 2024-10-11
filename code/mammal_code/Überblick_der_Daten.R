library(Virusparies)

file <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv")

file_filter <- VhgSubsetHittable(file,num_hits_min = 5,ViralRefSeq_E_criteria = 1e-5,
                                ViralRefSeq_ident_criteria = -90)




sra_runfil <- VhgRunsBarplot(file_filter,plot_text = 1,theme_choice = "linedraw_dotted",title = NULL,
               title_size = 11,axis_title_size = 12,xtext_size = 12,subtitle_size = 11)


path <- "output/mammals/"


ExportVirusPlot(plot = sra_runfil$plot,file_name = "gesamtuberblick_mammals.png",path=path,
                width = 8,height = 7,units = "in")


