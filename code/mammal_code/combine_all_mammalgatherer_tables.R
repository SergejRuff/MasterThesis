rm(list=ls())

library(Virusparies)


m1 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virusgatherer-cap3.tsv")
m2 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/hepevirga/virusgatherer-cap3.tsv")
m3 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Nido/virusgatherer-cap3.tsv")
m4 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Bunya/virusgatherer-cap3.tsv")
m5 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Mono_chu_08july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")
m6 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/orthomyxo_20july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")



combined_ga <- CombineHittables(m1,m2,m3,m4,m5,m6)

combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")

## Generate Gatherer Plots

sra <- VhgRunsBarplot(combined_ga,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
               plot_text = 1,title = "Mammal sequence data (Gatherer)\nDistribution of viral groups detected across query sequences")

boxplot_ <- VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_E",
           theme_choice = "linedraw_dotted",legend_position = "right",
           title = "Mammal sequence data (Gatherer)\nBoxplot of viral reference E-values for each group")

VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",
           theme_choice = "linedraw_dotted")

VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",contiglen_log10_scale = TRUE,
           theme_choice = "linedraw_dotted")

# VgConLenViolin(combined_ga)

# VhgIdentityScatterPlot(combined_ga,groupby = "ViralRefSeq_taxonomy",conlen_bubble_plot = TRUE,theme_choice = "linedraw_dotted")


# export

ExportVirusPlot(plot = boxplot_$boxp,"combined_boxplot.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)

ExportVirusPlot(plot = sra$plot,"combined_sra.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)