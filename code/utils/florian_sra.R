rm(list=ls())

library(Virusparies)
library(ggplot2)
library(tidyverse)
library(gt)


m1 <- ImportVirusTable("data/Florian_data/RNA/virusgatherer-cap3.tsv")
m2 <- ImportVirusTable("data/Florian_data/largedna/virusgatherer-cap3.tsv")
m3 <- ImportVirusTable("data/Florian_data/smalldna/virusgatherer-cap3.tsv")

m1 <- VhgPreprocessTaxa(m1,taxa_rank = "Family")
m2 <- VhgPreprocessTaxa(m2,taxa_rank = "Family")
m3 <- VhgPreprocessTaxa(m3,taxa_rank = "Family")


sra_m1 <- VhgRunsBarplot(m1,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
                      plot_text = 1,title = "Pfaff data (Gatherer)\nDistribution of viral groups detected in queries",xlabel = "RNA viruses")


sra_m2 <- VhgRunsBarplot(m2,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
                         plot_text = 1,title = "Pfaff data (Gatherer)\nDistribution of viral groups detected in queries",xlabel = "Large DNA viruses")


sra_m3 <- VhgRunsBarplot(m3,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
                         plot_text = 1,title = "Pfaff data (Gatherer)\nDistribution of viral groups detected in queries",xlabel = "Small DNA viruses")


ExportVirusPlot(plot = sra_m1$plot,"sra_m1.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m2$plot,"sra_m2.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m3$plot,"sra_m3.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)




##############################################


sra_m1 <- VhgBoxplot(m1,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted"
                     ,title = "Pfaff data (Gatherer)\nViral reference identity",xlabel = "RNA viruses",y_column = "ViralRefSeq_ident")


sra_m2 <- VhgBoxplot(m2,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",y_column = "ViralRefSeq_ident",title = "Pfaff data (Gatherer)\nViral reference identity",xlabel = "Large DNA viruses")


sra_m3 <- VhgBoxplot(m3,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",title = "Pfaff data (Gatherer)\nViral reference identity",xlabel = "Small DNA viruses",y_column = "ViralRefSeq_ident")

ExportVirusPlot(plot = sra_m1$boxp,"sra_m1_ident.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m2$boxp,"sra_m2_ident.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m3$boxp,"sra_m3_ident.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)


##############################################################

sra_m1 <- VhgBoxplot(m1,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted"
                     ,title = "Pfaff data (Gatherer)\nContig length",xlabel = "RNA viruses",y_column = "contig_len")


sra_m2 <- VhgBoxplot(m2,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",y_column = "contig_len",title = "Pfaff data (Gatherer)\nContig length",xlabel = "Large DNA viruses")


sra_m3 <- VhgBoxplot(m3,x_column = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",title = "Pfaff data (Gatherer)\nContig length",xlabel = "Small DNA viruses",y_column = "contig_len")

ExportVirusPlot(plot = sra_m1$boxp,"sra_m1_contig_len.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m2$boxp,"sra_m2_contig_len.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)

ExportVirusPlot(plot = sra_m3$boxp,"sra_m3_contig_len.png",path = "output/Florian/",
                units = "in",limitsize = FALSE,width = 8,height = 6)