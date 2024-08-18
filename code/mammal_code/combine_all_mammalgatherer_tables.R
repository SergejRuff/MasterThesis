rm(list=ls())

library(Virusparies)
library(ggplot2)


m1 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virusgatherer-cap3.tsv")
m2 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/hepevirga/virusgatherer-cap3.tsv")
m3 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Nido/virusgatherer-cap3.tsv")
m4 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Bunya/virusgatherer-cap3.tsv")
m5 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Mono_chu_08july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")
m6 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/orthomyxo_20july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")



combined_ga <- CombineHittables(m1,m2,m3,m4,m5,m6)

combined_ga <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv")

combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")


combined_ga$Category <- with(combined_ga, 
                             ifelse(grepl("viridae$", ViralRefSeq_taxonomy, ignore.case = TRUE), "Classified", 
                                    ifelse(grepl("^unclassified\\s+\\w+", ViralRefSeq_taxonomy, ignore.case = TRUE), "Unclassified + Phylum", 
                                           "Unclassified")))

combined_ga <- VhgSubsetHittable(combined_ga,num_hits_min = 4,ViralRefSeq_E_criteria = 1e-5,ViralRefSeq_ident_criteria = -90)

combined_ga <- VhgSubsetHittable(combined_ga,ViralRefSeq_E_criteria = 1e-5)

total_count <- nrow(combined_ga)

# Define the updated color palette
category_colors <- c("Classified" = "violet",  # Yellow
                     "Unclassified" = "#999999",  
                     "Unclassified + Phylum" = "#E69F00")  # Gray


# Set up the base theme
base_theme <- theme_linedraw()

# Add dotted grid lines
base_theme <- base_theme +
  theme(
    panel.grid.major = element_line(linewidth = 0.5, linetype = "dotted", colour = "grey50"),  # Dotted major grid lines
    panel.grid.minor = element_line(linewidth = 0.25, linetype = "dotted", colour = "grey75")   # Dotted minor grid lines
  )

# Create the bar plot with the customized theme
ggplot(combined_ga, aes(x = Category, fill = Category)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = category_colors) +  # Apply custom colors
  labs(title = "Mammal Data - Distribution of Classified vs. Unclassified Viral Observations",
       subtitle = paste("Total count:", total_count),
       x = "Category",
       y = "Count") +
  base_theme +  # Apply the customized theme
  geom_text(stat = 'count', aes(label = paste0(..count.., " (", round(..count.. / total_count * 100, 1), "%)")), vjust = -0.5)


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