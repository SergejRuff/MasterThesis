rm(list=ls())

library(Virusparies)
library(ggplot2)
library(tidyverse)
library(gt)

m1 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virusgatherer-cap3.tsv")
m2 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/hepevirga/virusgatherer-cap3.tsv")
m3 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Nido/virusgatherer-cap3.tsv")
m4 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Bunya/virusgatherer-cap3.tsv")
m5 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Mono_chu_08july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")
m6 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/orthomyxo_20july_RNAvirus_nofil_1/virusgatherer-cap3.tsv")
m7 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/NidoAstro_RdRp/virusgatherer-cap3.tsv")
m8 <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/mammals/Nido_NiRAN/virusgatherer-cap3.tsv")


combined_ga <- CombineHittables(m1,m2,m3,m4,m5,m6,m7,m8)

# combined_ga <- ImportVirusTable("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.tsv")

combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")


combined_ga$Category <- with(combined_ga, 
                             ifelse(grepl("viridae$", ViralRefSeq_taxonomy, ignore.case = TRUE), "Classified Families", 
                                    ifelse(grepl("^unclassified\\s+\\w+", ViralRefSeq_taxonomy, ignore.case = TRUE), "Unclassified + Phylum", 
                                           "Unclassified Families")))

# combined_ga <- VhgSubsetHittable(combined_ga,num_hits_min = 4,ViralRefSeq_E_criteria = 1e-5,ViralRefSeq_ident_criteria = -90)

combined_ga <- VhgSubsetHittable(combined_ga,ViralRefSeq_E_criteria = 1e-5)

total_count <- nrow(combined_ga)

# Define the updated color palette
category_colors <- c("Classified Families" = "violet",  # Yellow
                     "Unclassified Families" = "#999999",  
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
class_vs_noclass <- ggplot(combined_ga, aes(x = Category, fill = Category)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = category_colors) +  # Apply custom colors
  labs(
       subtitle = paste("Total count:", total_count),
       x = "Category",
       y = "Count") +
  base_theme +  # Apply the customized theme
  geom_text(stat = 'count', aes(label = paste0(..count.., " (", round(..count.. / total_count * 100, 2), "%)")), vjust = -0.5)


ggsave("class_v-unclassified_22sep.png",plot = class_vs_noclass,path ="output/mammals/",units = "px", 
       width = 600,height = 800 )
## Generate Gatherer Plots

# title = "Mammal sequence data (Gatherer)\nDistribution of viral groups detected across query sequences"
sra <- VhgRunsBarplot(combined_ga,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
                      plot_text = 1,title = NULL,group_unwanted_phyla = "rna",reorder_criteria = "phylum_max",xlabel = "RNA viruses")

# title = "Mammal sequence data (Gatherer)\nBoxplot of viral reference E-values for each group"
boxplot_ <- VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_E",
                       theme_choice = "linedraw_dotted",legend_position = "right",
                       title = NULL,group_unwanted_phyla = "rna",reorder_criteria = "phylum_median",xlabel = "RNA viruses")

boxplot_iden <- VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",
                           theme_choice = "linedraw_dotted",group_unwanted_phyla = "rna",reorder_criteria = "phylum_median",title = NULL,xlabel = "RNA viruses")

boxplot_con <-VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",contiglen_log10_scale = TRUE,
                         theme_choice = "linedraw_dotted",group_unwanted_phyla = "rna",reorder_criteria = "phylum_median",title = NULL,xlabel = "RNA viruses")

eval_stats <- VhgTabularRasa(boxplot_ $summary_stats,title = NULL,names_ = c("Viral reference taxonomy","Median","Q1","Q3",
                                                                             "Mean","SD","Min","Max"),col_everyrow = TRUE,cell_colour = "white")


iden_stats <- VhgTabularRasa(boxplot_iden$summary_stats,title = NULL,names_ = c("Viral reference taxonomy","Median","Q1","Q3",
                                                                                "Mean","SD","Min","Max"),col_everyrow = TRUE,cell_colour = "white")


con_stats <-VhgTabularRasa(boxplot_con$summary_stats,title = NULL,names_ = c("Viral reference taxonomy","Median","Q1","Q3",
                                                                              "Mean","SD","Min","Max"),col_everyrow = TRUE,cell_colour = "white")


ExportVirusGt(eval_stats,filename = "mammal_evalstats.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)
ExportVirusGt(iden_stats,filename = "mammal_idenstats.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)
ExportVirusGt(con_stats,filename = "mammal_constats.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)


# export

ExportVirusPlot(plot = boxplot_$boxp,"combined_boxplot.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)

ExportVirusPlot(plot = sra$plot,"combined_sra.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)

ExportVirusPlot(plot = boxplot_iden$boxp,"combined_boxplot_iden.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)

ExportVirusPlot(plot = boxplot_con$boxp,"combined_boxplot_con.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 10,height = 13)


subject <- VhgGetSubject(combined_ga,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE,group_unwanted_phyla = "rna")

ExportVirusDataFrame(subject,file_name = "subject_mammal.csv",dir_path = "output/mammals")

stats <- SummarizeViralStats(combined_ga,groupby = "ViralRefSeq_taxonomy",
                             metric = "ViralRefSeq_ident",metric_cutoff = 90,filter_cutoff = 1e-5,
                             show_total = TRUE,extra_stats = c("median","Q1","Q3"),group_unwanted_phyla = "rna",sort_by = "total")

ExportVirusDataFrame(stats,file_name = "identity_mammal.csv",dir_path = "output/mammals")

ExportVirusDataFrame(combined_ga,file_name = "combined_ga.tsv",dir_path = "data/RNAvirus_Mammals_newJan2023/mammals")





##########################################################

# virus_list <- sra[["sample_run"]]$cyl[1]
# virus_list <- strsplit(virus_list, ",\\s*")[[1]]
# 
# # Calculate the number of rows and columns for the 5x5 table
# n_rows <- 6
# n_cols <- ceiling(length(virus_list) / n_rows)
# 
# # Create a matrix for the table with NA for missing values
# matrix_data <- matrix(NA, nrow = n_rows, ncol = n_cols)
# matrix_data[1:length(virus_list)] <- virus_list
# 
# # Convert matrix to dataframe
# df <- as.data.frame(matrix_data, stringsAsFactors = FALSE)
# 
# # Create a gt table
# gt_table <- df %>%
#   gt() %>%
#   tab_header(
#     title = "Non-RNA-viruses found in Mammal data"
#   ) %>%
#   cols_label(
#     V1 = "Family 1-6",
#     V2 = "Family 7-12"
#   ) %>%
#   fmt_markdown(
#     columns = everything()
#   )%>% 
#   opt_align_table_header(align = "left") %>% 
#   tab_options(
#     # These were the ones we applied in the first chapter
#     data_row.padding = px(2),
#     summary_row.padding = px(3), # A bit more padding for summaries
#     row_group.padding = px(4)    # And even more for our groups
#   ) %>% opt_stylize(style = 4) %>% 
#   tab_options(
#     heading.title.font.size = px(20)
#   )
# 
# # Print the gt table
# print(gt_table)


virus_list <- sra[["sample_run"]]$cyl[1]
virus_list <- strsplit(virus_list, ",\\s*")[[1]]

# Calculate the number of rows and columns for the 6xN table
n_rows <- 6
n_cols <- ceiling(length(virus_list) / n_rows)

# Create a matrix for the table with NA for missing values
matrix_data <- matrix(NA, nrow = n_rows, ncol = n_cols)
matrix_data[1:length(virus_list)] <- virus_list

# Convert matrix to dataframe
df <- as.data.frame(matrix_data, stringsAsFactors = FALSE)

# Create a gt table
gt_table <- df %>%
  gt() %>%
  tab_header(
    title = NULL
  ) %>%
  cols_label(
    V1 = "Family 1-6",
    V2 = "Family 7-11"
  ) %>%
  fmt_markdown(
    columns = everything()
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = "—"  # Use a dash or another symbol to indicate missing data
  ) %>%
  tab_spanner(
    label = "Virus Families",
    columns = everything()
  ) %>%
  opt_align_table_header(align = "left") %>%
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4)
  ) %>%
  opt_stylize(style = 4) %>%
  tab_options(
    heading.title.font.size = px(20)
  )

# Print the gt table
print(gt_table)


gtsave(gt_table,"non-rna-viruses-in-mammal.png",path = "output/mammals")


########################################

library(dplyr)


# Your tibble data (assuming it's named 'data')
filtered_data <- stats %>%
  filter(!ViralRefSeq_taxonomy %in% c("unclassified", "Non-RNA-viruses")) %>%
  filter(ViralRefSeq_taxonomy != "Total") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10)


ident_gt <- VhgTabularRasa(filtered_data,title = NULL,
                           names_ = c("Viral taxonomy","Identity < 90%", "Identity ≥ 90%","Total","Median","Q1","Q3"),title_colour = "black")


ident_gt <- ident_gt %>%
  fmt_number(
    columns = c("Median", "Q1", "Q3"),
    decimals = 2
  ) %>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )


gtsave(ident_gt,filename = "top10_iden_mammal.png",path = "output/mammals/")


gsa_subset <- VhgSubsetHittable(combined_ga,group_column = "ViralRefSeq_taxonomy",
                                virus_groups = c("Arteriviridae", "Picobirnaviridae","Fiersviridae",
                                                 "Flaviviridae","Astroviridae","Tombusviridae","Steitzviridae",
                                                 "Partitiviridae","Picornaviridae","unclassified Pisuviricota"))

t <- VhgGetSubject(gsa_subset,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE)

# Assuming 't' is your data frame
grouped_data <- t %>%
  group_by(ViralRefSeq_taxonomy, Processed_ViralRefSeq_subject) %>%
  summarise(subject_count = sum(subject_count), .groups = 'drop') %>%
  arrange(desc(subject_count))

# Print the result
print(grouped_data)


gsa_subset_below90 <- VhgSubsetHittable(combined_ga,group_column = "ViralRefSeq_taxonomy",
                                virus_groups = c("Arteriviridae", "Picobirnaviridae","Fiersviridae",
                                                 "Flaviviridae","Astroviridae","Tombusviridae","Steitzviridae",
                                                 "Partitiviridae","Picornaviridae","unclassified Pisuviricota"),ViralRefSeq_ident_criteria = -90)
t90 <- VhgGetSubject(gsa_subset_below90,groupby = "ViralRefSeq_taxonomy",extract_brackets = TRUE)

# Assuming 't' is your data frame
grouped_data90 <- t90 %>%
  group_by(ViralRefSeq_taxonomy, Processed_ViralRefSeq_subject) %>%
  summarise(subject_count = sum(subject_count), .groups = 'drop') %>%
  arrange(desc(subject_count))

grouped_data <- grouped_data %>%
  left_join(
    grouped_data90 %>%
      select(ViralRefSeq_taxonomy, Processed_ViralRefSeq_subject, subject_count) %>%
      rename(subject_90 = subject_count), 
    by = c("ViralRefSeq_taxonomy", "Processed_ViralRefSeq_subject")
  ) %>%
  mutate(subject_90 = ifelse(is.na(subject_90), 0, subject_90)) # Replace NA in subject_90 with 0


grouped_data <- grouped_data%>%
  slice_head(n = 10)

viral_co<- VhgTabularRasa(grouped_data,title = "Top 10 Viral Subjects",
                          names_ = c("Viral taxonomy","Viral subjects","number of contigs","number of contigs below 90% identity"),title_colour = "black")


viral_co <- viral_co%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

gtsave(viral_co,filename = "viral_co_mammal.png",path = "output/mammals/")





################################################


# sra <- VhgRunsBarplot(combined_ga,groupby = "ViralRefSeq_taxonomy",theme_choice = "linedraw_dotted",
#                       plot_text = 1,title = NULL,group_unwanted_phyla = "rna",xlabel = "RNA viruses",legend_position = "none",reorder_criteria = NULL,
#                       xtext_size = 8,plot_text_size = 3.5)
# 
# 
# boxplot_iden <- VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "ViralRefSeq_ident",
#                            theme_choice = "linedraw_dotted",group_unwanted_phyla = "rna",xlabel = "",legend_position = "right",reorder_criteria = NULL,title = NULL,remove_group_labels = TRUE,
#                            xtext_size = 8)
# 
# boxplot_con <-VhgBoxplot(combined_ga,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len",contiglen_log10_scale = TRUE,
#                          theme_choice = "linedraw_dotted",group_unwanted_phyla = "rna",xlabel = "RNA viruses",legend_position = "right",reorder_criteria = NULL,title = NULL,
#                          xtext_size = 8)
# 
# 
# plot_list <- list(sra$plot,boxplot_iden$boxp,boxplot_con$boxp)
# 
# # Adjust y-axis title margin
# plot_list <- lapply(plot_list, function(p) {
#   p + theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
# })
# 
# 
# ExportVirusPlot("taubert_smalldna_grid.png",plot = plot_list,
#                 width = 12,height = 13,units = "in",limitsize = FALSE,ncol = 2,nrow = 2,greedy = TRUE,align="h",
#                 labels = c("A","B","C"))



# Define a function to import the table and add the source column
ImportVirusTableWithSource <- function(file_path, source_name) {
  df <- ImportVirusTable(file_path)
  df$source <- source_name
  return(df)
}

# Import the tables with the source information
m1 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/Flavi/virusgatherer-cap3.tsv", "Flavi")
m2 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/hepevirga/virusgatherer-cap3.tsv", "hepevirga")
m3 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/Nido/virusgatherer-cap3.tsv", "Nido")
m4 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/Bunya/virusgatherer-cap3.tsv", "Bunya")
m5 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/Mono_chu_08july_RNAvirus_nofil_1/virusgatherer-cap3.tsv", "Mono_chu_08july_RNAvirus_nofil_1")
m6 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/orthomyxo_20july_RNAvirus_nofil_1/virusgatherer-cap3.tsv", "orthomyxo_20july_RNAvirus_nofil_1")
m7 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/NidoAstro_RdRp/virusgatherer-cap3.tsv", "NidoAstro_RdRp")
m8 <- ImportVirusTableWithSource("data/RNAvirus_Mammals_newJan2023/mammals/Nido_NiRAN/virusgatherer-cap3.tsv", "Nido_NiRAN")

# Combine the tables
combined_ga <- CombineHittables(m1, m2, m3, m4, m5, m6, m7, m8)

# View the first few rows of the combined data frame to check the results
head(combined_ga)
combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")
ExportVirusDataFrame(combined_ga,file_name = "combined_ga_with_source.tsv",dir_path = "data/RNAvirus_Mammals_newJan2023/mammals")


###################################################

combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")


combined_ga <- VhgAddPhylum(combined_ga,"ViralRefSeq_taxonomy")

combined_ga <- VhgSubsetHittable(combined_ga,ViralRefSeq_E_criteria = 1e-5)

# Assuming combined_ga is your dataframe
result <- combined_ga %>%
  group_by(Phylum, ViralRefSeq_taxonomy) %>%
  summarise(num_observations = n(), .groups = 'drop') %>%
  group_by(Phylum) %>%
  mutate(total_observations = sum(num_observations)) %>%
  ungroup() %>%
  mutate(percentage = (num_observations / total_observations) * 100) %>%
  arrange(Phylum, desc(num_observations))

# View the result
print(result)



#########################

result_iden <- combined_ga %>%
  group_by(Phylum) %>%
  summarise(
    total_contigs = n(),
    contigs_below_90 = sum(ViralRefSeq_ident < 90, na.rm = TRUE),
    contigs_above_90 = sum(ViralRefSeq_ident >= 90, na.rm = TRUE),
    percent_below_90 = round((contigs_below_90 / total_contigs) * 100, 2),
    percent_above_90 = round((contigs_above_90 / total_contigs) * 100, 2)
  )

# Print the result
print(result_iden)

phyl_stats <-VhgTabularRasa(result_iden,title = NULL,names_ = c("Phylum","Number of contigs","Contigs with sequence identity < 90%","Contigs with sequence identity >= 90%",
                                                                              "Percentage < 90%","Percentage >= 90%"))


ExportVirusGt(phyl_stats,filename = "mammal_phyl_stats.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)
########################################



contiglarger1000 <- combined_ga %>% 
  group_by(Phylum) %>% 
  summarise(
    over1000 = sum(contig_len > 1000),
    percentage = round((over1000 / 1225) * 100, 2)
  )

# Assuming VhgSubsetHittable function modifies combined_ga based on criteria
combined_ga_90 <- VhgSubsetHittable(combined_ga, ViralRefSeq_ident_criteria = -90)

# For contiglarger1000_novel
contiglarger1000_novel <- combined_ga_90 %>% 
  group_by(Phylum) %>% 
  summarise(
    over1000 = sum(contig_len > 1000),
    percentage = round((over1000 / 1225) * 100, 2)
  )

contigcombined_data <- cbind(contiglarger1000, contiglarger1000_novel)

contigcombined_data <- contigcombined_data[ , -4]


contigcombined_data_stats <-VhgTabularRasa(contigcombined_data,title = NULL,names_ = c("Phylum","Contig length >= 1000 nt","% of Contigs >= 1000 nt",
                                                                "Contigs length >= 1000 nt & sequence identity < 90%","% of length >= 1000 nt & sequence identity < 90%"))

contigcombined_data_stats <- contigcombined_data_stats%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )

ExportVirusGt(contigcombined_data_stats,filename = "contigcombined_data_stats.png",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)


combined_ga_90_1000 <- VhgSubsetHittable(combined_ga, ViralRefSeq_ident_criteria = -90,contig_len_criteria = 1000)

phylum_smaller90_larger1000 <- combined_ga_90_1000  %>%
  group_by(Phylum, ViralRefSeq_taxonomy) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Phylum) %>%
  mutate(total_phylum = sum(count), 
         percentage = round((count / total_phylum) * 100, 2)) %>%
  ungroup() %>%
  arrange(Phylum, desc(percentage))

phylum_smaller90_larger1000 <- phylum_smaller90_larger1000[ , -4]

phylum_smaller90_larger1000 <-VhgTabularRasa(phylum_smaller90_larger1000,title = NULL,names_ = c("Phylum","Viral reference family","Contigs length >= 1000 nt & sequence identity < 90%",
                                                                                            "% of length >= 1000 nt & sequence identity < 90%"),col_everyrow = TRUE,cell_colour = "white")



phylum_smaller90_larger1000 <- phylum_smaller90_larger1000%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )


ExportVirusGt(phylum_smaller90_larger1000,filename = "phylum_smaller90_larger1000.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)




###################################################################



combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")


combined_ga <- VhgAddPhylum(combined_ga,"ViralRefSeq_taxonomy")

combined_ga <- VhgSubsetHittable(combined_ga,ViralRefSeq_E_criteria = 1e-5)


host_taxon_df <- as.data.frame(table(combined_ga$host_taxon)) %>%
  rename(host_taxon = Var1, count = Freq) %>%
  mutate(total = sum(count),
         percentage = round((count / total) * 100, 2)) %>%
  arrange(desc(percentage))

host_taxon_df <- host_taxon_df[,-3]


host_info_gt <-VhgTabularRasa(host_taxon_df,title = NULL,names_ = c("Host","Number of contigs","Percentage"),col_everyrow = TRUE,cell_colour = "white")

host_info_gt <- host_info_gt%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )



ExportVirusGt(host_info_gt,filename = "host_info_gt.docx",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)


