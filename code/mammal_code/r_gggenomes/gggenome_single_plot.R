rm(list=ls())


library(gggenomes)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(ggpattern) 
library(Polychrome)
library(colorspace)

# Set working directory
setwd("/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/code/mammal_code/r_gggenomes/")

GenoInfo <- read.delim("interproscan_data 1.tsv")

top_csv <- read.delim("/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/output/mammals/longest_con_gvtable.tsv")
top_csv <- top_csv %>%
  rename(seq_id = contig_id)

# Select all rows, no genus filtering
selected_rows <- GenoInfo

colnames(selected_rows)[colnames(selected_rows) == "seq_length"] <- "sequence_length"
colnames(selected_rows)[colnames(selected_rows) == "sequence_id"] <- "seq_id"

selected_rows <- selected_rows %>%
  filter(!(seq_id %in% c("SRR17345012_cap3_Contig-13", "SRR17345012_cap3_Contig-11","SRR17345013_cap3_Contig-21",
                         "SRR13364363_cap3_Contig-11")))



# Adjust start and end positions
selected_rows  <- selected_rows  %>%
  mutate(
    start = orf_start + (start *3)-1,
    end = orf_start + (end *3)-1
  )
# 
# # Adjust start and end positions
# selected_rows  <- selected_rows  %>%
#   mutate(
#     start = orf_start + start - 1,
#     end = orf_start + end - 1
#   )

selected_rows <- selected_rows %>%
  filter(interpro_description != "-")

# Convert e.value to numeric safely, without introducing NA
selected_rows <- selected_rows %>%
  mutate(e.value_numeric = suppressWarnings(as.numeric(e.value))) %>% # Convert to numeric safely
  filter(!is.na(e.value_numeric) & e.value_numeric <= 1e-5) %>% # Filter out NA and values above threshold
  select(-e.value_numeric)  # Remove the temporary numeric column if not needed



# source("group_all_together.R")

# selected_rows <- selected_rows %>%
#   mutate(end = start + (end - start) * 3)
# 
# # Function to adjust overlaps within each ORF
# adjust_overlaps_within_orfs <- function(data) {
#   # Sort by ORF and start position
#   data <- data %>%
#     arrange(orf_start, start)
#   
#   # Loop through each ORF group
#   data <- data %>%
#     group_by(orf_start, orf_end) %>%
#     mutate(
#       # Adjust start positions within each ORF to avoid overlap
#       start = ifelse(
#         row_number() > 1 & start <= lag(end), 
#         lag(end) + 1, 
#         start
#       )
#     ) %>%
#     ungroup()
#   
#   return(data)
# }
# 
# # Apply the function to adjust overlaps within ORFs
# selected_rows <- adjust_overlaps_within_orfs(selected_rows)

selected_rows <- selected_rows %>%
  filter(start <= end)

source("gggenome.R")

selected_rows <- selected_rows %>%
  group_by(seq_id, interpro_description) %>%
  filter(e.value == min(e.value)) %>%
  ungroup()

selected_rows <- selected_rows %>%
  filter(interpro_description !="PL")


selected_rows$alternate <- ifelse(selected_rows$accession %in% c("orf3211", "orf3253", "orf3272", "orf2627", "orf340"), 1, 0)

# Calculate sequence length based on ORF start and end positions
s0 <- as.data.frame(selected_rows[, "seq_id"], stringsAsFactors = FALSE)
names(s0) <- "seq_id" 
s0 <- s0 %>%
  left_join(top_csv %>% select(seq_id, contig_len), by = "seq_id")
s0 <- s0[!duplicated(s0$seq_id), ]
s0 <- s0 %>%
  rename(length  = contig_len)



# ORF information, ensure no duplicates
orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end","alternate")]
orf0 <- orf0[!duplicated(orf0), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)




# Gene and feature information
g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "interpro_description","start","end","alternate")]


top_csv <- read.delim("/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/output/mammals/longest_con_gvtable.tsv")

top_csv <- top_csv %>%
  rename(seq_id = contig_id)

# Change "unclassified Pisuviricota" to "Unclass. Pisuviricota"
top_csv$ViralRefSeq_taxonomy <- gsub("unclassified Pisuviricota", "Unclass. Pisuviricota", top_csv$ViralRefSeq_taxonomy, ignore.case = TRUE)
top_csv$ViralRefSeq_ident <- paste0("Identity: ", round(top_csv$ViralRefSeq_ident, 2), "%")
# top_csv$host_taxon <- paste0("host: ",top_csv$host_taxon)
#top_csv$contig_len <- paste0("length: ", top_csv$contig_len, " nt")

g0 <- g0 %>%
  left_join(top_csv %>% select(seq_id, Phylum,ViralRefSeq_taxonomy,ViralRefSeq_ident,host_taxon), by = "seq_id")

s0 <- s0 %>%
  left_join(top_csv %>% select(seq_id, Phylum,ViralRefSeq_taxonomy,ViralRefSeq_ident,host_taxon), by = "seq_id")

orf0 <- orf0 %>%
  left_join(top_csv %>% select(seq_id, Phylum,ViralRefSeq_taxonomy,ViralRefSeq_ident,host_taxon), by = "seq_id")

g0 <- g0 %>%
  mutate(seq_id = gsub("_cap3_", "_", seq_id)) %>%
  mutate(seq_id = gsub("Contig", "Con", seq_id))

s0 <- s0 %>%
  mutate(seq_id = gsub("_cap3_", "_", seq_id)) %>%
  mutate(seq_id = gsub("Contig", "Con", seq_id))

orf0 <- orf0  %>%
  mutate(seq_id = gsub("_cap3_", "_", seq_id)) %>%
  mutate(seq_id = gsub("Contig", "Con", seq_id))




##################################
# Plot #################

# custom_colors <- createPalette(length(unique(g0$interpro_description)), c("#ff0000", "#00ff00", "#0000ff"), target = "normal") # Polychrome
# names(custom_colors) <- levels(unique(g0$interpro_description))


s0 <- s0[order(s0$length,decreasing = TRUE), ]
combined_vector_names<- paste(s0$seq_id, s0$ViralRefSeq_taxonomy,s0$ViralRefSeq_ident,s0$host_taxon,s0$contig_len, sep = "\n")
combined_vector_names_orf<- paste(orf0$seq_id, orf0$ViralRefSeq_taxonomy,orf0$ViralRefSeq_ident,orf0$host_taxon,orf0$contig_len ,sep = "\n")
combined_vector_names_go <- paste(g0$seq_id, g0$ViralRefSeq_taxonomy,g0$ViralRefSeq_ident,g0$host_taxon,g0$contig_len ,sep = "\n")

s0$seq_id <- combined_vector_names
orf0$seq_id <- combined_vector_names_orf
g0$seq_id <- combined_vector_names_go

g0 <- g0 %>%
  filter(interpro_description != "-")



plot <- gggenomes(seqs = s0, genes = g0, feats = orf0)+
  geom_seq() +
  geom_feat(aes(y=y+0.12*alternate, yend=y+0.12*alternate), position="identity", size=5,alpha=0.9) +
  geom_gene(aes(y=y+0.12*alternate,fill = interpro_description))+
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     breaks = c(0, 5000, 10000, 15000, 20000, 25000,30000) )+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 14),  # Adjust the size of the legend title
    legend.text = element_text(size = 10)
  ) +
  labs(fill = "Protein domains")+
  theme(
    panel.grid.major = element_line(linewidth = 0.5, linetype = "dotted", colour = "grey50"),  # Dotted major grid lines
    panel.grid.minor = element_line(linewidth = 0.25, linetype = "dotted", colour = "grey75")   # Dotted minor grid lines
  )+
  geom_gene_tag(aes(label=interpro_description), check_overlap = TRUE,angle = 45,vjust = -2,size = 3,hjust = -0.3)+
  theme(legend.position = "none")+geom_bin_label(size = 3, expand_left = .3)


#+
  # scale_fill_manual(values = custom_colors)

virus_name <- "complete_2509"

# Save the plot
ggsave(filename = paste0("genome_map_", gsub("/", "_", virus_name), ".png"), width = 9, height = 11, plot = plot)

rds_filename <- paste0("genome_map_", gsub("/", "_", virus_name), ".rds")

# Save the plot object as an RDS file
saveRDS(plot, file = rds_filename)