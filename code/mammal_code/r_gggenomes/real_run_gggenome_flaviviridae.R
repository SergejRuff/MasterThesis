# # library(gggenomes)
# 
# GenoInfo <- read.table("gggenome_test.txt", header = TRUE, sep = "\t")
# # GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)
# GenoInfo
# 
# s0 <- GenoInfo[, c("seq_id", "sequence_length")]
# s0 <- s0[!duplicated(s0$seq_id), ]
# s0
# colnames(s0)[colnames(s0) == "sequence_length"] <- "length"
# 
# g0 <- GenoInfo[, c("seq_id", "start", "end", "interpro_description")]
# g0
# 
# 
# gggenomes(seqs = s0, genes = g0) +
#   geom_seq() +
#   geom_bin_label(color="red") +
#   geom_gene(aes(fill = interpro_description))

# ===
library(gggenomes)
setwd("D:/ChongDesktop/chong2022_2025/myPhD/17. Flaviviridae/gggenome")

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

grouping_file <- read.table("fullTree_genus_grouping.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(grouping_file) <- c("seq_id", "genus")

GenoInfo <- merge(GenoInfo, grouping_file, by = "seq_id", all.x = TRUE)
GenoInfo <- GenoInfo[complete.cases(GenoInfo$genus), ]

s0 <- GenoInfo[, c("seq_id", "sequence_length")]
s0 <- s0[order(-s0$sequence_length), ]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- GenoInfo[, c("seq_id", "start", "end", "interpro_description", "genus")]
g0$start <- g0$start * 3
g0$end <- g0$end * 3


# Define the number of rows and columns for each page
rows_per_page <- 60
cols_per_page <- 1

# Calculate the total number of pages needed
total_pages <- ceiling(nrow(s0) / rows_per_page)

pdf("genome_map_flaviviridae_fullTree.pdf", width = 8.27, height = 11.69)

# Color mapping for interpro_description
color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

# Color mapping for genus
genus_color_mapping <- c(
  "Pegivirus" = "#97F8FF",
  "Hepacivirus" = "#B7FFB7",
  "Pestivirus" = "#F599FF",
  "Orthoflavivirus" = "#FF9999"
)

for (page in 1:total_pages) {
  # Subset data for the current page
  start_row <- (page - 1) * rows_per_page + 1
  end_row <- min(page * rows_per_page, nrow(s0))
  current_s0 <- s0[start_row:end_row, , drop = FALSE]
  current_g0 <- g0[g0$seq_id %in% current_s0$seq_id, , drop = FALSE]
  
  plot <- gggenomes(seqs = current_s0, genes = current_g0) +
    geom_seq() +
    geom_gene(aes(fill = interpro_description)) +
    scale_fill_manual(values = color_mapping) +
    geom_bin_label(size = 2) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    theme(legend.position = "none",  axis.text.x = element_text(size = 10)) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))
  
  print(plot)
}

dev.off()


# === Selected === 
# Pestivirus
# NC_025677.1; SRR6478578_Supercontig_1

library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_025677.1", "SRR6478578_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "start", "end", "interpro_description")]
g0$start <- g0$start * 3
g0$end <- g0$end * 3

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

# gggenomes(seqs = s0, genes = g0) +
#   geom_seq() +
#   geom_bin_label(size = 4, color="#E809FF") +
#   geom_gene(aes(fill = interpro_description)) +
#   scale_fill_manual(values = color_mapping) + 
#   theme(legend.position = "none")

p<- gggenomes(seqs = s0, genes = g0) +
  geom_seq() +
  geom_bin_label(size = 4, color="#E809FF") +
  geom_gene(aes(fill = interpro_description)) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none",  axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p

# Hepacivirus
# NC_009827.1; SRR11262332_Supercontig_1

library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_009827.1", "OX394150.1", "SRR6426018_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "start", "end", "interpro_description")]
g0$start <- g0$start * 3
g0$end <- g0$end * 3

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)


p<- gggenomes(seqs = s0, genes = g0) +
  geom_seq() +
  geom_bin_label(size = 4, color="#00E600") +
  geom_gene(aes(fill = interpro_description)) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none",  axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p

# Pegivirus
# NC_001710.1; ERR1331708_Supercontig_1

library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_001710.1", "ERR1331708_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "start", "end", "interpro_description")]
g0$start <- g0$start * 3
g0$end <- g0$end * 3

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

p<- gggenomes(seqs = s0, genes = g0) +
  geom_seq() +
  geom_bin_label(size = 4, color="#00D1D6") +
  geom_gene(aes(fill = interpro_description)) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none",  axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p

# Flavivirus
# NC_001474.2; SRR8076048_Supercontig_4
library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_001474.2", "SRR8076048_Supercontig_4")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "start", "end", "interpro_description")]
g0$start <- g0$start * 3
g0$end <- g0$end * 3

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

p<- gggenomes(seqs = s0, genes = g0) +
  geom_seq() +
  geom_bin_label(size = 4, color="#FF3333") +
  geom_gene(aes(fill = interpro_description)) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none",  axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p




# testing - working as expected (17032024)
library(gggenomes)
library(dplyr)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

# selected_seq_ids <- c("NC_038882.1", "NC_001710.1", "ERR1331708_Supercontig_1", "DRR026823_Supercontig_2")
selected_seq_ids <- c("NC_038882.1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2
g0$end <- g0$end * 3

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=5)+
  geom_rect(aes(xmin=x, xmax=xend, ymin=y+.05, ymax=y+-.05), data=feats(), fill="grey", color='black', linetype = "dotted")+
  geom_seq() +
  geom_bin_label(size = 4, color = "#00D1D6") +
  geom_gene(aes(fill = interpro_description)) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p





# === Pegivirus ===
library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_030291.1", "SRR1758976_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2 + g0$orf_start
g0$end <- g0$end * 3 + g0$orf_start

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0$seq_id), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=9, color="grey")+
  geom_seq() +
  geom_seq_label(size = 4, color = "#00D1D6") +
  geom_bin_label(size = 4, color = "#00D1D6") +
  geom_gene(aes(fill = interpro_description), size = 6) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p


# === Hepacivirus ===
library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_009825.1", "SRR3131110_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2 + g0$orf_start
g0$end <- g0$end * 3 + g0$orf_start

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0$seq_id), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=9, color="grey")+
  geom_seq() +
  geom_seq_label(size = 4, color = "#00E600") +
  geom_bin_label(size = 4, color = "#00E600") +
  geom_gene(aes(fill = interpro_description), size = 6) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p


# === Pestivirus ===
library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_025677.1", "SRR6478578_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2 + g0$orf_start
g0$end <- g0$end * 3 + g0$orf_start

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0$seq_id), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=9, color="grey")+
  geom_seq() +
  geom_seq_label(size = 4, color = "#E809FF") +
  geom_bin_label(size = 4, color = "#E809FF") +
  geom_gene(aes(fill = interpro_description), size = 6) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p


# === Orthoflavivirus ===
library(gggenomes)

GenoInfo <- read.csv("csv_for_genome_map.csv", header = TRUE)

selected_seq_ids <- c("NC_001672.1", "SRR6286363_Supercontig_2")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2 + g0$orf_start
g0$end <- g0$end * 3 + g0$orf_start

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0$seq_id), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=9, color="grey")+
  geom_seq() +
  geom_seq_label(size = 4, color = "#FF3333") +
  geom_bin_label(size = 4, color = "#FF3333") +
  geom_gene(aes(fill = interpro_description), size = 6) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p






# testing 
library(dplyr)
GenoInfo <- read.csv("csv_for_genome_map_testing.csv", header = TRUE)

selected_seq_ids <- c("SRR7138162_Supercontig_1")
selected_rows <- GenoInfo[GenoInfo$seq_id %in% selected_seq_ids, ]

s0 <- selected_rows[, c("seq_id", "sequence_length")]
s0 <- s0[!duplicated(s0$seq_id), ]
colnames(s0)[colnames(s0) == "sequence_length"] <- "length"

g0 <- selected_rows[, c("seq_id", "orf_start", "orf_end", "start", "end", "interpro_description")]
g0$start <- (g0$start * 3) - 2 + g0$orf_start
g0$end <- g0$end * 3 + g0$orf_start

color_mapping <- c(
  "C" = "#F8766D",
  "prM" = "#EF67EB",
  "E" = "#DB8E00",
  "NS1" = "#AEA200", 
  "NS2" = "#64B200", 
  "NS3" = "#00C1A7", 
  "NS4" = "#00BADE", 
  "NS5" = "#0097E2",
  "RdRp" = "#B687FF", 
  "Npro" = "#FF0000"
)

orf0 <- selected_rows[, c("seq_id", "orf_start", "orf_end")]
orf0 <- orf0[!duplicated(orf0), ]
orf0 <- orf0 %>% rename(start = orf_start, end = orf_end)

p <- gggenomes(seqs = s0, genes = g0, feats = orf0) +
  geom_feat(size=9, color="grey")+
  geom_seq() +
  # geom_seq_label(size = 4, color = "#FF3333") +
  geom_bin_label(size = 4, color = "#FF3333") +
  geom_gene(aes(fill = interpro_description), size = 6) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

p



