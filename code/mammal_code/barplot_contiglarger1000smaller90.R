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

combined_ga <- VhgPreprocessTaxa(combined_ga,taxa_rank = "Family")

combined_ga <- VhgAddPhylum(combined_ga,"ViralRefSeq_taxonomy")

combined_ga <- VhgSubsetHittable(combined_ga,ViralRefSeq_E_criteria = 1e-5)


# Define the phylum values to keep
target_phyla <- c("Lenarviricota", "Kitrinoviricota", "Pisuviricota", "Duplornaviricota", "Negarnaviricota")

# Filter combined_ga to keep only rows with specified Phylum values
filtered_data <- combined_ga %>% filter(Phylum %in% target_phyla)

# Subset 1: contig_len > 1000
subset1 <- filtered_data %>% filter(contig_len > 1000)

# Subset 2: contig_len > 1000 and ViralRefSeq_ident < 90
subset2 <- filtered_data %>% filter(contig_len > 1000, ViralRefSeq_ident < 90)

# Create new columns to label the subsets for plotting
subset1 <- subset1 %>% mutate(group = "> 1000 nt")
subset2 <- subset2 %>% mutate(group = "> 1000 nt & < 90% Identity")

# Combine subsets
combined_data <- bind_rows(subset1, subset2)

# Calculate the total count for contig_len > 1000 for each Phylum
total_counts <- subset1 %>%
  group_by(Phylum) %>%
  summarise(total_count = n())

# Calculate the percentage for subset2
subset2_counts <- subset2 %>%
  group_by(Phylum) %>%
  summarise(count = n()) %>%
  left_join(total_counts, by = "Phylum") %>%
  mutate(percentage = round((count / total_count) * 100, 1),  # Percentage calculation
         label = paste0(total_count, " (", percentage, "%)"))  # Combined label for display

# Order Phylum by total count of contig_len > 1000 in descending order
combined_data$Phylum <- factor(combined_data$Phylum, levels = total_counts %>%
                                 arrange(desc(total_count)) %>%
                                 pull(Phylum))

# Calculate counts and percentages for specific phyla with contig_len > 1000 and ViralRefSeq_ident < 90
target_phyla_counts <- subset2 %>%
  filter(Phylum %in% c("Pisuviricota", "Kitrinoviricota", "Lenarviricota")) %>%
  summarise(count = n()) %>%
  pull(count)

# Calculate the percentage in comparison to 444
target_percentage <- round((target_phyla_counts / 444) * 100, 1)

# Combine text for annotation
annotation_text <- paste0("Pisuviricota, Kitrinoviricota, and Lenarviricota:\n", 
                          target_phyla_counts, " contigs (", target_percentage, "% of 444)")


# Plot
con_plot <- ggplot(combined_data, aes(x = Phylum, fill = group)) +
  geom_bar(position = "identity", alpha = 0.7) +  # Overlapping bars with transparency
  scale_fill_manual(values = c("> 1000 nt" = "lightblue",
                               "> 1000 nt & < 90% Identity" = "steelblue")) +
  # Add combined count and percentage text as number (%)
  geom_text(data = subset2_counts, aes(x = Phylum, y = total_count, label = label),
            vjust = -0.5, color = "black", inherit.aes = FALSE) +
  labs(x = "Phylum", y = "Count", fill = "Contig Criteria") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(linewidth = 0.5, linetype = "dotted", colour = "grey50"),  
        panel.grid.minor = element_line(linewidth = 0.25, linetype = "dotted", colour = "grey75"))  +  # Flip coordinates for horizontal bars
  # Add annotation box at the top right
  annotate("label", x = Inf, y = Inf, label = annotation_text, 
           hjust = 1.1, vjust = 1.1, color = "black", 
           fill = "white", size = 3.5,              # Adjust text size
           label.size = 0.4,
           label.padding = unit(c(1, 1, 1, 1), "lines"),  # Increase padding
           fill = "lightgray")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 

# Save the plot to output/mammals/
ggsave("output/mammals/contig_analysis_plot.png", plot = con_plot, width = 10, height = 8, dpi = 300)