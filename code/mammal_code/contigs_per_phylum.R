library(ggplot2)
library(Virusparies)

# Create the data frame with the given data
contig_data <- data.frame(
  Phylum = c("Pisuviricota", "Lenarviricota", "Kitrinoviricota", "Negarnaviricota", "Duplornaviricota"),
  Total_Contigs = c(2553, 773, 552, 281, 124),
  Aligned_Contigs = c(2212, 686, 453, 209, 115),
  Percentage = c("86.64%", "88.75%", "86.78%", "74.37%", "92.74%")
)

# Order the Phylum by Total_Contigs
contig_data$Phylum <- factor(contig_data$Phylum, 
                             levels = contig_data$Phylum[order(contig_data$Total_Contigs, decreasing = TRUE)])

con_per_phylum <- ggplot(contig_data, aes(x = Phylum)) +
  geom_bar(aes(y = Total_Contigs, fill = "Total Contigs"), stat = "identity", width = 0.7, color = "black") +
  geom_bar(aes(y = Aligned_Contigs, fill = "E-value < 1e-5"), stat = "identity", width = 0.7, color = "black") +
  
  # Add the percentage labels above total contigs
  geom_text(aes(y = Total_Contigs, label = paste(Total_Contigs, "(", Percentage, ")", sep = "")), 
            vjust = -0.5, color = "black", size = 3.5) +
  
  # Labels and theme
  labs(x = "Phylum", y = "Number of Contigs") +
  scale_fill_manual(name = "Contig Type", values = c("Total Contigs" = "lightblue", "E-value < 1e-5" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.5, linetype = "dotted", colour = "grey50"),  # Dotted major grid lines
    panel.grid.minor = element_line(linewidth = 0.25, linetype = "dotted", colour = "grey75")   # Dotted minor grid lines
  )+
  scale_y_continuous(limits = c(0, 2600))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ExportVirusPlot(plot = con_per_phylum,"con_per_phylum.png",path = "output/mammals/",
                units = "in",limitsize = FALSE,width = 7,height = 8)