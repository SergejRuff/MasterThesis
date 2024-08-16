# This script is here to filter the hittable Chris sent me.
# I try to extract observations based on the following criteria:
# 1) Treffer mit 'best_query' gegen eine bestimmte Virusgruppe (siehe unten)
# 2) 'num_hits' > 4
# 3) 'ViralRefSeq_E' < 1e-5
# 4) 'ViralRefSeq_ident' < 90
#
# Ausserdem hier ein Ranking nach 'best_query' (also der Reihenfolge nach abarbeiten):
#  
# 1) Flavi_RdRp
# 2) Hepe-Virga_RdRp
# 3) Nido_RdRp
# 4) Negative_Bunya-Arena_RdRp, Negative_Mono-Chu_RdRp, Negative_Orthomyxo_RdRp
#
# after extracting groups bases on the criteria, I analyse them with vhvg.



# --- empty global env --- #
rm(list=ls())

# --- load packages --- #
library(tidyverse)


# --- import data --- #
mammal_raw <- 
  readxl::read_xlsx("data/RNAvirus_Mammals_newJan2023/excel_with_alldata_unfiltered/RNAviruses-Mammals-RNA-newJan2023.xlsx")


# --- filter group based on criteria --- #

Flavi_RdRp <- mammal_raw %>% 
  filter(best_query=="Flavi_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

Hepe_Virga_RdRp <- mammal_raw %>% 
  filter(best_query=="Hepe-Virga_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)


Nido_RdRp <- mammal_raw %>% 
  filter(best_query=="Nido_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

Negative_Bunya_Arena_RdRp <- mammal_raw %>% 
  filter(best_query=="Negative_Bunya-Arena_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

Negative_Mono_Chu_RdRp <- mammal_raw %>% 
  filter(best_query=="Negative_Mono-Chu_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

Negative_Orthomyxo_RdRp <- mammal_raw %>% 
  filter(best_query=="Negative_Orthomyxo_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

Nido_NiRAN <- mammal_raw %>% 
  filter(best_query=="Nido_NiRAN",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

NidoAstro_RdRp <- mammal_raw %>% 
  filter(best_query=="NidoAstro_RdRp",num_hits > 4,ViralRefSeq_E < 1e-5,ViralRefSeq_ident < 90)

# --- extract SRA_run --- #

flavi_sra <- unique(Flavi_RdRp$SRA_run)
hepevirga_sra <- unique(Hepe_Virga_RdRp$SRA_run)
nido_sra <- unique(Nido_RdRp$SRA_run)
bunya <- unique(Negative_Bunya_Arena_RdRp$SRA_run)
mono_chu <- unique(Negative_Mono_Chu_RdRp$SRA_run)
orthomyxo <- unique(Negative_Orthomyxo_RdRp$SRA_run)
Nido_NiRAN <- unique(Nido_NiRAN$SRA_run)
NidoAstro_RdRp <- unique(NidoAstro_RdRp$SRA_run)

# --- export --- #

fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/flavi_sra.txt")
writeLines(flavi_sra, fileConn)
close(fileConn)


fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/hepevirga_sra.txt")
writeLines(hepevirga_sra, fileConn)
close(fileConn)

fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/nido_sra.txt")
writeLines(nido_sra, fileConn)
close(fileConn)

fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/bunya.txt")
writeLines(bunya, fileConn)
close(fileConn)


fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/mono_chu.txt")
writeLines(mono_chu, fileConn)
close(fileConn)


fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/orthomyxo.txt")
writeLines(orthomyxo, fileConn)
close(fileConn)

fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/Nido_NiRAN.txt")
writeLines(Nido_NiRAN, fileConn)
close(fileConn)

fileConn<-file("data/RNAvirus_Mammals_newJan2023/sra_ids/NidoAstro_RdRp.txt")
writeLines(NidoAstro_RdRp, fileConn)
close(fileConn)