rm(list=ls())

library(Virusparies)
library(tidyverse)


top_csv <- read.delim("/media/sergej/EXTERNAL_US/git_clone/RscriptsMasterthesis/output/mammals/longest_con_gvtable.tsv")

top_csv <- top_csv %>%
  rename(seq_id = contig_id)


top_csv <- top_csv %>%
  filter(!(seq_id %in% c("SRR17345012_cap3_Contig-13", "SRR17345012_cap3_Contig-11","SRR17345013_cap3_Contig-21",
                         "SRR13364363_cap3_Contig-11","SRR13364364_cap3_Contig-12")))

subjects <- VhgGetSubject(top_csv,groupby = ViralRefSeq_taxonomy,extract_brackets = TRUE)

subjects <- subjects %>% select(-subject_count)


subjects  <- subjects  %>%
  left_join(top_csv %>% select(SRA_run,Phylum,ViralRefSeq_ident,contig_len,ViralRefSeq_taxonomy,host_taxon), by = "ViralRefSeq_taxonomy")


subjects <- subjects %>%
  select(SRA_run,Phylum,ViralRefSeq_taxonomy,Processed_ViralRefSeq_subject,host_taxon,contig_len,ViralRefSeq_ident)


subjects <- subjects[order(subjects$contig_len,decreasing = TRUE), ]

sub_table <- VhgTabularRasa(subjects,title = NULL,names_ = c("SRA experiment","Phylum","Viral reference taxonomy","viral subject","Host","contig length",
                                                "Sequence Identity [%]"))


sub_table <- sub_table%>%
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  # Apply predefined styling
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )


ExportVirusGt(sub_table,filename = "sub_table.png",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)