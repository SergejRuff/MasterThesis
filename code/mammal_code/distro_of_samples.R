rm(list=ls())


library(Virusparies)
library(gt)



mammalian_dis <- data.frame(
  Host = c("Bos taurus (Cattle)", "Sus scrofa domesticus (Domestic Pig)", "Chlorocebus sabaeus (Green monkey)", 
              "Sus scrofa (Wild Boar)", "Otomops harrisoni (Harrison's large-eared giant mastiff bat)"),
  Num_of_contigs = c(2401, 1073, 1059, 317, 186),
  Percentage = c(38.67, 17.28, 17.06, 5.11, 3.00)
)


mammal_gt <- VhgTabularRasa(mammalian_dis,title = NULL,names_ = c("Host","Number of contigs","Percentage"))


mammal_gt <- mammal_gt %>% 
  opt_align_table_header(align = "left") %>%
  opt_stylize(style = 4) %>%  
  tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3),
    row_group.padding = px(4),
    heading.title.font.size = px(20)
  )%>%
  # Rename the third column by adding an asterisk
  cols_label(
    Percentage = "Percentage*"
  )  %>%
  tab_footnote(
    footnote = "*Percentage relative to the total number of contigs (6209).",
   
  )


ExportVirusGt(mammal_gt,filename = "mammal_distro.png",export_gt_obj = TRUE,path = "output/mammals/statistics",create.dir = TRUE)