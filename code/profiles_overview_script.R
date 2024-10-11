rm(list=ls())


library(Virusparies)
library(gt)


# Import the text files
df_RNA <- read.table("data/profile names/Domain_RNA.txt", header = FALSE, sep = " ")
df_largeDNA <- read.table("data/profile names/Domain_largeDNA.txt", header = FALSE, sep = " ")
df_smallDNA <- read.table("data/profile names/Domain_smallDNA.txt", header = FALSE, sep = " ")


# Add an index column as the first column
df_RNA <- cbind(Index = seq_len(nrow(df_RNA)), df_RNA)
df_largeDNA <- cbind(Index = seq_len(nrow(df_largeDNA)), df_largeDNA)
df_smallDNA <- cbind(Index = seq_len(nrow(df_smallDNA)), df_smallDNA)

df_RNA <- df_RNA[,-c(2:3)]
df_largeDNA <- df_largeDNA[,-c(2:3)]
df_smallDNA <- df_smallDNA[,-c(2:3)]

colnames(df_RNA) <- c("Index","RNA virus profile")
colnames(df_largeDNA) <- c("Index","Large DNA virus profile")
colnames(df_smallDNA) <- c("Index","Small DNA virus Profile")

rnagt <- VhgTabularRasa(df_RNA,title = NULL,cell_colour = "White",col_everyrow = TRUE)
largegt <- VhgTabularRasa(df_largeDNA,title = NULL,cell_colour = "White",col_everyrow = TRUE)
smallgt <- VhgTabularRasa(df_smallDNA,title = NULL,cell_colour = "White",col_everyrow = TRUE)
# 
# rnagt <- rnagt%>%
#   opt_align_table_header(align = "left") %>%
#   opt_stylize(style = 4) %>%  # Apply predefined styling
#   tab_options(
#     data_row.padding = px(2),
#     summary_row.padding = px(3),
#     row_group.padding = px(4),
#     heading.title.font.size = px(20)
#   )
# 
# largegt <- largegt%>%
#   opt_align_table_header(align = "left") %>%
#   opt_stylize(style = 4) %>%  # Apply predefined styling
#   tab_options(
#     data_row.padding = px(2),
#     summary_row.padding = px(3),
#     row_group.padding = px(4),
#     heading.title.font.size = px(20)
#   )
# 
# smallgt <- smallgt%>%
#   opt_align_table_header(align = "left") %>%
#   opt_stylize(style = 4) %>%  # Apply predefined styling
#   tab_options(
#     data_row.padding = px(2),
#     summary_row.padding = px(3),
#     row_group.padding = px(4),
#     heading.title.font.size = px(20)
#   )


ExportVirusGt(rnagt,filename = "rnagt.docx",export_gt_obj = TRUE,path = "data/profile names/",create.dir = TRUE)

ExportVirusGt(largegt,filename = "largegt.docx",export_gt_obj = TRUE,path = "data/profile names/",create.dir = TRUE)

ExportVirusGt(smallgt,filename = "smallgt.docx",export_gt_obj = TRUE,path = "data/profile names/",create.dir = TRUE)