rm(list = ls())

"

this script is meant for the generation of timelines used in my
defense - see powerpoint prasentation

page:


"


## load libraries
library(timevis)
library(timelineS)
library(ggplot2)
library(ggrepel)
library(ggbreak)
library(viridis)

# generate data for plot

"
First Plot is suppose to show a timeline with all pandemics and epedemics
caused by viruses

source used for data:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8525686/ - Pandemics Throughout the History
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8483091/ - Viral Pandemics in the Past Two Decades: An Overview

"


pandem_time <- data.frame(
  id = 1:12,
  content = c("Russian flu","Spanish flu","Asian flu",
              "Hong Kong flu","HIV/AIDS","SARS CoV","Swine flu","MERS","Ebola","Zika Virus","SARS CoV-2 (COVID19)","Mpox epidemic"),
  start = paste0(c("1889","1918","1957","1968","1981","2002","2009","2012","2014","2015","2019","2023"), "-01-01"),
  end = c("1890-12-31","1920-12-31","1958-12-31",NA,NA,"2004-12-31","2010-12-31",NA,"2016-12-31",NA,NA,NA)
)

pandem_time$year <- as.numeric(sub("^(\\d{1,4})-.*$", "\\1", pandem_time$start))
######### Generate plot


# Set the heights we will use for our milestones.
positions <- c(0.5, -0.5, 1.0, -1.0, 1.25, -1.25, 1.5, -1.5)

# Set the directions we will use for our milestone, for example above and below.
directions <- c(1, -1)

line_pos <- data.frame(
  "start"=unique(pandem_time$start),
  "position"=rep(positions, length.out=length(unique(pandem_time$start))),
  "direction"=rep(directions, length.out=length(unique(pandem_time$start))))


pandem_time <- merge(x=pandem_time, y=line_pos, by="start", all = TRUE)

pandem_time$start <- as.Date(pandem_time$start, format = "%Y-%m-%d")
pandem_time$end <- as.Date(pandem_time$end, format = "%Y-%m-%d")

# Create a year dataframe with a sequence of years from 100 to 2030
year_df <- data.frame(
  year_date_range = seq(1880, 2030, by = 20),  # Numeric years
  year_format = as.character(seq(1880, 2030, by = 20))  # Year labels
)

pandem_time$start <- as.Date(pandem_time$start)

# Add a new column for colors based on the start year
pandem_time <- pandem_time %>%
  mutate(color = case_when(
    start < as.Date("1900-01-01") ~ "19th Century",  # Color for years before 1900
    start >= as.Date("1900-01-01") & start < as.Date("2000-01-01") ~ "20th Century",  # Color for years between 1900 and 2000
    start >= as.Date("2000-01-01") ~ "21st Century"     # Color for years after 2000
  ))


# Lets offset the labels 0.2 away from scatter points
text_offset <- 0.2 

# Let's use the absolute value since we want to add the text_offset and increase space away from the scatter points 
absolute_value<-(abs(pandem_time$position)) 
text_position<- absolute_value + text_offset

# Let's keep the direction above or below for the labels to match the scatter points
pandem_time$text_position<- text_position * pandem_time$direction 

# Create the initial timeline plot with specified colors
timeline_plot <- ggplot(pandem_time, aes(x = as.numeric(format(start, "%Y")),position, label = content, color = color)) + 
  theme_classic() + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) +
  geom_segment(aes(yend = 0, xend = as.numeric(format(start, "%Y"))), size = 0.2) +
  geom_point(size = 3) +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom") +
  geom_text(data = pandem_time, 
            aes(x = as.numeric(format(start, "%Y")), y = text_position, label = paste0(content,"\n",year)), 
            size = 3.5, 
            color = 'black') +  # Adjust y position for text labels
  labs(color = "Century") +  # Add title for the legend
  guides(color = guide_legend(override.aes = list(size = 5)))  # Ensure legend shows larger points

# Create the legend data for custom color mapping
legend_data <- data.frame(
  label = c("19th Century", "20th Century", "21st Century"),
  color = c("darkorange", "deepskyblue", "seagreen")
)

# Update the pandem_time dataframe to map colors correctly
pandem_time$color <- factor(pandem_time$color, levels = legend_data$label)

# Update the timeline_plot to use scale_color_manual
timeline_plot <- timeline_plot + 
  scale_color_manual(values = setNames(legend_data$color, legend_data$label))

timeline_plot <- timeline_plot + geom_text(data=year_df, aes(x=year_date_range,y=-0.1,label=year_format),size=3.5, color='black') 

timeline_plot <- timeline_plot + ggtitle("Major Viral Pandemics and Epidemics Since the Industrial Era")+ 
  theme(plot.title = element_text(
    hjust = 0.5,                  # Center-align the title
    size = 16,                   # Font size
    face = "bold",               # Bold font
    color = "#4B4B4B",              # Title color changed to blue
    family = "Arial",             # Font family
    margin = margin(b = 10)      # Bottom margin for spacing
  ))
  
# Print the plot
print(timeline_plot)