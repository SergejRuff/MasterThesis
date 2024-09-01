rm(list = ls())

# remove old install of Virusparies and install new Github version
remove.packages("Virusparies") # remove old version before installing new
library(remotes)
remotes::install_github("SergejRuff/Virusparies")

# set sys env for brave executable for gtsave 
Sys.setenv(
  CHROMOTE_CHROME = "C:/Program Files/BraveSoftware/Brave-Browser/Application/brave.exe"
)
