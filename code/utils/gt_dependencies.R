library(Virusparies)

dep <- read.csv("misc/Vhvg_dependencies/dependencies.csv")

tv <- VhgTabularRasa(dep)

ExportVirusGt(tv,"dependencies.png",path = "misc/Vhvg_dependencies/")