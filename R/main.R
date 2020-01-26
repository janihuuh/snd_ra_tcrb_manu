
me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/ra_sn_tcrb/"))

require(ggplot2)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(data.table)

## Run all fun_* codes
for(code in list.files("src/sn_ra/R/", "fun", full.names = T, recursive = T)){

  print(code)
  source(code)

}

## Define visualizations
theme_set(theme_classic())

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette  <- colorRampPalette(brewer.pal(9, "Pastel1"))
facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = "black"))
