#------------------------------------------------------------------------#
# organizing data to modeling (in batch)
# run this will save data into the folder "organized_data", 
# each source will provide a map

# carregar os pacotes
source ('R/packagesR.R')

# open a pallette for maps
par(mfrow=c(2,2),mar=c(0,0,1,1))

# run code to organize and save GBif data
source("data/deteccoes/Gbif/codigos_gbif.R")

# run code to organize and save wikiaves data
source("data/deteccoes/wikiaves/codigos_wiki.R")

# run code to organize and save EBird data
source("data/deteccoes/ebird/codigos_ebird.R")

