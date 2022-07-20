
#### 1) Instalar pacotes ####
#install.packages(c("here", "doBy", "reshape", "rgdal", "raster",
#                   "sp", "maps", "maptools", "rgeos", "ggplot2",
#                   "gridExtra", "grid", "lattice", "jagsUI", "R2WinBUGS",
#                   "corrplot", "vegan", "unmarked", "pROC", "dplyr"), dependencies = T)

# Uma vez instalados os pacotes, precisamos carregar cada um:

#### 2) carregar os pacotes ####
# pacote para transitar entre pastas
require(here)

# pacote para trabalhar com tabelas dinamicas (como excel)
require(doBy)
require (reshape)
require(dplyr)

## pacote para trabalhar com dados espaciais (shapefiles, rasters)
require(rgdal)
require(raster)
library(sp)
library(maps)
library(maptools)
require(rgeos)

# pacote para plots
require(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

# pacotes para carregar output MCMC
require(jagsUI)
require(R2WinBUGS)

# pacote para correlation plot
require(corrplot)

# pacote para transformacao e padronizacao de dados
require(vegan)

# modelos hierarquicos no unmarked
require(unmarked)

# curvas AUC
require(pROC)
