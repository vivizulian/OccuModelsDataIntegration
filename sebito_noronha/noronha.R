# load data gbif

#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################


library(rgbif)
library(scrubr)
library(maps)
library(maditr)


# load noronha shapefile (world heritage sites)
require(rgdal)
require(here)
require(raster)
require(ggplot2)
library(rgeos)
library(maptools)
library(ggsn)

#noronha <- readOGR(here("sebito_noronha", "FernandoNoronha.shp"))
noronha <- readOGR("FernandoNoronha.shp")
plot(noronha)

#noronha <- fortify(noronha)

#crs(noronha) <- crs("+proj=longlat +datum=WGS84 +no_defs")

#hs <- readOGR(here("WorldMarineHeritageSites"), "2013_02_WorldHeritageMarineProgramme")
#noronha<- hs [grep ("Fernando de Noronha",hs$Subarea),]


#Make a grid over the study area
make_grid <- function(x, cell_diameter, cell_area, clip = TRUE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons") #define the extension of the grid
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.05, 0.05))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter) 
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE) #intersect and clip the grid using the study area
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

hex <- make_grid(noronha, cell_diameter = 0.005) #cell diameter-you can change to bigger or smaller
crs(hex) <- crs(noronha)

plot(noronha)
plot(hex, add=T)



# studied species
myspecies <- c("Vireo gracilirostris", "Vireo olivaceus gracilirostris", "Vireo gracilirostris Sharpe, 1890")


# help : https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/
# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_data <- occ_data (locality = "Fernando de Noronha", hasCoordinate = TRUE, limit = 20000)#(scientificName = myspecies, hasCoordinate = TRUE, limit = 20000)
# get the DOIs for citing these data properly:
gbif_citation(gbif_data)

# avian data
require(dplyr)
gbif_data <- gbif_data$data %>% filter (class == "Aves")

# gbif
gbif_data <- gbif_data %>% #select (decimalLatitude, decimalLongitude,scientificName) %>% 
  rename (decimallatitude = decimalLatitude ,
          decimallongitude = decimalLongitude,
          scientificname = scientificName) %>%
  mutate (base = "gbif")


# records within noronha
sp_points_gbif <- data.frame (x = gbif_data$decimallongitude,
                                 y = gbif_data$decimallatitude)

# into sp points
# pointss
coordinates(sp_points_gbif) <- ~x + y
crs (sp_points_gbif) <- crs (noronha)

# rasterize to count the number of papers per cell
overlap_pts_gbif <- over(sp_points_gbif, noronha)

# filter
gbif_data <- gbif_data[is.na(overlap_pts_gbif$NM_MUNICIP) !=T,]

# ------------------------------------------------

library(rinat)
# help
# https://www.r-bloggers.com/2021/12/mapping-inaturalist-data-in-r/

inat_obs_df <- get_inat_obs(bounds = noronha,
                            taxon_name = "Aves",
                              maxresults = 10000)

# change names

inat_obs_df <- inat_obs_df %>% #  select (latitude, longitude,scientific_name) %>% 
  rename (decimallatitude = latitude ,
          decimallongitude = longitude,
          scientificname = scientific_name) %>%
  mutate (base = "inat")


# records within noronha
sp_points_inat <- data.frame (x = inat_obs_df$decimallongitude,
                                 y = inat_obs_df$decimallatitude)

# into sp points
# pointss
coordinates(sp_points_inat) <- ~x + y
crs (sp_points_inat) <- crs (noronha)

# rasterize to count the number of papers per cell
overlap_pts_inat <- over(sp_points_inat, noronha)

# filter
inat_obs_df<-inat_obs_df[is.na(overlap_pts_inat$NM_MUNICIP) !=T,]



# ------------------------------------
  
# vertnet

require(rvertnet)

vert_data <- searchbyterm(stateprovince = "pernambuco", limit = 20000)
vert_data_local <- vert_data$data[is.na(vert_data$data$decimallatitude) != T,]
vert_data_local$decimallongitude <- as.numeric(vert_data_local$decimallongitude)
vert_data_local$decimallatitude <- as.numeric(vert_data_local$decimallatitude)

# avian data
vert_data_local <- vert_data_local %>% filter (class == "Aves")

# mutate data
vert_data_local <- vert_data_local %>% 
  mutate (base = "vertnet")


# records within noronha
sp_points_vertnet <- data.frame (x = vert_data_local$decimallongitude,
                               y = vert_data_local$decimallatitude)

# into sp points
# pointss
coordinates(sp_points_vertnet) <- ~x + y
crs (sp_points_vertnet) <- crs (noronha)

# rasterize to count the number of papers per cell
overlap_pts_vertnet <- over(sp_points_vertnet, noronha)

# filter
vert_data_local<-vert_data_local[is.na(overlap_pts_vertnet$NM_MUNICIP) !=T,]


# ---------------------- #_---------------------------

# load eBird data (downloaded by V Zulian)
ebird_data <- read.csv ("ebd_BR-PE_relJun-2022.csv", sep = ";")

#eBird data processing:
#read table BR
ebBR <- read.delim("ebd_BR-PE_relJun-2022.txt", header=T, quote="")
ebBR <- ebBR[,-50]

### Filtering data of Vireo gracilirostris 

recordsBR <- subset(ebBR, SCIENTIFIC.NAME %in% c("Vireo gracilirostris")) #check number of vireo records #easy to change the sp

ebBR <- within(ebBR, ID <- NA) #create a collum to put the unique sample event number
ebBR <- transform(ebBR, ID=match(SAMPLING.EVENT.IDENTIFIER, unique(SAMPLING.EVENT.IDENTIFIER))) #put the same number for each unique sample event

ebBR <- within(ebBR, NSPECIES <- 1) #add a column to put the number of species of each list
ebBR <- within(ebBR, Vireo_gracilirostris <- 0) #add a column to put 0 or 1 to Vireo in each list

vr.grac <- subset(ebBR, SCIENTIFIC.NAME %in% c("Vireo gracilirostris")) #extract the recorders of the species
all.else <- subset(ebBR, SCIENTIFIC.NAME != "Vireo gracilirostris")

for (i in 1:nrow(vr.grac)) {   #put 1 on line of the species record
  vr.grac$Vireo_gracilirostris[i] = 1
}

ebBR <- rbind(vr.grac,all.else) #join both tables again

ebBR2 <- ebBR[!duplicated(ebBR$SAMPLING.EVENT.IDENTIFIER),] #extract the unique lists - the same length that ID
dim(ebBR2[ebBR2$Vireo_gracilirostris==1,]) #just to check

test2 <- aggregate(ebBR$NSPECIES, by=list(ebBR$ID), FUN=sum) #summ of the number of species by list

ebBR2$NSPECIES <- test2$x[match(ebBR2$ID, test2$Group.1)] #join the total number of species by list

virgrac <- ebBR2

save(virgrac, file="VIRGRAC_DetNDet.RData")

#150 registros com protocolo Traveling, 71 Stationary, 32 como Historical e 30 como Incidental


load("VIRGRAC_DetNDet.RData")
ebird_data <- virgrac

# records within noronha
sp_points_ebird <- data.frame(x = ebird_data$LONGITUDE,
                               y = ebird_data$LATITUDE )

# into sp points
require(sp)
require(rgeos)

# pointss
coordinates(sp_points_ebird) <- ~x + y
crs(sp_points_ebird) <- crs(noronha)

# rasterize to count the number of papers per cell
overlap_pts_ebird <- over(sp_points_ebird, noronha)

# filter
ebird_data <- ebird_data[which(is.na(overlap_pts_ebird$CD_GEOCMU)!= T),]

# ebird
ebird_data <- ebird_data %>% #select (LATITUDE, LONGITUDE, SCIENTIFIC.NAME) %>% 
  rename (decimallatitude = LATITUDE ,
          decimallongitude = LONGITUDE,
          scientificname = SCIENTIFIC.NAME) %>%
  mutate (base = "eBird")



# --------------------------------------------------------------------



# create a grid to have sites
res_raster<- 0.01 # latlong degrees

# based on the extent of extracted data
grd_df <- expand.grid(x = seq(from = (extent(noronha))[1]-0.01,                                           
                              to = (extent(noronha))[2]+0.01,
                              by = res_raster),
                      y = seq(from = (extent(noronha))[3]-0.01,                                           
                              to = (extent(noronha))[4]+0.01,
                              by = res_raster))  # expand points to grid

# points object
coordinates(grd_df) <- ~x + y

# Sp points into raster
grd_raster <- (raster(grd_df,resolution = res_raster))

crs(grd_raster) <- crs(noronha)

# define values
values(grd_raster) <- rnorm (length(values(grd_raster)))


# load elevation map
# https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
require(elevatr)
elevation <- get_elev_raster(noronha, z = 9,expand=1)

# find site elevation
cont.elev.sp <- raster::extract(elevation, 
                                coordinates(hex), 
                                fun=mean,
                                cell=T) 
# the sites
rownames(cont.elev.sp) <- seq (1,nrow(cont.elev.sp))


# all site of analysis (all raster cells)
all_sites <- rownames(cont.elev.sp)

# maps
plot(hex, add=T)
plot(noronha)
plot(elevation, add=T)
plot(noronha, add=T)
points(inat_obs_df[ , c("decimallongitude", "decimallatitude")], pch = 19, col = "red")
points(gbif_data[ , c("decimallongitude", "decimallatitude")], pch = 19, col = "black")
points(vert_data_local[ , c("decimallongitude", "decimallatitude")], 
       pch = 19, col = "green")
points(ebird_data[ , c("decimallongitude", "decimallatitude")], 
       pch = 19, col = "cyan")



# --------------------------------------------------- #

overlap_sites <- lapply(list(inat_obs_df, gbif_data, vert_data_local, ebird_data), function (i) {
      
    # overlap with noronha's shape to have the sites
    # sp points
    sp_points <- data.frame(x = as.numeric(i$decimallongitude), 
                             y = as.numeric(i$decimallatitude))
    # pointss
    coordinates(sp_points) <- ~x + y
    crs(sp_points) <- crs(hex) # crs
    
    # rasterize to count the number of papers per cell
    overlap_pts <- raster::extract(hex, sp_points, cell=TRUE )
    
    # bind the cell identity (our site)
    sites <- overlap_pts[,"poly.ID"]
    length(unique(overlap_pts[,"poly.ID"])) # nsites
    
    # find points over noronha
    noronha_points <- over(sp_points,noronha)
    overlap_noronha <- noronha_points$Subarea
    ;
    sites
    
})


# bind sites in the datasets
inat_obs_df$sites <- overlap_sites[[1]]
gbif_data$sites <- overlap_sites[[2]]
vert_data_local$sites <- overlap_sites[[3]]
ebird_data$sites <- overlap_sites[[4]]



# --------------------------------------------------------------# 
# species per site - iNaturalist

#iNaturalist has records from 2003 to 2022 - think about the time scale!

inat_obs_df <- within(inat_obs_df, NSPECIES <- 1) #add a column to put the number of species per grid
inat_obs_df <- within(inat_obs_df, Vireo_gracilirostris <- 0) #add a column to put 0 or 1 to Vireo in each list

inat_obs_df$Vireo_gracilirostris <- ifelse(inat_obs_df$scientificname=="Vireo gracilirostris", 1, 0)

inat_obs_df2 <- inat_obs_df[!duplicated(inat_obs_df$sites),] #extract the unique lists - the same length that ID

dim(inat_obs_df2[inat_obs_df2$Vireo_gracilirostris==1,]) #just to check

nspec <- aggregate(inat_obs_df$NSPECIES, by=list(inat_obs_df$sites), FUN=sum) #summ of the number of species by list
vir <- aggregate(inat_obs_df$Vireo_gracilirostris, by=list(inat_obs_df$sites), FUN=sum)

inat_obs_df2$NSPECIES <- nspec$x[match(inat_obs_df2$sites, nspec$Group.1)] #join the total number of species by list
inat_obs_df2$Vireo_gracilirostris <- vir$x[match(inat_obs_df2$sites, vir$Group.1)] #join the detections



# --------------------------------------------------------------# 
# species per site - GBIF

#iNaturalist has records from 2008 to 2018 - think about the time scale!

gbif_data <- within(gbif_data, NSPECIES <- 1) #add a column to put the number of species per grid
gbif_data <- within(gbif_data, Vireo_gracilirostris <- 0) #add a column to put 0 or 1 to Vireo in each list

gbif_data$Vireo_gracilirostris <- ifelse(gbif_data$scientificname=="Vireo gracilirostris Sharpe, 1890", 1, 0)

gbif_data2 <- gbif_data[!duplicated(gbif_data$sites),] #extract the unique lists - the same length that ID

dim(gbif_data2[gbif_data2$Vireo_gracilirostris==1,]) #just to check

nspec <- aggregate(gbif_data$NSPECIES, by=list(gbif_data$sites), FUN=sum) #summ of the number of species by list
vir <- aggregate(gbif_data$Vireo_gracilirostris, by=list(gbif_data$sites), FUN=sum)

gbif_data2$NSPECIES <- nspec$x[match(gbif_data2$sites, nspec$Group.1)] #join the total number of species by list
gbif_data2$Vireo_gracilirostris <- vir$x[match(gbif_data2$sites, vir$Group.1)] #join the detections
gbif_data2$Vireo_gracilirostris <- ifelse(gbif_data2$Vireo_gracilirostris>0, 1,0)



# --------------------------------------------------------------# 
# species per site - VertNet

#VetNet only has only 5 records. All 5 from V. gracilirostris and from 1973.
#I don't think that make sense to use this data set. 



# --------------------------------------------------------------# 
# species per site - eBird - It's ready in the object ebird_data

# Filter the lists with Stationary and Traveling protocols:

ebird_data <- ebird_data[which(ebird_data$PROTOCOL.TYPE == "Traveling" 
                        | ebird_data$PROTOCOL.TYPE == "Stationary"),]

ebird_data$EFFORT.DISTANCE.KM <- ifelse(is.na(ebird_data$EFFORT.DISTANCE.KM), 0, ebird_data$EFFORT.DISTANCE.KM)



# ----------------------------------------------------------------------

# elevation (site occupancy variable)

modeling_data <- list(site_covs = cont.elev.sp,
                       ebird_detection = ebird_data$Vireo_gracilirostris,
                       ebird_duration = ebird_data$DURATION.MINUTES,
                       ebird_distance = ebird_data$EFFORT.DISTANCE.KM,
                       ebird_nsp = ebird_data$NSPECIES,
                       inat_detection = inat_obs_df2$Vireo_gracilirostris,
                       inat_effort = inat_obs_df2$NSPECIES,
                       gbif_detection = gbif_data2$Vireo_gracilirostris,
                       gbif_effort = gbif_data2$NSPECIES,
                       EBsite = ebird_data$sites,
                       GBsite = gbif_data2$sites,
                       INsite = inat_obs_df2$sites, site = noronha,
                       grid = hex)
                       


# save 
save (modeling_data,
      file = "modeling_dataVZ.RData")






#####################################
#Data processing Andre

# species per site
require("data.table")
inat_obs_df <- setDT(inat_obs_df)
gbif_data <- setDT(gbif_data)
vert_data_local <- setDT(vert_data_local)

data_modeling <- lapply(list(inat_obs_df, gbif_data, vert_data_local), function (i) {

  # pivot table
  detection_matrix <- dcast(sites~scientificname,
                           data = inat_obs_df,
                           drop=F)
    
  # imputing missing sites  
  to_input <-data.frame(matrix(NA, 
                                nrow = length(all_sites[which(as.numeric(all_sites) %in% (detection_matrix$sites) !=T)]), # cells not in the table
                                ncol = ncol (detection_matrix)))
  colnames(to_input) <- colnames(detection_matrix)        
  to_input$sites <- all_sites[which(all_sites %in% detection_matrix$sites !=T)] 
  
  # bind & order
  detection_matrix_complete <- rbind (detection_matrix, to_input) # bind
  detection_matrix_complete <- detection_matrix_complete[order(as.numeric(detection_matrix_complete$sites)),]# order
  
  # detection matrix of my species
  focus_sp_mat <- detection_matrix_complete[,which(colnames(detection_matrix_complete) %in% myspecies)]
  focus_sp_mat <- ifelse(focus_sp_mat>1,1,focus_sp_mat)
  
  # effort matrix
  effort_mat <- detection_matrix_complete[,-which(colnames(detection_matrix_complete) == "sites")]
  effort_mat[effort_mat>1] <- 1
  effort_mat <- apply(effort_mat,1,sum,na.rm=T)
  
  # list of data
  res <- list(det = focus_sp_mat,
               effort = effort_mat)
  res
  
})
  

# data modeling ebird

# pivot table for each species
# all_lists <- unique(ebird_data$SAMPLING.EVENT.IDENTIFIER)
detection_matrix <- dcast(sites~SAMPLING.EVENT.IDENTIFIER,
                          data = ebird_data[which(ebird_data$scientificname %in% myspecies
                                                  &
                                                  ebird_data$PROTOCOL.TYPE == "Traveling" 
                                                  | ebird_data$PROTOCOL.TYPE == "Stationary"),],
                           drop=F,
                           value = "scientificname")

# imput sites & lists
to_input <- data.frame(matrix(NA, 
                               nrow = length(all_sites[which(as.numeric(all_sites) %in% (detection_matrix$sites) !=T)]), # cells not in the table
                               ncol = ncol (detection_matrix))) # sites
colnames(to_input) <- colnames(detection_matrix)    # set colnames
to_input$sites <- all_sites[which(all_sites %in% detection_matrix$sites !=T)]  # set site names

# bind & order
detection_matrix_complete <- rbind(detection_matrix, to_input) # bind
detection_matrix_complete <- detection_matrix_complete[order(as.numeric(detection_matrix_complete$sites)),]# order


# -----------

# To think!: should we input lists [I don't think so]? There are 199 lists, of which 112 in V. graci dataset

# ------------

# effort data
# duration 
ebird_data$DURATION.MINUTES <- as.numeric(ebird_data$DURATION.MINUTES)
duration_matrix <- dcast (sites~SAMPLING.EVENT.IDENTIFIER,
                          data = ebird_data[which(ebird_data$scientificname %in% myspecies
                                                  &
                                                    ebird_data$PROTOCOL.TYPE != "Historical"),], # my focus spp / avoiding historical data
                          drop=F,
                          value.var = "DURATION.MINUTES",
                          fun.aggregate = sum,
                          na.rm=T)
# imput sites & lists
to_input <-data.frame (matrix (NA, 
                               nrow = length(all_sites[which(as.numeric(all_sites) %in% (duration_matrix$sites) !=T)]), # cells not in the table
                               ncol = ncol (duration_matrix))) # sites
colnames(to_input) <- colnames(duration_matrix)         # set colnames
to_input$sites <- all_sites[which(all_sites %in% duration_matrix$sites !=T)]  # set site names

# bind & order
duration_matrix <- rbind (duration_matrix, to_input) # bind
duration_matrix  <- duration_matrix[order(as.numeric(duration_matrix$sites)),]# order
duration_matrix[is.na(duration_matrix)] <- 0 # no effort == zero duration



# duration 
ebird_data$EFFORT.DISTANCE.KM <- as.numeric(ebird_data$EFFORT.DISTANCE.KM)
distance_matrix <- dcast (sites~SAMPLING.EVENT.IDENTIFIER,
                          data = ebird_data[which(ebird_data$scientificname %in% myspecies
                                                  &
                                                    ebird_data$PROTOCOL.TYPE != "Historical"),],
                          drop=F,
                          value.var = "EFFORT.DISTANCE.KM",
                          fun.aggregate = sum,
                          na.rm=T)
# imput sites & lists
to_input <-data.frame (matrix (NA, 
                               nrow = length(all_sites[which(as.numeric(all_sites) %in% (distance_matrix$sites) !=T)]), # cells not in the table
                               ncol = ncol (distance_matrix))) # sites
colnames(to_input) <- colnames(distance_matrix)         # set colnames
to_input$sites <- all_sites[which(all_sites %in% distance_matrix$sites !=T)]  # set site names

# bind & order
distance_matrix <- rbind (distance_matrix, to_input) # bind
distance_matrix  <- distance_matrix[order(as.numeric(distance_matrix$sites)),]# order
distance_matrix[is.na(distance_matrix)] <- 0 # no effort == zero duration


# --------------------------------------------------------------------------------

# historical data

detection_matrix_historical <- dcast (sites~SAMPLING.EVENT.IDENTIFIER,
                           data = ebird_data[which(ebird_data$scientificname %in% myspecies
                                                   &
                                                     ebird_data$PROTOCOL.TYPE == "Historical"),],
                           drop=F,
                           value = "scientificname")

# imput sites & lists
to_input <-data.frame (matrix (NA, 
                               nrow = length(all_sites[which(as.numeric(all_sites) %in% (detection_matrix_historical$sites) !=T)]), # cells not in the table
                               ncol = ncol (detection_matrix_historical))) # sites
colnames(to_input) <- colnames(detection_matrix_historical)         # set colnames
to_input$sites <- all_sites[which(all_sites %in% detection_matrix_historical$sites !=T)]  # set site names

# bind & order
detection_matrix_historical <- rbind (detection_matrix_historical, to_input) # bind
detection_matrix_historical  <- detection_matrix_historical[order(as.numeric(detection_matrix_historical$sites)),]# order

# det historical
historical_frequency <- apply (is.na(detection_matrix_historical[,-1])!=T,1,sum,na.rm=T) 
historical_frequency<-ifelse(historical_frequency>0,1,historical_frequency)

# ----------------------------------------------------------------------

# elevation (site occupancy variable)

modeling_data <- list (site_covs = cont.elev.sp,
                       hist_data_covs = historical_frequency,
                       hist_data_det = detection_matrix_historical,
                       ebird_detection = detection_matrix,
                       ebird_duration = duration_matrix,
                       ebird_distance = distance_matrix,
                       inat_detection = data_modeling[[1]]$det,
                       inat_effort = data_modeling[[1]]$effort,
                       gbif_detection = data_modeling[[2]]$det,
                       gbif_effort = data_modeling[[2]]$effort,
                       vertnet_detection = data_modeling[[3]]$det,
                       vertnet_effort = data_modeling[[3]]$effort,
                       # map
                       noronha,
                       # grid
                       grd_raster)
                       

# save 
save (modeling_data,
      file = "modeling_data.RData")


