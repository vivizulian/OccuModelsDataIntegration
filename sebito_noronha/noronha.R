# load data gbif

#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################


library(rgbif)
library(scrubr)
library(maps)



# load noronha shapefile (world heritage sites)
require(rgdal)
require(here)
require(raster)
hs <- readOGR(here ("WorldMarineHeritageSites"), "2013_02_WorldHeritageMarineProgramme")
noronha<- hs [grep ("Fernando de Noronha",hs$Subarea),]

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



# ------------------------------------------------

library(rinat)
# help
# https://www.r-bloggers.com/2021/12/mapping-inaturalist-data-in-r/

inat_obs_df <- get_inat_obs(bounds= noronha,
                            taxon_name = "Aves",
                              maxresults = 10000)

# change names

inat_obs_df <- inat_obs_df %>% #  select (latitude, longitude,scientific_name) %>% 
  rename (decimallatitude = latitude ,
          decimallongitude = longitude,
          scientificname = scientific_name) %>%
  mutate (base = "inat")



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
crs (sp_points_vertnet) <- crs ("+proj=longlat +datum=WGS84 +no_defs")

# rasterize to count the number of papers per cell
overlap_pts_vertnet <- over(sp_points_vertnet,noronha)

# filter
vert_data_local<-vert_data_local[is.na(overlap_pts_vertnet$Full_Name) !=T,]


# ---------------------- #_---------------------------

# load eBird data (downloaded by V Zulian)
ebird_data <- read.csv ("ebd_BR-PE_relJun-2022.csv", sep = ";")

# records within noronha
sp_points_ebird <- data.frame (x = ebird_data$LONGITUDE,
                               y = ebird_data$LATITUDE )

# into sp points
require(sp)
require(rgeos)

# pointss
coordinates(sp_points_ebird) <- ~x + y
crs (sp_points_ebird) <- crs ("+proj=longlat +datum=WGS84 +no_defs")

# rasterize to count the number of papers per cell
overlap_pts_ebird <- over(sp_points_ebird,noronha)

# filter
ebird_data <- ebird_data[which(is.na(overlap_pts_ebird$Latitude) != T),]

# ebird
ebird_data <- ebird_data %>% #select (LATITUDE, LONGITUDE, SCIENTIFIC.NAME) %>% 
  rename (decimallatitude = LATITUDE ,
          decimallongitude = LONGITUDE,
          scientificname = SCIENTIFIC.NAME) %>%
  mutate (base = "eBird")



# --------------------------------------------------------------------



# create a grid to have sites
res_raster<- 0.05 # latlong degrees

# based on the extent of extracted data
grd_df <- expand.grid(x = seq(from = (extent(noronha))[1]-0.05,                                           
                              to = (extent(noronha))[2]+0.05,
                              by = res_raster),
                      y = seq(from = (extent(noronha))[3]-0.05,                                           
                              to = (extent(noronha))[4]+0.05,
                              by = res_raster))  # expand points to grid

# points object
coordinates(grd_df) <- ~x + y

# Sp points into raster
grd_raster <- (raster(grd_df,resolution = res_raster))
crs(grd_raster) <- "+proj=longlat +datum=WGS84"

# define values
values(grd_raster) <- rnorm (length(values(grd_raster)))


# load elevation map
# https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
require(elevatr)
elevation <- get_elev_raster(noronha, z = 9,expand=1)

# find site elevation
cont.elev.sp <- raster::extract(elevation, 
                                coordinates(grd_raster), 
                                fun=mean,
                                cell=T) 
# the sites
rownames(cont.elev.sp) <- seq (1,nrow(cont.elev.sp))


# all site of analysis (all raster cells)
all_sites <- rownames(cont.elev.sp)

# maps
plot(grd_raster)
plot(noronha,add=T)
plot(elevation,add=T)
plot(noronha,add=T)
points(inat_obs_df[ , c("decimallongitude", "decimallatitude")], pch = 19,col = "red")
points(gbif_data[ , c("decimallongitude", "decimallatitude")], pch = 1,cex=2)
points(vert_data_local[ , c("decimallongitude", "decimallatitude")], 
       pch = 19,col = "green")
points(ebird_data[ , c("decimallongitude", "decimallatitude")], 
       pch = 19,col = "cyan")




# --------------------------------------------------- #



overlap_sites <- lapply (list(inat_obs_df, gbif_data, vert_data_local, ebird_data), function (i) {
      
    # overlap with noronha's shape to have the sites
    # sp points
    sp_points <- data.frame (x=as.numeric(i$decimallongitude), 
                             y=as.numeric(i$decimallatitude))
    # pointss
    coordinates(sp_points) <- ~x + y
    crs (sp_points) <- crs (grd_raster) # crs
    
    # rasterize to count the number of papers per cell
    overlap_pts <- raster::extract(grd_raster,sp_points,cell=TRUE )
    
    # bind the cell identity (our site)
    sites <- overlap_pts[,"cells"]
    length(unique(overlap_pts[,"cells"])) # nsites
    
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


# species per site

data_modeling <- lapply (list(inat_obs_df, gbif_data, vert_data_local), function (i) {

  # pivot table
  detection_matrix <- dcast (sites~scientificname,
                           data = i,
                           drop=F)
    
  
  # imputing missing sites  
  to_input <-data.frame (matrix (NA, 
                                 nrow = length(all_sites[which(as.numeric(all_sites) %in% (detection_matrix$sites) !=T)]), # cells not in the table
                                ncol = ncol (detection_matrix)))
  colnames(to_input) <- colnames(detection_matrix)        
  to_input$sites <- all_sites[which(all_sites %in% detection_matrix$sites !=T)] 
  
  # bind & order
  detection_matrix_complete <- rbind (detection_matrix, to_input) # bind
  detection_matrix_complete  <- detection_matrix_complete[order(as.numeric(detection_matrix_complete$sites)),]# order
  
  # detection matrix of my species
  focus_sp_mat <- detection_matrix_complete[,which(colnames(detection_matrix_complete) %in% myspecies)]
  focus_sp_mat<-ifelse (focus_sp_mat>1,1,focus_sp_mat)
  
  # effort matrix
  effort_mat <- detection_matrix_complete[,-which(colnames(detection_matrix_complete) == "sites")]
  effort_mat [effort_mat>1] <- 1
  effort_mat <- apply (effort_mat,1,sum,na.rm=T)
  
  # list of data
  res <- list (det = focus_sp_mat,
               effort = effort_mat)
  res
  
})
  

# data modeling ebird

# pivot table for each species
# all_lists <- unique(ebird_data$SAMPLING.EVENT.IDENTIFIER)
detection_matrix <- dcast (sites~SAMPLING.EVENT.IDENTIFIER,
                           data = ebird_data[which(ebird_data$scientificname %in% myspecies
                                                   &
                                                     ebird_data$PROTOCOL.TYPE != "Historical"),],
                           drop=F,
                           value = "scientificname")

# imput sites & lists
to_input <-data.frame (matrix (NA, 
                               nrow = length(all_sites[which(as.numeric(all_sites) %in% (detection_matrix$sites) !=T)]), # cells not in the table
                               ncol = ncol (detection_matrix))) # sites
colnames(to_input) <- colnames(detection_matrix)         # set colnames
to_input$sites <- all_sites[which(all_sites %in% detection_matrix$sites !=T)]  # set site names

# bind & order
detection_matrix_complete <- rbind (detection_matrix, to_input) # bind
detection_matrix_complete  <- detection_matrix_complete[order(as.numeric(detection_matrix_complete$sites)),]# order


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


