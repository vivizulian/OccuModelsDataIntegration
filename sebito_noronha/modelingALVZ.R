# modeling
rm(list=ls())
require(here)
require(raster)
# load data
load(here("sebito_noronha", "modeling_dataVZ.RData"))


# ------------------------------------------------------
# modelo 4
# com integracao e deteccao

# write model
sink("model_integration.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################

    # PRIORS
    BETA0 ~ dunif(-10,10)
    BETA.ALT ~ dunif(-10,10)
    
    
    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
  
      mu[i] <- BETA0 + BETA.ALT * altitude[i]

      # Keeping Winbugs on the track
      mu.lim[i] <- min(10, max(-10, mu[i]))  
      logit(psi[i]) <- mu.lim[i]

    }    

    ###################################
    ###### OBSERVATION MODELS #########
    ###################################

    ### GBIF 
    
    for (l in 1:nObsGB) {


        e1[l] <- ALPHA.GBIF * nSP.gbif[l]
        P1[l] <- 1-pow((1-0.5), e1[l])
        zP1[l] <- P1[l] * z[siteGB[l]]
        y.gbif[l] ~ dbern(zP1[l])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.GBIF ~ dnorm(0,0.0001)I(0,10000)# truncated in positive values
    
    
    ### INATURALIST 
    
    for (r in 1:nObsIN) {


        e3[r] <- ALPHA.INAT * nSP.inat[r]
        P3[r] <- 1-pow((1-0.5), e3[r])
        zP3[r] <- P3[r] * z[siteIN[r]]
        y.inat[r] ~ dbern(zP3[r])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.INAT ~ dnorm(0,0.0001)I(0,10000)# truncated in positive values
    
    
    ## EBIRD
  
    for (n in 1:nObsEB) {

        e4[n] <- ALPHA.DIST * dist [n] + ALPHA.DURATION * duration[n] + ALPHA.SP * nSP.eb[n]
        P4[n] <- 1-pow((1-0.5), e4[n])
        zP4[n] <- P4[n] * z[siteEB[n]]
        yEB[n] ~ dbern(zP4[n])

    }
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,10000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,10000)
    ALPHA.SP ~ dnorm(0,0.0001)I(0,10000)
  

    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    #compute the mean detection probability of each dataset: 
    muP1 <- mean(P1[])
    #muP2 <- mean(P2[])
    muP3 <- mean(P3[])
    muP4 <- mean(P4[])
    
    
    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    }## end of the model
    ",fill = TRUE)
sink()

# bound data
str(data <- list(nsite = ncell(modeling_data$grid),
                     y.gbif = modeling_data$gbif_detection,
                     nSP.gbif = modeling_data$gbif_effort, 
                     nObsGB = length(modeling_data$gbif_detection),
                     siteGB = modeling_data$GBsite,
                     y.inat = modeling_data$inat_detection,
                     nSP.inat = modeling_data$inat_effort,
                     nObsIN = length(modeling_data$inat_detection),
                     siteIN = modeling_data$INsite,
                     yEB = modeling_data$ebird_detection,
                     dist = as.numeric(modeling_data$ebird_distance),
                     duration = as.numeric(modeling_data$ebird_duration),
                     nSP.eb = modeling_data$ebird_nsp,
                     nObsEB = length(modeling_data$ebird_detection),
                     siteEB = modeling_data$EBsite,
                     altitude = as.numeric(scale((modeling_data$site_covs[,2])))))

# initial values
#set.seed(32)
inits <- function() {list(z = rep(1, data$nsite))}


params <- c("BETA0", 
  "BETA.ALT", "ALPHA.GBIF",
  "ALPHA.DIST", "ALPHA.DURATION","ALPHA.SP",
  "ALPHA.INAT",
  "muP1","muP3","muP4",
  "z","psi","fs.z")

# settings
ni<-5000; nc<-3; na<-500; nb <- 2500; nt<-5
require(jagsUI)
require(R2WinBUGS)
model_output <- bugs(data = data, 
                     parameters.to.save = params, 
                     model.file = "model_integration.txt", 
                     inits = inits, 
                     n.chains = nc, 
                     n.thin = nt, 
                     n.iter = ni, 
                     n.burnin = nb,
                     codaPkg = F, 
                     DIC = TRUE, 
                     debug = T,
                     bugs.directory = "C:/Program Files/WinBUGS14/", 
                     program = "WinBUGS")



# save it
save (model_output, file=here("output", "model_output.RData"))



### plot the map ####

require(rgdal)
require(here)
require(raster)
require(ggplot2)
library(rgeos)
library(maptools)
library(ggsn)

noronha <- readOGR("FernandoNoronha.shp")
hex = modeling_data$grid

z <- model_output$mean$z
#z <- data.frame(seq=seq(1:3701), out2$mean$psi, out2$sd$psi, out2$mean$z, out2$sd$z)
rownames(z) <- hex@plotOrder

f.mun <- fortify(hex)
#f.mun$id<-as.numeric(as.factor(f.mun$id))

intervalos <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
cortes <- cut(z, intervalos, include.lowest=TRUE)

levels(cortes)[levels(cortes)== "[0,0.2]"]  <- "[0 - 0.2["
levels(cortes)[levels(cortes)== "(0.2,0.4]"]  <- "[0.2 - 0.4["
levels(cortes)[levels(cortes)== "(0.4,0.6]"]  <- "[0.4 - 0.6["
levels(cortes)[levels(cortes)== "(0.6,0.8]"]  <- "[0.6 - 0.8["
levels(cortes)[levels(cortes)== "(0.8,1]"]  <- "[0.8 - 1.0]"

z<-cbind(z, cores=cortes)

f.mun$cores <- z[match(f.mun$id, as.numeric(rownames(z))),2]

(p <- ggplot() + geom_polygon(data=noronha, aes(x=long, y=lat, group = group),colour="gray70",fill="gray90" ,alpha = 0.5))

(b <- p +  geom_polygon(data=f.mun, aes(x=long, y=lat, group=group, fill=factor(cores)), colour = "black", size=0.5) + 
  scale_fill_brewer(palette="YlOrRd", name="Mean Z", labels = c("[0 - 0.2[", "[0.2 - 0.4[", "[0.4 - 0.6[",
                                                                "[0.6 - 0.8[", "[0.8 - 1.0]")))
    