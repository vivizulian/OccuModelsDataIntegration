# modeling

# load data
load("modeling_data.RData")


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
    BETA.ALT ~ dunif(-10, 10)
    
    
    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
  
      mu[i] <- BETA0 + BETA.ALT * altitude [i]

     # Keeping Winbugs on the track
     mu.lim[i] <- min(10, max(-10, mu[i]))  
     logit(psi[i]) <- mu.lim[i]

    }    

    ###################################
    ###### OBSERVATION MODELS #########
    ###################################

    ### GBIF 
    
    for (i in 1:nsite) {


        e1[i] <- ALPHA.GBIF * nSP.gbif [i]
        P1[i] <- 1-pow((1-0.5), e1[i])
        zP1[i] <- P1[i] * z[i]
        y.gbif[i] ~ dbern(zP1[i])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.GBIF ~ dnorm(0,0.0001)I(0,1000)# truncated in positive values
    
    
    ### VERTNET 
    
    for (i in 1:nsite) {


        e2[i] <- ALPHA.VERTNET * nSP.vertnet [i]
        P2[i] <- 1-pow((1-0.5), e2[i])
        zP2[i] <- P2[i] * z[i]
        y.vertnet[i] ~ dbern(zP2[i])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.VERTNET ~ dnorm(0,0.0001)I(0,1000)# truncated in positive values
    
    
    
    ### INATURALIST 
    
    for (i in 1:nsite) {


        e3[i] <- ALPHA.VERTNET * nSP.inat [i]
        P3[i] <- 1-pow((1-0.5), e3[i])
        zP3[i] <- P3[i] * z[i]
        y.inat[i] ~ dbern(zP3[i])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.INAT ~ dnorm(0,0.0001)I(0,1000)# truncated in positive values
    
    
    ## EBIRD
  
    for (i in 1:nsite) {
      for (j in 1:nrepEB) { # n checklists

        e4[i,j] <- ALPHA.DIST * dist [i,j] +
                   ALPHA.DURATION * duration[i,j] 
        P4[i,j] <- 1-pow((1-0.5), e4[i,j])
        zP4[i,j] <- P4[i,j] * z[i]
        yEB[i,j] ~ dbern(zP4[i,j])

      }
    }
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,1000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,2000)
    

    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    #compute the mean detection probability of each dataset: 
    muP1 <- mean(P1[])
    muP2 <- mean(P2[])
    muP3 <- mean(P3[])
    muP4 <- mean(P4[,])
    
    
    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    }## end of the model
    ",fill = TRUE)
sink()

# bound data
str(win.data <- list(nsite= length(modeling_data$site_covs[,1]),

                     y.gbif = modeling_data$gbif_detection,
                     nSP.gbif = modeling_data$gbif_effort, 
                     
                     y.vertnet = modeling_data$vertnet_detection,
                     nSP.vertnet = modeling_data$vertnet_effort,
                     
                     y.inat = modeling_data$inat_detection,
                     nSP.inat = modeling_data$inat_effort, 
                     
                     yEB = as.matrix(modeling_data$ebird_detection[,-1]),
                     dist = as.matrix(modeling_data$ebird_duration[,-1]),
                     duration = as.matrix(modeling_data$ebird_distance[,-1]),
                     nrepEB = ncol(modeling_data$ebird_detection[,-1]),
                     
                     altitude=as.numeric(scale (modeling_data$site_covs[,2]))))

# initial values
set.seed(32)
inits <- function() {list(z = rep (1, nrow(modeling_data$ebird_detection)),
                          ALPHA.DIST = rnorm (1,mean=2),
                          ALPHA.DURATION = rnorm (1,mean=20))}


params <- c(#"BETA0", 
  #"BETA.CAMPO","BETA.AGRI",
  #"ALPHA.GBIF","ALPHA.DIST", "ALPHA.DURATION","ALPHA.OBSERVER",
  #"ALPHA.PICT", "ALPHA.SONG", 
  #"muP1","muP2","muP5",
  "z","psi","fs.z"
)

# settings
ni<-50; nc<-3; na<-5; nb <- 25; nt<-1
require(R2WinBUGS)
model_output <- bugs(data = win.data, 
               parameters.to.save = params, 
               model.file = "model_integration.txt", 
               inits = inits, 
               n.chains = nc, 
               n.thin = nt, 
               n.iter = ni, 
               n.burnin = nb,
               codaPkg=F, 
               DIC=TRUE, 
               debug=T,
               bugs.directory="C:/Program Files/WinBUGS14/", 
               program= "WinBUGS"
)

# save it
save (model4,file=here("output", "model4_bugs.RData"))
