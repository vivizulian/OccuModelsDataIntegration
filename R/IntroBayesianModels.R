
##### Introduction to Bayesian Analysis #####
## Code from Kery and Royle 2016, chapters 5 and 10 ##

#In this script you will find:
#Presence-absence data simulation without covariates (1A) and with covariates (1D)
#Occupancy models analyzed using unmarked package (1B and 1E)
#Occupancy models analyzed using jags, in a Bayesian framework (1C and 1F)
#Count data simulation (2A and 2D)
#Abundance models analyzed using unmarked package (2B and 2E)
#Abundance models analyzed using jags, in a Bayesian framework (2C and 2F)

#The same data set is analyzed using unmarked package and jagsUI package, in a bayesian framework.

#If you do not have the jagsUI installed in your R, run:
install.packages("jagsUI")


####___________1) Occupancy models (Chapter 10)_______________####


#### 1A) Presence-absence data simulation ####

# Choose sample sizes and prepare observed data array y
set.seed(24) # So we all get same data set
M <- 100 # Number of sites
J <- 2 # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8 # Probability of occupancy or presence
p <- 0.5 # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}


#### 1B) Run using unmarked package ####

library(unmarked)
umf <- unmarkedFrameOccu(y = y) # Create unmarked data frame
summary(umf) # Summarize data frame
(fm1 <- occu(~1 ~1, data = umf)) # Fit model

# Get estimates on probability scale
backTransform(fm1, "state") 
backTransform(fm1, "det")


#### 1C) Run using jagsUI package ####

# Bundle data and summarize data bundle
str(data <- list(y = y, M = nrow(y), J = ncol(y)) )

# Specify model in BUGS language
sink("model1.txt")
cat("
model {
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) { # Loop over sites
    z[i] ~ dbern(psi) # State model
    for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p) # Observation model (only JAGS !)
      # y[i,j] ~ dbern(mu[i]) # For WinBUGS define 'straw man'
    }
    # mu[i] <- z[i]*p # Only WinBUGS
  }
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, 1, max) # Avoid data/model/inits conflict
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "p")

# MCMC settings
ni <- 50 ; nt <- 2 ; nb <- 10 ; nc <- 3; na=2

# Call JAGS and summarize posteriors
library(jagsUI)
out1 <- jags(data, inits, params, "model1.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb)
print(out1, dig = 3)

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out1, params = c("psi", "p"))


#### 1D) Simulate more complex presence-absence data with covariates ####

# Choose sample sizes and prepare obs. data array y
set.seed(1) # So we all get same data set
M <- 100 # Number of sites
J <- 3 # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # Sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0 # Logit-scale intercept
beta1 <- 3 # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability
# plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship

# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi) # True presence/absence

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2 # Logit-scale intercept
alpha1 <- -3 # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
# plot(p ~ wind, ylim = c(0,1)) # Look at relationship

# Take J [ 3 presence/absence measurements at each site
for(j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j])
}


#### 1E) Run using package unmarked ####

# Load unmarked, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(y = y, # Pres/Abs measurements
                         siteCovs = data.frame(vegHt = vegHt), # site-specific covs.
                         obsCovs = list(wind = wind)) # obs-specific covs.
summary(umf)

# Fit model and extract estimates
# Detection covariates follow first tilde, then occupancy covariates
summary(fm2 <- occu(~wind ~vegHt, data=umf))


#### 1F) Ran using the package jagsUI ####

# Bundle and summarize data set
str(data <- list(y = y, vegHt = vegHt, wind = wind, M = nrow(y), J = ncol(y)))
# Specify model in BUGS language
sink("model2.txt")
cat("
model {
  # Priors
  mean.p ~ dunif(0, 1) # Detection intercept on prob. scale
  alpha0 <- logit(mean.p) # Detection intercept
  alpha1 ~ dunif(-20, 20) # Detection slope on wind
  mean.psi ~ dunif(0, 1) # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi) # Occupancy intercept
  beta1 ~ dunif(-20, 20) # Occupancy slope on vegHt
  
  # Likelihood
  # True state model for the partially observed true state
  for (i in 1:M) {
    z[i] ~ dbern(psi[i]) # True occupancy z at site i
    logit(psi[i]) <- beta0 + beta1 * vegHt[i]
    
    # Observation model for the actual observations
    for (j in 1:J) {
      y[i,j] ~ dbern(p.eff[i,j]) # Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j] # 'straw man' for WinBUGS
      logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
    }
  }
  # Derived quantities
  N.occ <- sum(z[]) # Number of occupied sites among sample of M
  psi.fs <- N.occ/M # Proportion of occupied sites among sample of M
  
}
",fill = TRUE)
sink()

# Initial values: must give for same quantities as priors given !
zst <- apply(y, 1, max) # Avoid data/model/inits conflict
inits <- function(){list(z = zst, mean.p = runif(1), alpha1 = runif(1), mean.psi
                         = runif(1), beta1 = runif(1))}

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs", "p", "z", 
            "psi")

# MCMC settings
ni <- 2500 ; nt <- 10 ; nb <- 2000 ; nc <- 3
#ni <- 25000 ; nt <- 10 ; nb <- 2000 ; nc <- 3

# Call WinBUGS from R (ART 2 min) and summarize posteriors
out2 <- jags(data, inits, params, "model2.txt", n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)
print(out2, dig = 3)

## Plot the chains to check convergence:
library(mcmcOutput)
mcmcOutput::diagPlot(out2, params = c("alpha1","beta1", "beta0", "alpha0"))



####_____________2) Abundance models___________________####


#### 2A) Simulate the data (Chapter 5) ####

# Choose sample sizes and prepare observed data array C
set.seed(24) # So we all get same data set
M <- 150 # Number of sites
J <- 2 # Number of abu. measurements per site (rep. counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
lambda <- 2.5 # Expected abundance
p <- 0.4 # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}


#### 2B) Run using unmarked package ####

library(unmarked) 
# Create the data frame
umf <- unmarkedFramePCount(y = C) 
summary(umf) # Summarize

# Fit model: get estimates on link scale
(fm3 <- pcount(~1 ~1, data = umf)) 

# Get estimates on natural scale
backTransform(fm3, "state") 
backTransform(fm3, "det")


#### 2C) Run using JAGS ####

# Bundle and summarize data set
win.data <- list(C = C, M = nrow(C), J = ncol(C))
str(win.data) # Look at data

# Specify model in BUGS language
sink("model3.txt")
cat("
model {
  # Priors
  lambda ~ dgamma(0.001, 0.001)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {
    N[i] ~ dpois(lambda) # State model
    for (j in 1:J) {
      C[i,j] ~ dbin(p, N[i]) # Observation model
    }
  }
}
",fill = TRUE)
sink()

# Initial values
Nst <- apply(C, 1, max) # Avoid data/model/inits conflict
inits <- function(){list(N = Nst)}

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 2500 ; nt <- 10 ; nb <- 2000 ; nc <- 3  #no convergence
#ni <- 25000 ; nt <- 20 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 1 min) and summarize posteriors
library(jagsUI)
out3 <- jags(win.data, inits, params, "model3.txt", n.chains = nc, n.thin = nt, n.iter = ni,
            n.burnin = nb)
print(out3, dig = 3)

## Plot the chains to check convergence:
library(mcmcOutput)
mcmcOutput::diagPlot(out3, params = c("lambda", "p"))


#### 2D) Simulate more complex count data with covariates ####

# Choose sample sizes and prepare observed data array y
set.seed(1) # So we all get same data set
M <- 100 # Number of sites
J <- 3 # Number of repeated abundance measurements
C <- matrix(NA, nrow = M, ncol = J) # to contain the observed data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for abundance model and compute lambda
beta0 <- 0 # Log-scale intercept
beta1 <- 2 # Log-scale slope for vegHt
lambda <- exp(beta0 + beta1 * vegHt) # Expected abundance
plot(vegHt, lambda, type = "l", lwd = 3) # Expected abundance

# Draw local abundance and look at data so far
N <- rpois(M, lambda)
points(vegHt, N) # Add realized abundance to plot

# Create a covariate called wind - detection covariate
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2 # Logit-scale intercept
alpha1 <- -3 # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
plot(p ~ wind, ylim = c(0,1)) # Look at relationship

# Take J [ 3 abundance measurements at each site
for(j in 1:J) {
  C[,j] <- rbinom(M, N, p[,j])
}

# Expected (lambda) and realized abundance (N) and measurements (C)
cbind(lambda=round(lambda,2), N=N, C1=C[,1], C2=C[,2], C3=C[,3])


#### 2E) Run using unmarked ####

# Load unmarked, format data in unmarked data frame and summarize
library(unmarked)
umf <- unmarkedFramePCount(y = C, # Counts matrix
                           siteCovs = data.frame(vegHt = vegHt), # Site covariates
                           obsCovs = list(wind = wind)) # Observation covs
summary(umf)

# Fit model and extract estimates
# linear model for p follows first tilde, then comes linear model for lambda
summary(Nmix1 <- pcount(~wind ~vegHt, data=umf, control=list(trace=T, REPORT=1)))


#### 2F) Run using JAGS ####

# Bundle data
data <- list(C = C, M = nrow(C), J = ncol(C), wind = wind, vegHt = vegHt)
str(data)

# Specify model in BUGS language
cat(file = "model4.txt", "
model {
  # Priors
  beta0 ~ dunif(-10, 10) # Abundance intercepts
  beta1 ~ dunif(-10, 10) # Abundance slopes
  alpha0 ~ dunif(-10, 10) # Detection intercepts
  alpha1 ~ dunif(-10, 10) # Detection slopes
  
  # Likelihood
  # Ecological model for true abundance
  for (i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1 * vegHt[i]
    
      # Some intermediate derived quantities
      critical[i] <- step(2-N[i]) # yields 1 whenever N is 2 or less
      z[i] <- step(N[i]-0.5) # Indicator for occupied site
    
      # Observation model for replicated counts
      for (j in 1:J){
        C[i,j] ~ dbin(p[i,j], N[i])
        logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
      }
  }
  
  # Derived quantities: functions of latent variables and predictions
  Nocc <- sum(z[]) # Number of occupied sites among sample of M
  Ntotal <- sum(N[]) # Total population size at M sites combined
  N.critical <- sum(critical[]) # Number of populations with critical size
  
}
")

# Initial values
Nst <- apply(C, 1, max)+1 # Important to give good inits for latent N
inits <- function() list(N = Nst)
# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "Nocc", "Ntotal", 
            "N.critical", "lambda", "p") 

# MCMC settings
nc <- 3 ; ni <- 22000 ; nb <- 2000 ; nt <- 10

# Call JAGS and summarize posteriors
library(jagsUI)
out4 <- jags(data, inits, params, "model4.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

traceplot(out4, param = c('alpha0', 'alpha1', 'beta0', 'beta1', 'Nocc', 'Ntotal',
                         'Nhab', 'N.critical'))
print(out4, 2)





