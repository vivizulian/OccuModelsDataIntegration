
#####################################################################
### Workshop: Mapeamento probabilistico da distribuicao  de especies 
###           baseado na integracao de dados de ciencia cidada

### Ministrantes: Viviane Zulian e Andre Luza ###
#####################################################################

### Codigo em R e BUGS para pratica ###

# -------------------------------------------------
#           Modelos de distribuicao de especies
#     Lista de implementacoes
#     1 - modelos sem integracao, e sem consideracao de esforco amostral (eBird) - GLM
#     2 - modelos sem integracao, e consideracao de esforco amostral (eBird) - HM 'unmarked'
#     2 - modelos sem integracao, e consideracao de esforco amostral (eBird) - HM BUGS
#     3 - modelos integracao (agregacao), e sem esforco amostral (eBird) - GLM
#     3 - modelos integracao (agregacao), e sem esforco amostral (eBird) - BUGS
#     4 - modelos com integracao, e consideracao de esforco amostral (eBird, GBIF, WikiAves) - BUGS

# carregar os pacotes
source ('R/packagesR.R')
source ('R/functions.R')

# carregar dados
load(here("data","organized_data", "input_GBIF.RData")) ## gbif
load(here("data","organized_data", "input_ebird.RData")) ## eBird
load(here("data","organized_data", "input_wikiaves.RData")) # WikiAves



#-------------------------
# carregar dados espaciais
#--------------------------



# carregar shapefile RS
# carregar em uma projecao Lambert, para trabalhar adequadamente com areas
load(here ("data","shape_munRS","shapeRS_lambert.RData"))

## calcular municipality_area, km
area_mun <- gArea (shape_RS,byid=T)/1000.00

#----------------------------
# covariaveis de sitio (municipio)
#----------------------------

load(here("data","organized_data", "dados_covariaveis.RData"))

# campo
campo <- usos_tabela$`Campo seco.area` + 
                      usos_tabela$`Campo umido.area` +
                      usos_tabela$`Campo de feixe de restinga.area`+
                      usos_tabela$`Campo em regeneracao.area`
# agricultura
agricultura_uso <- usos_tabela$`Agricultura de sequeiro.area`
# correlacao
cor(campo,agricultura_uso)


# ---------------------------------------------------------- #

# padronizacoes
## area de habitat campestre e agricultura por municipio
## multiplicar a cobertura de habitat pela area do municipio
area_campo <- area_mun*campo # para plot
habitat_campo <- area_mun*campo # para padronizar e ir no modelo 
area_lavoura <- area_mun*agricultura_uso # para plot
habitat_lavoura <- area_mun*agricultura_uso # para padronizar e ir no modelo

## raiz quadrada, e entao padronizar pela media e sd
habitat_campo <- decostand (sqrt(habitat_campo),"standardize")
habitat_lavoura <- decostand (sqrt(habitat_lavoura),"standardize")


# CRIAR DIRETORIO PARA HOSPEDAR RESULTADos
dir.create("output")

# salvar covariaveis
save(habitat_campo, area_campo,
     habitat_lavoura,area_lavoura,
     file=here("output", "covariaveis.RData"))


# colar area de campo e agricultura no shape
shape_RS@data$campo <- area_campo # campo
shape_RS@data$agricultura <- area_lavoura # agricultura

# mapa campo
mapa_campo <- funcao_mapa_variaveis(shape= shape_RS,variavel = shape_RS$campo,
                      titulo = "Area de Campos Naturais", 
                      cor.min = "white", 
                      cor.max = "darkgreen",
                      cor.na = "gray50")

# mapa agricultura
mapa_agricultura <- funcao_mapa_variaveis(shape= shape_RS,variavel = shape_RS$agricultura,
                      titulo = "Area de Agricultura", 
                      cor.min = "white", 
                      cor.max = "darkred",
                      cor.na = "gray50")



# bind maps and save
pdf (here ("output", "maps_land_use.pdf"),heigh=7,width=5)
grid.arrange(mapa_campo,
             mapa_agricultura,
             ncol=2)

dev.off()

# -------------------------------------------------------------------------

# MAPAS DAS OBSERVACOES


# carregar mapa do RS (em latlong)
shape_RS <- readOGR(dsn=here("data","shape_munRS"), 
                    layer="43MUE250GC_SIR",
                    encoding = "UTF-8",use_iconv = T)


## obter os lagos para pinta-los com cores diferentes depois
lagos <- shape_RS [c(96,250),]

## remover os lagos
shape_RS <- shape_RS [-c(96,250),]

## shape south america
southAme<- readOGR(dsn=here ("data","South_America"),encoding="latin1", layer="South_America")
BR_AR_URU<- southAme [southAme@data$COUNTRY == "Paraguay" | southAme@data$COUNTRY == "Brazil" | southAme@data$COUNTRY == "Argentina" | southAme@data$COUNTRY == "Uruguay", ]
crs(BR_AR_URU)<-crs(shape_RS)


# mapa de base (south america and lakes)

## colocar o shape da america do sul = comum a todos os mapas
a <- ggplot() + geom_polygon (data=BR_AR_URU, aes(x=long, 
                                                  y=lat, 
                                                  group=group),
                              size = 0.1, fill="gray90", 
                              colour="gray75",alpha=1) +
  coord_fixed (xlim = c(-57.5, -49),  
               ylim = c(-34, -27), 
               ratio = 1) 

## inserir os lagos = comum a todos os mapas
b <- a + geom_polygon (data=lagos, aes(x=long, 
                                       y=lat, 
                                       group=group), 
                       fill="lightcyan",
                       colour = "lightcyan", 
                       size=1)


## abrir deteccoes do eBird
load (here("data","organized_data",  "input_ebird.RData"))## gbif

# observacoes
ebird_rhea <- apply (y.ebird, 1, max,na.rm=T) # 
map_ebird <- funcao_mapa_observacoes (shape = shape_RS, observacoes = ebird_rhea, titulo="eBird")

## abrir deteccoes GIBF
load (here("data","organized_data", "input_GBIF.RData"))## gbif
map_gbif <- funcao_mapa_observacoes (shape = shape_RS, 
                         observacoes = dados_det_ema_gbif$det, 
                         titulo="GBIF")


## abrir deteccoes wikiaves
load (here("data","organized_data", "input_wikiaves.RData"))
map_wikiaves <- funcao_mapa_observacoes (shape = shape_RS, 
                         observacoes = dados_wikiaves$RHAMERICANA, 
                         titulo="WikiAves")

# dados agregados
cores_agregado <- data.frame(eBird=ebird_rhea, 
                             GBIF=dados_det_ema_gbif$det, 
                             WikiAves=dados_wikiaves$RHAMERICANA)
# agregar
agg_data <- apply(cores_agregado, 1, max,na.rm=T)## municipios com deteccoes
apply(cores_agregado, 2, function (i) sum (i > 0,na.rm=T))## numero de deteccoes por base
table(apply (cores_agregado,1,max,na.rm=T)>0) # numero de municipios com deteccao e sem deteccao

# mapa
map_agg <- funcao_mapa_observacoes (shape = shape_RS, 
                         observacoes = agg_data, 
                         titulo="Agregado")

# uma unica legenda
legenda_comum_data <- get_legend (map_agg)

######### arranjar o painel
pdf(file=here ("output","maps_observations.pdf"),width = 7,height = 5,family="serif")

grid.arrange(legenda_comum_data,
             map_agg + theme(legend.position = "none"), 
             map_ebird+ theme(legend.position = "none"),
             map_gbif+ theme(legend.position = "none"),
             map_wikiaves+ theme(legend.position = "none"),
             ncol=5,nrow=6,
             layout_matrix = rbind(c(1,1,1,1,1), 
                                   c(2,2,2,2,2),
                                   c(2,2,2,2,2),
                                   c(2,2,2,2,2),
                                   c(NA,3,4,5,NA),
                                   c(NA,3,4,5,NA))) 
dev.off()


#################################################################################
# ----------------------------------------------------------------------------- #

# MODEL 1
#     1 - modelos sem integracao, e sem consideracao de esforco amostral (eBird) - GLM
# ------ relacionando a 'presenca' da ema com covariaveis de ambiente (area de campo e agri, por municipio)
# ------ dados do eBird

# agregar dados das listas do eBird
agg_ebird <- apply (y.ebird,1,max,na.rm=T)
agg_ebird [is.infinite(agg_ebird)] <-NA # se nao teve amostra, colocar NA
# colar covariaveis
agg_ebird<-data.frame(agg_ebird,
                    campo=habitat_campo,
                   agri=habitat_lavoura)

# modelo para municipios com dados
agg_ebird_subset <- agg_ebird[which(is.na(agg_ebird$agg_ebird)!=T),]
model1<- glm (agg_ebird ~ campo + agri,
              data=agg_ebird_subset,
              family=binomial(link = "logit"))
# proporcao de municipios ocupados
plogis (model1$coefficients[1])

# save it
save (model1,file=here("output", "out_model1_glm.RData"))

# --------------------------------------------------------------------
# MODEL 2
#     2 - modelos sem integracao, e consideracao de esforco amostral (eBird) - HM 'unmarked'

# organize os dados para adequar-se ao framework
umf <- unmarkedFrameOccu(y= y.ebird,
                         siteCovs=  cbind(data.frame (agri = habitat_lavoura, 
                                                      campo = habitat_campo)), 
                         obsCovs=
                           list (duracao = dist.duration,
                                 distancia = dist.ebird,
                                 observadores=dist.observers)								    
)
summary(umf)

# modelo
# nulo
null <- occu(~duracao + distancia + observadores # deteccao
               ~ 1 , # ocupacao 
               data=umf)
# alternativa 1
model2 <- occu(~duracao + distancia + observadores # deteccao
               ~agri + campo , # ocupacao 
               data=umf)
# alternativa 2
model2b <- occu(~duracao + distancia + observadores # deteccao
               ~ campo , # ocupacao 
               data=umf)
models <- list(null=null,m1=model2, m2=model2b) # lista de modelos
fmList <- fitList(fits = models) # resultados
modSel(fmList, nullmod="null") # selecao de modelos

# save it
save (model2,
      file=here("output", "out_model2_unmarked.RData"))

# proporcao de municipios ocupados
plogis(coef(model2,"state")[1])

# ------------------------------------------------------ #
# modelo 2 em linguagem bugs


# GLOBAL MCMC settings
## short form
# na <- 30; nb <- 40; ni <- 50; nc <- 3; nt <- 1
na <- 3000; nb <- 4000; ni <- 5000; nc <- 3; nt <- 1



# escrever modelo
sink("model2_NoIntegration_Detection.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################
    
    BETA0 ~ dunif(-10,10)
    BETA.CAMPO ~ dunif(-10,10)
    BETA.AGRI ~ dunif(-10,10)
    
    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
      mu[i] <- BETA0 + BETA.CAMPO * campo[i] + # occ model
                      BETA.AGRI * agri [i] 
                      
      mu.lim[i] <- min(10, max(-10, mu[i]))  
      logit(psi[i]) <- mu.lim[i]

    }    

    ###################################
    ###### OBSERVATION MODEL #########
    ###################################

   
    for (i in 1:nsite) {
      for (j in 1:nrepEB) {

        
        yEB [i,j] ~ dbern(muY[i,j])
        muY[i,j] <- p[i,j] * z [i]
        logit(p[i,j]) <-  ALPHA + ALPHA.DIST * dist [i,j] + 
                                  ALPHA.DURATION * duration [i,j] + 
                                  ALPHA.OBSERVER * observer [i,j]

      }
    }
    
    # priors 
    ALPHA ~ dnorm (0,0.001)
    ALPHA.DIST ~ dnorm (0,0.001)
    ALPHA.DURATION ~ dnorm (0,0.001)
    ALPHA.OBSERVER ~ dnorm (0,0.001)
    
    # ================================= 
    # fit statistics
    # Computation of fit statistic (for Bayesian p-value)
    # assess modelo fit using chi-squared discrepancy
    # compute fit statistics for observed data
    
     for (i in 1:nsite) {
      for (j in 1:nrepEB) {
        eval[i,j]<-p[i,j] * z [i]
        E[i,j] <- pow((yEB [i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
        y.new[i,j]~dbern(z[i]) 
        E.new[i,j] <- pow((y.new [i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
        
      }
     }
    
    # aggregated fit statistics
    fit<-sum(E[,])# Discrepancy for actual data set
    fit.new<-sum(E.new[,]) # Discrepancy for replicate data set
    test<-step(fit.new - fit) # Test whether new data set more extreme
    bpvalue<-mean(test) # Bayesian p-value
   
    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    } ## end of the model
    ",fill = TRUE)
sink()

# define a subset of ebird data -- too much missing lists
n_so_ebird <- round (mean(unlist(
  lapply (seq(1,nrow(y.ebird)), function (i) 
    sum(is.na(y.ebird[i,])!=T)))
)
)

# histograma
hist(unlist(
  lapply (seq(1,nrow(y.ebird)), function (i) 
    sum(is.na(y.ebird[i,])!=T))),
  xlab="Number of checklists",main="")
abline(v=n_so_ebird,lwd=2)

## initial values
set.seed(32)
zint <- apply (y.ebird,1,max,na.rm=T)
zint[is.infinite(zint)]<- 0
inits = function() {list(z = zint,
                         ALPHA.DIST = rnorm (1,mean=2),
                         ALPHA.DURATION = rnorm (1,mean=20),
                         ALPHA.OBSERVER = rnorm (1,mean=2))}

# bound data
str(winbugs.data <- list(nsite= dim(y.ebird)[1],
                      yEB = y.ebird [,1:n_so_ebird],
                      dist = dist.ebird [,1:n_so_ebird],
                      duration = dist.duration [,1:n_so_ebird],
                      observer = dist.observers [,1:n_so_ebird],
                      nrepEB = ncol(dist.observers [,1:n_so_ebird]),
                      campo=habitat_campo[,1],
                      agri=habitat_lavoura[,1]
))

# parametros para monitorar nas MCMC
params <- c("ALPHA.DIST", "ALPHA.DURATION","ALPHA.OBSERVER",
            "BETA0","BETA.CAMPO","BETA.AGRI",
            "z","psi","fs.z","p",
            "bpvalue", "fit", "fit.new"
    )


# run model
out_model2 <- bugs(data = winbugs.data, 
                  parameters.to.save = params, 
                  model.file = "model2_NoIntegration_Detection.txt", 
                  inits = inits, 
                  n.chains = nc, 
                  n.thin = nt, 
                  n.iter = ni, 
                  n.burnin = nb,
                  debug=F,
                  codaPkg=F, 
                  DIC=TRUE, 
                  bugs.directory="C:/Program Files/WinBUGS14/", 
                  program= "WinBUGS")

# save it
save (out_model2,file=here("output", "out_model2_bugs.RData"))

out_model2$mean$fs.z
out_model2$mean$bpvalue

# high density interval
hpd_vals <- HPDinterval(as.mcmc(out_model2$sims.list$ALPHA.DIST), prob=0.95)

ggplot () + 
  geom_density(data = data.frame (y=out_model2$sims.list$ALPHA.DIST),
               aes (x=y),fill="red",alpha =0.3)+
  geom_segment(aes(x = hpd_vals[1], 
                   y = 0, 
                   xend = hpd_vals[2], 
                   yend = 0),
               size=5,col="#5A8F7B") + 
  geom_segment(aes(x = out_model2$mean$ALPHA.DIST, 
                   y = 0, 
                   xend = out_model2$mean$ALPHA.DIST, 
                   yend = 0.018),
               size=3,col="#5A8F7B") 
  


# =======================================================

# modelo 2B, com modelo de deteccao alternativo (parecido com IDM)


# escrever modelo
sink("model2B_NoIntegration_Detection.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################
    
    BETA0 ~ dunif(-10,10)
    BETA.CAMPO ~ dunif(-10,10)
    BETA.AGRI ~ dunif(-10,10)
    
    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
      mu[i] <- BETA0 + BETA.CAMPO * campo[i] + # occ model
                      BETA.AGRI * agri [i] 
                      
      mu.lim[i] <- min(10, max(-10, mu[i]))  
      logit(psi[i]) <- mu.lim[i]

    }    

    ###################################
    ###### OBSERVATION MODEL #########
    ###################################

   
    for (i in 1:nsite) {
      for (j in 1:nrepEB) {

        e2[i,j] <- ALPHA.DIST * dist [i,j] +
                   ALPHA.DURATION * duration [i,j] +
                   ALPHA.OBSERVER * observer [i,j]
 
        P2[i,j] <- 1-pow((1-0.5), e2[i,j])
        zP2[i,j] <- P2[i,j] * z [i]
        yEB [i,j] ~ dbern(zP2[i,j])
        
      

      }
    }
    
    # priors 
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,1000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,2000)
    ALPHA.OBSERVER  ~ dnorm(0,0.0001)I(0,1000)


    # ================================= 
    # fit statistics
    # Computation of fit statistic (for Bayesian p-value)
    # assess modelo fit using chi-squared discrepancy
    # compute fit statistics for observed data
    
     for (i in 1:nsite) {
      for (j in 1:nrepEB) {
        eval[i,j]<-P2[i,j] * z [i]
        E[i,j] <- pow((yEB [i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
        y.new[i,j]~dbern (zP2[i,j]) 
        E.new[i,j] <- pow((y.new [i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
        
      }
     }
    
    # aggregated fit statistics
    fit<-sum(E[,])# Discrepancy for actual data set
    fit.new<-sum(E.new[,]) # Discrepancy for replicate data set
    test<-step(fit.new - fit) # Test whether new data set more extreme
    bpvalue<-mean(test) # Bayesian p-value
   
    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    #compute the mean detection probability of each dataset: 
    muP <- mean(P2[,])

    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    } ## end of the model
    ",fill = TRUE)
sink()


## initial values
set.seed(32)
zint <- apply (y.ebird,1,max,na.rm=T)
zint[is.infinite(zint)]<- 0
inits = function() {list(z = zint,
                         ALPHA.DIST = rnorm (1,mean=2),
                         ALPHA.DURATION = rnorm (1,mean=20),
                         ALPHA.OBSERVER = rnorm (1,mean=2))}

# bound data
str(winbugs.data <- list(nsite= dim(y.ebird)[1],
                         yEB = y.ebird [,1:n_so_ebird],
                         dist = dist.ebird [,1:n_so_ebird],
                         duration = dist.duration [,1:n_so_ebird],
                         observer = dist.observers [,1:n_so_ebird],
                         nrepEB = ncol(dist.observers [,1:n_so_ebird]),
                         campo=habitat_campo[,1],
                         agri=habitat_lavoura[,1]
))

# parametros para monitorar nas MCMC
params <- c("ALPHA.DIST", "ALPHA.DURATION","ALPHA.OBSERVER",
            "BETA0","BETA.CAMPO","BETA.AGRI",
            "muP",
            "z","psi","fs.z","P2",
            "bpvalue", "fit", "fit.new"
)

# run model
out_model2B <- bugs(data = winbugs.data, 
                   parameters.to.save = params, 
                   model.file = "model2B_NoIntegration_Detection.txt", 
                   inits = inits, 
                   n.chains = nc, 
                   n.thin = nt, 
                   n.iter = ni, 
                   n.burnin = nb,
                   debug=T,
                   codaPkg=F, 
                   DIC=TRUE, 
                   bugs.directory="C:/Program Files/WinBUGS14/", 
                   program= "WinBUGS")

# save it
save (out_model2B,file=here("output", "out_model2B_bugs.RData"))

out_model2B$mean$fs.z
out_model2B$mean$bpvalue


# -----------------------------------------------------
# modelo 3
# com integracao (agregacao de dados), sem deteccao

dados_agregados <- cbind (WikiAves=dados_wikiaves$RHAMERICANA,
                          eBird = apply (y.ebird,1,max,na.rm=T),
                           GBIF= dados_det_ema_gbif$det)
dados_agregados<- apply(dados_agregados,1,max,na.rm=T)
dados_agregados <- data.frame (rhea=dados_agregados,
                               campo = habitat_campo,
                               agri=habitat_lavoura)

# modelo 3
model3 <- glm (rhea ~ campo+agri,
               data=dados_agregados,
               family=binomial(link="logit"))

# salvar
save (model3, file=here("output","model3_glm.RData"))
summary(model3)

# ----------------------------------------------------------------



# modelo 3 em BuGS

# escrever modelo
sink("model3_aggregation_NoDetection.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################

    # priors
    BETA0 ~ dnorm(0,0.001)
    BETA.CAMPO ~ dnorm(0,0.001)
    BETA.AGRI ~ dnorm(0,0.001)
    
    # LIKELIHOOD
    for (i in 1:nsite) {

      y1[i] ~ dbern(psi[i]) # Assuming y = z at site i
      mu[i] <- BETA0 + BETA.CAMPO * campo[i] + # occ model
                       BETA.AGRI * agri [i] 
      mu.lim[i] <- min(10, max(-10, mu[i]))  
      logit(psi[i]) <- mu.lim[i]

    }    

    # ======================================
    # goodness of fit
    # posterior predictive checks
    for (i in 1:nsite) {
    
        residual[i] <- y1[i]-psi[i] # Residuals for observed data
        predicted[i] <- psi[i] # Predicted values
        sq[i] <- pow(residual[i], 2) # Squared residuals for observed data
       
        # Generate replicate data and compute fit stats for them
        y.new[i] ~ dbern(psi[i]) # one new data set at each MCMC iteration
        sq.new[i] <- pow(y.new[i]-predicted[i], 2) # Squared residuals for new data
    }
        
    fit <- sum(sq[]) # Sum of squared residuals for actual data set
    fit.new <- sum(sq.new[]) # Sum of squared residuals for new data set
    test <- step(fit.new - fit) # Test whether new data set more extreme
    bpvalue <- mean(test) # Bayesian p-value
    
    # proportion of occupied municipalites
    pm <- sum(y1[])/nsite

    }## end of the model
    ",fill = TRUE)
sink()


## initial values
set.seed(32)
inits = function() {list(BETA0=runif(1),
                         BETA.CAMPO=runif(1),
                         BETA.AGRI=runif(1))}


# bound data
str(win.data <- list(nsite= nrow (dados_agregados),
                     y1 = dados_agregados$rhea,
                     campo = dados_agregados$campo,
                     agri = dados_agregados$agri
))


params <- c("BETA0", 
            "BETA.CAMPO","BETA.AGRI",
            "psi","pm",
            "fit",
            "fit.new",
            "bpvalue"
)

model3_bugs <- bugs(data = win.data, 
               parameters.to.save = params, 
               model.file = "model3_aggregation_NoDetection.txt", 
               inits = NULL, 
               n.chains = nc, 
               n.thin = nt, 
               n.iter = ni, 
               n.burnin = nb,
               codaPkg=F, 
               debug=T,
               DIC=TRUE, 
               bugs.directory="C:/Program Files/WinBUGS14/", 
               program= "WinBUGS"
)

# salvar
save(model3_bugs,
     file=here("output","model3_bugs.RData"))

# model fit
model3_bugs$mean$pm
model3_bugs$mean$bpvalue


# ------------------------------------------------------
# modelo 4
# com integracao e deteccao

# write model
sink("model4_integration_grassland_agriculture.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################

    # PRIORS
    BETA0 ~ dunif(-10,10)
    BETA.CAMPO ~ dunif(-10, 10)
    BETA.AGRI ~ dunif(-10, 10)

    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
  
      mu[i] <- BETA0 + BETA.CAMPO * campo [i] + BETA.AGRI * agri [i]

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
        zP1[i] <- P1[i] * z [i]
        y.gbif [i] ~ dbern(zP1[i])


    }

    # PRIOR FOR RICHNESS EFFECT
    ALPHA.GBIF ~ dnorm(0,0.0001)I(0,1000)# truncated in positive values
    
    ## EBIRD
  
    for (i in 1:nsite) {
      for (j in 1:nrepEB) { # n checklists

        e2[i,j] <- ALPHA.DIST * dist [i,j] +
                   ALPHA.DURATION * duration [i,j] +
                   ALPHA.OBSERVER * observer [i,j]
 
        P2[i,j] <- 1-pow((1-0.5), e2[i,j])
        zP2[i,j] <- P2[i,j] * z [i]
        yEB [i,j] ~ dbern(zP2[i,j])

      }
    }
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,1000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,2000)
    ALPHA.OBSERVER  ~ dnorm(0,0.0001)I(0,1000)

    ### WIKIAVES
    
    for (i in 1:nsite) {
    
        e5[i] <- ALPHA.PICT * nPIC.wikiaves [i] + ALPHA.SONG * nSONG.wikiaves [i]
        P5[i] <- 1-pow((1-0.5), e5[i])
        zP5[i] <- P5[i] * z [i]
        y.wikiaves [i] ~ dbern (zP5[i])
    
    }
    
    # PRIOR FOR RICHNESS EFFECT
    ALPHA.PICT ~ dnorm(0,0.0001)I(0,30000)# truncated in positive values
    ALPHA.SONG ~ dnorm(0,0.0001)I(0,10000)# truncated in positive values
    

    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    #compute the mean detection probability of each dataset: 
    muP1 <- mean(P1[])
    muP2 <- mean(P2[,])
    muP5 <- mean(P5[])

    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    }## end of the model
    ",fill = TRUE)
sink()

# bound data
str(win.data <- list(nsite= dim(dados_det_ema_gbif)[1],
                     y.gbif = dados_det_ema_gbif[,"det"],
                     nSP.gbif = dados_det_ema_gbif [,"riqueza_aves"], 
                     y.wikiaves = dados_wikiaves [, "RHAMERICANA"],
                     nPIC.wikiaves = dados_wikiaves [, "NPIC"],
                     nSONG.wikiaves = dados_wikiaves [, "NSONG"],
                     yEB = y.ebird [,1:n_so_ebird],
                     dist = dist.ebird [,1:n_so_ebird],
                     duration = dist.duration [,1:n_so_ebird],
                     observer = dist.observers [,1:n_so_ebird],
                     nrepEB = ncol(dist.observers [,1:n_so_ebird]),
                     campo=habitat_campo[,1],
                     agri = habitat_lavoura[,1]
))

# initial values
set.seed(32)
inits <- function() {list(z = rep (1, nrow(y.ebird)),
                         ALPHA.DIST = rnorm (1,mean=2),
                         ALPHA.DURATION = rnorm (1,mean=20),
                         ALPHA.OBSERVER = rnorm (1,mean=2))}


params <- c("BETA0", 
            "BETA.CAMPO","BETA.AGRI",
            "ALPHA.GBIF","ALPHA.DIST", "ALPHA.DURATION","ALPHA.OBSERVER",
            "ALPHA.PICT", "ALPHA.SONG", 
            "muP1","muP2","muP5",
            "z","psi","fs.z"
)

model4 <- bugs(data = win.data, 
                              parameters.to.save = params, 
                              model.file = "model4_integration_grassland_agriculture.txt", 
                              inits = inits, 
                              n.chains = nc, 
                              n.thin = nt, 
                              n.iter = ni, 
                              n.burnin = nb,
                              codaPkg=F, 
                              DIC=TRUE, 
                              debug=F,
                              bugs.directory="C:/Program Files/WinBUGS14/", 
                              program= "WinBUGS"
                              )

# save it
save (model4,file=here("output", "model4_bugs.RData"))


model4$mean$fs.z

# Pstar
source ('R/4_Pstar.R')



# -------------------------------------------------------------------------------



# example using the sppOccupancy package of R

# bound data

y <- list (dados_det_ema_gbif[,"det"],
           dados_wikiaves [, "RHAMERICANA"],
           #ifelse (apply (y.ebird,1,sum,na.rm=T)>0,1,apply (y.ebird,1,sum,na.rm=T))
           y.ebird [which(rowSums(is.na(y.ebird[,1:11])) < 11),1:11]
           )
# list of sites           
sites <- list (seq (1,length(rownames(dados_det_ema_gbif))),
               seq (1,length(rownames(dados_det_ema_gbif))),
               which(rowSums(is.na(y.ebird[,1:11])) < 11))

# detection covariates
det.covs <- list ()
det.covs [[1]] <- list (nSP.gbif = as.numeric(scale (dados_det_ema_gbif [,"riqueza_aves"])))
det.covs [[2]] <- list (nPIC.wikiaves = as.numeric(scale (dados_wikiaves [, "NPIC"])),
                        nSONG.wikiaves = as.numeric(scale (dados_wikiaves [, "NSONG"])))
#det.covs [[3]] <- list (duration = as.numeric(scale (apply (dist.duration,1,sum) )))
det.covs [[3]] <- list (duration = (scale ( (dist.duration) ))[which(rowSums(is.na(y.ebird[,1:11])) < 11),1:11],
                        dist = (scale ( (dist.ebird) ))[which(rowSums(is.na(y.ebird[,1:11])) < 11),1:11],
                        observer = (scale ( (dist.observers) ))[which(rowSums(is.na(y.ebird[,1:11])) < 11),1:11])
# site cov

site.cov <- as.matrix (cbind (habitat_campo[,1], habitat_lavoura[,1]))
colnames (site.cov)<-c("grassland", "agriculture")

# list of data
data.list <- list (y = y,
                   occ.covs = site.cov,
                   det.covs = det.covs,
                   sites = sites)

# initial values
set.seed(32)
inits <- list(z = rep (1, length(sites[[1]])),
              beta = list (2),
              alpha = list (rnorm (1,mean=2),
                            rnorm (1,mean=20),
                            rnorm (1,mean=2))
              )

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72)))
# using the spOccupancy package of R

model4_spOcc <- intPGOcc(occ.formula = ~ grassland+agriculture, 
                         det.formula = list (f.1 = ~ nSP.gbif,
                                             f.2 = ~ nPIC.wikiaves+nSONG.wikiaves,
                                             f.3 = ~ duration+dist+observer), 
                         data = data.list, 
                         inits = inits, 
                         priors = prior.list, 
                         n.samples =ni, 
                         n.omp.threads = nc, 
                         verbose = TRUE, 
                         n.burn = nb, 
                         n.thin = nt, 
                         n.chains = nc)

summary (model4_spOcc)

# save it
save (model4_spOcc,file=here("output", 
                             "model4_spOcc.RData"))



# -----------------------
# model com implementacao de validacao cruzada para avaliacao de ajuste do modelo
##### Dividir dados em conjunto de teste/treino(validacao)

## Validation taking out part of the munis from all data sets: Model 4
### Sample 20% of the municipalities ###
sampled <- sample_n(ID_MUN, size=100, replace=F) #aprox 20% of the muni, sample without replacement
#for each data set, take out the samples from munis that are in 'sampled' object:

# GBIF
datGBOut <- dados_det_ema_gbif[dados_det_ema_gbif$mun %in% sampled$NM_MUNICIP, ]
datGBIn <- dados_det_ema_gbif[-which(dados_det_ema_gbif$mun %in% sampled$NM_MUNICIP), ]

# WIKIAVES
datWAOut <- dados_wikiaves[dados_wikiaves$MUNI %in% sampled$NM_MUNICIP, ]
datWAIn <- dados_wikiaves[-which(dados_wikiaves$MUNI %in% sampled$NM_MUNICIP), ]

# EBIRD
row.names(y.ebird) <- ID_MUN$NM_MUNICIP
datEBOut <- y.ebird[rownames(y.ebird) %in% sampled$NM_MUNICIP, ]
datEBIn <- y.ebird[-which(rownames(y.ebird) %in% sampled$NM_MUNICIP),]

# Ebird effort covariates:
# distancia
row.names(dist.ebird) <- ID_MUN$NM_MUNICIP
dist.ebirdOUT <- dist.ebird[rownames(dist.ebird) %in% sampled$NM_MUNICIP, ]
dist.ebird <- dist.ebird[-which(rownames(dist.ebird) %in% sampled$NM_MUNICIP),]
# duracao
row.names(dist.duration) <- ID_MUN$NM_MUNICIP
dist.durationOUT <- dist.duration[rownames(dist.duration) %in% sampled$NM_MUNICIP, ]
dist.duration <- dist.duration[-which(rownames(dist.duration) %in% sampled$NM_MUNICIP),]
# n observadores
row.names(dist.observers) <- ID_MUN$NM_MUNICIP
dist.observersOUT <- dist.observers[rownames(dist.observers) %in% sampled$NM_MUNICIP, ]
dist.observers <- dist.observers[-which(rownames(dist.observers) %in% sampled$NM_MUNICIP),]


# write model
sink("model4_integration_grassland_agriculture_cross.txt")
cat("

    model {
    
    ###################################
    ###### ECOLOGICAL MODEL ###########
    ###################################

    # PRIORS
    BETA0 ~ dunif(-10,10)
    BETA.CAMPO ~ dunif(-10, 10)
    BETA.AGRI ~ dunif(-10, 10)

    # LIKELIHOOD
    for (i in 1:nsite) {

      z[i] ~ dbern(psi[i]) # True occupancy z at site i
  
     mu[i] <- BETA0 + BETA.CAMPO * campo [i] + BETA.AGRI * agri [i]

     # Keeping Winbugs on the track
     mu.lim[i] <- min(10, max(-10, mu[i]))  
     logit(psi[i]) <- mu.lim[i]

    }    

    ###################################
    ###### OBSERVATION MODELS #########
    ###################################

    ### GBIF 
    
    #treino - IN
    for (i in 1:nsiteIN1) {

        e1[i] <- ALPHA.GBIF * nSP.gbif[i]
        P1[i] <- 1-pow((1-0.5), e1[i])
        zP1[i] <- P1[i] * z [i]
        y.gbif [i] ~ dbern(zP1[i])

    }
    
    #teste - OUT
    for (i in 1:nsiteOUT1) {

        e6[i] <- ALPHA.GBIF * nSP.gbif6[i]
        P6[i] <- 1-pow((1-0.5), e6[i])
        zP6[i] <- P6[i] * z [i]
        y.gbif6 [i] ~ dbern(zP6[i])

    }


    # PRIOR FOR RICHNESS EFFECT
    ALPHA.GBIF ~ dnorm(0,0.0001)I(0,1000)# truncated in positive values
    
    ### WIKIAVES
    #treino IN:
    
    for (i in 1:nsiteIN4) {
    
        e5[i] <- ALPHA.PICT * nPIC.wikiaves[i]+ALPHA.SONG*nSONG.wikiaves[i]
        P5[i] <- 1-pow((1-0.5), e5[i])
        zP5[i] <- P5[i] * z [i]
        y.wikiaves[i] ~ dbern(zP5[i])#+(1-zP5[i])*q[i])
        
    }
    
    #teste OUT:
    
    for (i in 1:nsiteOUT4) {
    
        e9[i] <- ALPHA.PICT * nPIC.wikiaves9[i]+ALPHA.SONG*nSONG.wikiaves9[i]
        P9[i] <- 1-pow((1-0.5), e9[i])
        zP9[i] <- P9[i] * z [i]
        y.wikiaves9[i] ~ dbern(zP9[i])#+(1-zP9[i])*q[i])
        
    }
    
    # PRIOR FOR RICHNESS EFFECT
    ALPHA.PICT ~ dnorm(0,0.0001)I(0,30000)# truncated in positive values
    ALPHA.SONG ~ dnorm(0,0.0001)I(0,10000)# truncated in positive values
  
  ## EBIRD
  
    for (i in 1:nsiteIN5) {
      for (j in 1:nrepEBin) {

        e2[i,j] <- ALPHA.DIST * dist[i,j] +
                   ALPHA.DURATION * duration[i,j] +
                   ALPHA.OBSERVER * observer[i,j]
 
        P2[i,j] <- 1-pow((1-0.5), e2[i,j])
        zP2[i,j] <- P2[i,j] * z [i]
        yEB [i,j] ~ dbern(zP2[i,j])

      }
    }
    
    for (i in 1:nsiteOUT5) {
      for (j in 1:nrepEBout) {

        e10[i,j] <- ALPHA.DIST * dist10[i,j] +
                   ALPHA.DURATION * duration10[i,j] +
                   ALPHA.OBSERVER * observer10[i,j]
 
        P10[i,j] <- 1-pow((1-0.5), e10[i,j])
        zP10[i,j] <- P10[i,j] * z [i]
        yEB10 [i,j] ~ dbern(zP10[i,j])

      }
    }
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,1000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,2000)
    ALPHA.OBSERVER  ~ dnorm(0,0.0001)I(0,1000)

    ##############################
    #### DERIVED PARAMETERS ######
    ##############################
    
    #compute the mean detection probability of each dataset: 
    muP1 <- mean(P1[])
    muP2 <- mean(P2[,])
    muP5 <- mean(P5[])

    ## Number of occupied municipalities (finite sample size)
    fs.z <- sum(z[])/nsite
    

    }## end of the model
    ",fill = TRUE)
sink()

## initial values
inits = function() {list(z = rep (1, nrow(y.ebird)),
                         ALPHA.GBIF= rnorm (1,mean=10),
                         ALPHA.DIST = rnorm (1,mean=2),
                         ALPHA.DURATION = rnorm (1,mean=20),
                         ALPHA.OBSERVER = rnorm (1,mean=2),
                         ALPHA.PICT = rnorm (1,mean=10),
                         ALPHA.SONG = rnorm (1,mean=10))}

## parametros a serem monitorados
params <- c("BETA0", "BETA.CAMPO", "ALPHA.GBIF","ALPHA.DIST", 
                   "ALPHA.DURATION","ALPHA.OBSERVER",
                   "ALPHA.PICT", "ALPHA.SONG", 
                   "muP1","muP2","muP5", 
                   "z","psi","fs.z", "y.gbif6",
                   "y.wikiaves9", "yEB10")

# dados
(win.data <- list(nsite = dim(dados_det_ema_gbif)[1],
                          y.gbif = datGBIn[,"det"],
                          nsiteIN1 = dim(datGBIn)[1],
                          nsiteOUT1 = dim(datGBOut)[1],
                          nSP.gbif = datGBIn[,"riqueza_aves"],
                          nSP.gbif6 = datGBOut[,"riqueza_aves"],
                          nsiteIN4 = dim(datWAIn)[1],
                          nsiteOUT4 = dim(datWAOut)[1],
                          y.wikiaves = datWAIn[, "RHAMERICANA"],
                          nPIC.wikiaves = datWAIn[, "NPIC"],
                          nSONG.wikiaves = datWAIn[, "NSONG"],
                          nPIC.wikiaves9 = datWAOut[, "NPIC"],
                          nSONG.wikiaves9 = datWAOut[, "NSONG"],
                          yEB = datEBIn [,1:n_so_ebird],
                          nsiteOUT5 = dim(datEBOut)[1],
                          nsiteIN5 = dim(datEBIn)[1],
                          dist = dist.ebird [,1:n_so_ebird],
                          duration = dist.duration [,1:n_so_ebird],
                          observer = dist.observers [,1:n_so_ebird],
                          dist10 = dist.ebirdOUT [,1:n_so_ebird],
                          duration10 = dist.durationOUT [,1:n_so_ebird],
                          observer10 = dist.observersOUT [,1:n_so_ebird],
                          nrepEBin = ncol(dist.observers [,1:n_so_ebird]),
                          nrepEBout = ncol(dist.observersOUT [,1:n_so_ebird]),
                          campo=habitat_campo[,1],
                          agri = habitat_lavoura[,1]
  ))
  
# run model
model4_cross <- bugs(data = win.data, 
                     parameters.to.save = params, 
                     model.file = "model4_integration_grassland_agriculture_cross.txt", 
                     inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
                     codaPkg=F, DIC=TRUE, debug=T,
                     bugs.directory="C:/Program Files/WinBUGS14/", program= "WinBUGS")
save(model4_cross,
     sampled, # é bom salvar a porcao dos dados utilizados para a validacao do modelo
     file=here ("output", "model4_cross.RData"))

# -----------------------
#### Comparando performance de diferentes modelos ####

## 1) Calcular a Deviance
# separar os dados observados (teste) - GBIF
Ytruth <- datGBOut[,"det"]
# separar os dados preditos/estimados pelo modelo - GBIF
Yhat <- model4_cross$mean$y.gbif6

# separar os dados observados (teste) - Wikiaves
Ytruth2 <- datWAOut[, "RHAMERICANA"]
# separar os dados preditos/estimados pelo modelo - Wikiaves
Yhat2 <- model4_cross$mean$y.wikiaves9

# separar os dados observados (teste) - eBird
Ytruth3 <- datEBOut[,1:n_so_ebird]
# separar os dados preditos/estimados pelo modelo - eBird
Yhat3 <- model4_cross$mean$yEB10

## Deviance em cada base de dados:
# GBIF
likhood <- (Yhat^Ytruth)*((1-Yhat)^(1-Ytruth))
DEV <- -(2*(sum(log(likhood))))

# Wikiaves
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
lokLik2<- log(likhood2)
lokLik2 <- lokLik2[is.infinite(lokLik2) == F]
DEV2 <- -(2*(sum(lokLik2)))

# eBird
likhood3 <- (Yhat3^Ytruth3)*((1-Yhat3)^(1-Ytruth3))
lokLik3<- log(likhood3)
lokLik3 <- lokLik3[is.na(lokLik3) == F]
DEV3 <- -(2*(sum(lokLik3)))

# Deviance total:
DEVtotal <- DEV + DEV2 + DEV3 # valor que deve ser comparado entre modelos.
                              # modelo com o menor valor de deviance tem 
                              # o melhor ajuste aos dados.


## 2) Plotar curvas AUC

# gbif
roc(Ytruth, Yhat)
# wikiaves
roc(Ytruth2, Yhat2)
# eBird
roc(as.numeric(Ytruth3), as.numeric(Yhat3))

# plot GBIF
plot(roc(Ytruth, Yhat), main = "Full Model")
mtext(text = paste ("AUC=",round (roc(estpred[,1], estpred[,2])$auc,2)), side=1)
# plot WikiAves
plot(roc(Ytruth2, Yhat2), main = "Full Model")
mtext(text = paste ("AUC=",round (roc(Ytruth2, Yhat2)$auc,2)), side=1)

# plot eBird
plot(roc(as.numeric(Ytruth3), as.numeric(Yhat3)), main = "Full Model")
mtext(text = paste ("AUC=",round (roc(as.numeric(Ytruth3), as.numeric(Yhat3))$auc,2)), side=1)


rm(list=ls())
