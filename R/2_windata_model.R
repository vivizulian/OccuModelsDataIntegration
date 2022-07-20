
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
## area de habitat campestre e agricutlura por municipio
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

# bind grassland data into the shapefile
shape_RS@data$campo <- area_campo
# plot
cores_campo <- data.frame (cores= sqrt(area_campo),
                         NM_MUNICIP=shape_RS$NM_MUNICIP)

# fortify
f.mun<-fortify(shape_RS, region="NM_MUNICIP")
f.mun_campo<- cbind (f.mun, 
                   Nespecies = cores_campo [match (f.mun$id, cores_campo$NM_MUNICIP),]$cores)

# grassland cover
# inserir vals
c_campo <-   ggplot() + geom_polygon(data=f.mun_campo, aes(x=long, y=lat, group=group, 
                                                color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Area de campo, por municipio") +
  scale_fill_gradient2 (low='white', high='#206A5D', na.value = "white",
                        limits=c(0,max(cores_campo$cores)), 
                        breaks=seq(0,max(cores_campo$cores,na.rm=T),by=500),
                        name=expression(sqrt("Area (km2)"))) ## para continuo

(c_campo <- c_campo + theme_classic() + 
    theme (axis.text = element_text(size=6),
           axis.title = element_text(size=8),
           legend.text = element_text(size=8),
           legend.title = element_text(size=9))+
    xlab("Longitude") + 
    ylab("Latitude")) 

## agricultura
# bind grassland data into the shape
shape_RS@data$agricultura <- area_lavoura
# plot
cores_agri <- data.frame (cores= sqrt(area_lavoura),
                         NM_MUNICIP=shape_RS$NM_MUNICIP)

# fortify
f.mun_agri<- cbind (f.mun, 
                   Nespecies= cores_agri [match (f.mun$id, cores_agri$NM_MUNICIP),]$cores)

## inserir vals
c_agri <-   ggplot() + geom_polygon(data=f.mun_agri, aes(x=long, y=lat, group=group, 
                                                       color=Nespecies, fill=Nespecies), 
                                   colour = NA, size=1) + 
  labs (title= "Area de agricultura, por municipio") +
  scale_fill_gradient2 (low='white', high='#E48900', na.value = "white",
                        limits=c(0,max(cores_agri$cores)), 
                        breaks=seq(0,max(cores_agri$cores,na.rm=T),by=500),
                        name=expression(sqrt("Area (km2)"))) ## para continuo

(c_agri<-c_agri + theme_classic() + 
    theme (axis.text = element_text(size=6),
           axis.title = element_text(size=8),
           legend.text = element_text(size=8),
           legend.title = element_text(size=9))+
  xlab("Longitude") + 
  ylab("Latitude"))

# bind maps and save
pdf (here ("output", "maps_land_use.pdf"),heigh=7,width=5)
grid.arrange(c_campo,
             c_agri)

dev.off()

# -------------------------------------------------------------------------
# MAPAS DAS OBSERVACOES

# carregar mapa do RS
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

# fortify RS map
f.mun<-fortify(shape_RS, region="NM_MUNICIP") # fortify mapa do RS, comum a todos os mapas

## colocar o shape da america do sul = comum a todos os mapas
a <- ggplot() + geom_polygon (data=BR_AR_URU, aes(x=long, y=lat, group=group),size = 0.1, fill="gray90", colour="gray75",alpha=1) +
  coord_fixed (xlim = c(-57.5, -49),  ylim = c(-34, -27), ratio = 1) 

## inserir os lagos = comum a todos os mapas
b <- a + geom_polygon (data=lagos,aes(x=long, y=lat, group=group), 
                       fill="lightcyan",colour = "lightcyan", size=1)

## abrir deteccoes do eBird
load (here("data","organized_data",  "input_ebird.RData"))## gbif
ebird_rhea <- apply (y.ebird, 1, max,na.rm=T) # 

## mapas das deteccoes
cores_ebird <- data.frame (cores= ebird_rhea,
                           NM_MUNICIP=shape_RS$NM_MUNICIP)
# escala discreta
cores_ebird$cores <- factor(cores_ebird$cores)
levels(cores_ebird$cores) [which(levels (cores_ebird$cores) == 0)] <- "Not detected"
levels(cores_ebird$cores) [which(levels (cores_ebird$cores) == 1)] <- "Detected"
levels(cores_ebird$cores) [which(levels (cores_ebird$cores) == -Inf)] <- "Not sampled"

# fortify
f.mun_ebird<- cbind (f.mun, 
                     Nespecies= cores_ebird [match (f.mun$id, 
                                                    cores_ebird$NM_MUNICIP),]$cores)

## inserir deteccoes do eBird
c_ebird <-   b + geom_polygon(data=f.mun_ebird, aes(x=long, y=lat, group=group, 
                                                    color=Nespecies, fill=Nespecies), 
                              colour=NA,size=1) + 
  labs (title= "eBird") +
  scale_fill_manual("Observation data",
                    values = c("Not detected" = "white",
                               "Detected" = "darkred",
                               "Not sampled" = "gray85"))

# anotar
e_ebird <- c_ebird + ggsn::scalebar(f.mun_ebird, dist = 100, st.dist=0.03,st.size=2.2, height=0.02, 
                                    transform = TRUE, dist_unit = "km",
                                    model = 'WGS84', location = "bottomright")

## plot para extrair a legenda
f_ebird_legend <- e_ebird +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        
        legend.title = element_text(size=9),
        legend.text = element_text(size=7),
        legend.key.width=unit(0.85,"cm"),
        legend.key.size = unit(0.40,"cm"),
        legend.position = "top",
        legend.justification = 0.5,
        legend.direction="horizontal",
        legend.box="horizontal",
        axis.text=element_text(size=3),
        axis.text.x = element_text(size=3),
        axis.title.x = element_text(size = 5),
        axis.text.y = element_text(size=3),
        axis.title.y = element_text(size = 5),
        plot.title = element_text(size=8),
        plot.margin = unit(c(0.1, -0.1,-0.1, 0.2), "lines")) 
#        legend.margin = margin (0,0,0,0),
#       legend.box.margin = margin(1,0,0,0)) +

f_ebird_legend

## funcao para capturar legenda
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## extrair legenda usando a funcao anterior
legenda_comum_data <- get_legend(f_ebird_legend)

## plot para o painel (sem a legenda)  
f_ebird <- e_ebird + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        plot.title = element_text(size=9),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.size = unit(0.5,"cm"),
        axis.text=element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size = 8),
        plot.margin = unit(c(-0.1,-0.1, -3.9, 0.05), "lines"))
# top, right, bottom, and left margins

# anota o norte
f_north_ebird <- f_ebird + ggsn::north(f.mun_ebird, symbol=1,scale = 0.2,location = "bottomleft") +
  theme(legend.text=element_text(size=7),
        legend.title=element_text(size=8))

f_north_ebird

## abrir deteccoes GIBF
load (here("data","organized_data", "input_GBIF.RData"))## gbif

## mapas das deteccoes
cores_gbif <- data.frame (cores= dados_det_ema_gbif$det,
                          NM_MUNICIP=shape_RS$NM_MUNICIP)
cores_gbif$cores [which (dados_det_ema_gbif$riqueza_aves == 0)] <- NA

# fortify
f.mun_gbif<- cbind (f.mun, 
                    Nespecies= cores_gbif [match (f.mun$id, 
                                                  cores_gbif$NM_MUNICIP),]$cores)

## inserir deteccoes do GBIF
c_gbif <-   b + geom_polygon(data=f.mun_gbif, aes(x=long, y=lat, group=group, 
                                                  color=Nespecies, fill=Nespecies), 
                             colour = NA, size=1) + 
  labs (title= "GBIF") +
  scale_fill_gradient2 (low='white', 
                        high='darkred', 
                        midpoint= 0.20,na.value = "gray85",
                        limits=c(0,1), 
                        #breaks=seq(0,max(cores_50km$cores,na.rm=T),by=0.2),
                        name="Detected") ## para continuo


f_gbif <- c_gbif + 
  xlab("") + ylab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        plot.title = element_text(size=9),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=3),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 4),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 4),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.2,-0.1, -3, -0.1), "lines")) # top, right, bottom, and left margins


## abrir deteccoes wikiaves
load (here("data","organized_data", "input_wikiaves.RData"))

## mapas das deteccoes
cores_wiki <- data.frame (cores= dados_wikiaves$RHAMERICANA,
                          NM_MUNICIP=shape_RS$NM_MUNICIP)
cores_wiki$cores [which (dados_wikiaves$NSPECIES == 0)] <- NA

# fortity
f.mun_wiki <- cbind (f.mun, 
                     Nespecies= cores_wiki [match (f.mun$id, 
                                                   cores_wiki$NM_MUNICIP),]$cores)

## inserir deteccoes do GBIF
c_wiki <-   b + geom_polygon(data=f.mun_wiki, aes(x=long, y=lat, group=group, 
                                                  color=Nespecies, fill=Nespecies), 
                             colour = NA, size=1) + 
  labs (title= "Wikiaves") +
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "gray85",
                        limits=c(0,1), 
                        #breaks=seq(0,max(cores_50km$cores,na.rm=T),by=0.2),
                        name="Detected") ## para continuo

f_wiki <- c_wiki +
  xlab("") + ylab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        plot.title = element_text(size=9),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=3),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 4),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 4),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.2,0, -3, -0.1), "lines"))# top, right, bottom, and left margins

## concenso entre as bases

cores_agregado <- data.frame(eBird=ebird_rhea, 
                             GBIF=cores_gbif$cores, 
                             WikiAves=cores_wiki$cores)

apply(cores_agregado, 2, function (i) sum (i > 0,na.rm=T))## numero de deteccoes por base
table(apply (cores_agregado,1,max,na.rm=T)>0) # numero de municipios com deteccao e sem deteccao

# cores para o mapa agregado
cores_agregado <- data.frame (cores= apply (cores_agregado,1,max,na.rm=T),
                              NM_MUNICIP=shape_RS$NM_MUNICIP)

#fortify
f.mun_con <- cbind (f.mun, 
                    Nespecies= cores_agregado [match (f.mun$id, 
                                                      cores_agregado$NM_MUNICIP),]$cores)

# mapa
c_con <-   b + geom_polygon(data=f.mun_con, aes(x=long, y=lat, group=group, 
                                                color=Nespecies, fill=Nespecies), 
                            colour = NA, size=1) + 
  labs (title= "Agregado das bases") +
  scale_fill_gradient2 (low='white', 
                        high='darkred', 
                        midpoint= 0.20,na.value = "gray85",
                        limits=c(0,1), 
                        #breaks=seq(0,max(cores_50km$cores,na.rm=T),by=0.2),
                        name="Detected") ## para continuo
f_con<- c_con +
  xlab("") + ylab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        plot.title = element_text(size=9),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=3),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 4),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 4),
        axis.ticks = element_blank(),
        plot.margin = unit(c(-5,0, 0, -0.1), "lines"))# top, right, bottom, and left margins

######### arranjar o painel
pdf(file=here ("output","maps_observations.pdf"),width = 7,height = 5,family="serif")

grid.arrange(f_north_ebird, f_gbif,f_wiki,f_con,
             legenda_comum_data,
             ncol=4,nrow=2,
             layout_matrix = rbind(c(1,2,3,4), 
                                   c(5,5,5,5))) 
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

# ------------------------------------------------------ #
# modelo 2 em linguagem bugs

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

    ## EBIRD
  
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
    
    ALPHA.DIST ~ dnorm(0,0.0001)I(0,1000)
    ALPHA.DURATION ~ dnorm(0,0.0001)I(0,2000)
    ALPHA.OBSERVER  ~ dnorm(0,0.0001)I(0,1000)

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
            "muP","z","psi","fs.z","P2"
    )

# GLOBAL MCMC settings
## short form
na <- 30; nb <- 40; ni <- 50; nc <- 3; nt <- 1
# na <- 3000; nb <- 4000; ni <- 5000; nc <- 3; nt <- 1

# run model
out_model2 <- bugs(data = winbugs.data, 
                  parameters.to.save = params, 
                  model.file = "model2_NoIntegration_Detection.txt", 
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
save (out_model2,file=here("output", "out_model2_bugs.RData"))


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

# --------
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
            "psi","pm"
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
inits = function() {list(z = rep (1, nrow(y.ebird)),
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

# Pstar
source ('R/4_Pstar.R')

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
