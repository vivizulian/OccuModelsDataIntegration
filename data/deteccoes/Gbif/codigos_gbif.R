# --------------------------------------------------------------------------- #
#                     Dados GBIF
# --------------------------------------------------------------------------- #

# carregar os pacotes
source ('R/packagesR.R')

# carregar
load(here ("data","deteccoes","Gbif","GBIF_aves_RS_filtrado.RData"))

## remover dados antigos  
## range EBIRD "2008-01-25" "2018-12-31"
dados_aves_filtrado <- dados_aves_filtrado [which (as.Date (as.character(dados_aves_filtrado$eventDate)) >= as.Date ("2008-01-01")),]
## remover dados muito recentes - para fechar com ebird 
dados_aves_filtrado <- dados_aves_filtrado [which (as.Date (as.character(dados_aves_filtrado$eventDate)) <= as.Date ("2018-12-31")),]

## carregar mapa do RS
shape_RS <- readOGR(dsn=here("data","shape_munRS"), layer="43MUE250GC_SIR",
                    encoding = "UTF-8",use_iconv = T)
shape_RS <- shape_RS [-c(96,250),] ## remover os lagos

### remover coordenadas fora do RS
dados_aves_pais <- dados_aves_filtrado[which(dados_aves_filtrado$decimalLatitude > extent(shape_RS)[3] & dados_aves_filtrado$decimalLatitude < extent(shape_RS)[4]),]
dados_aves_pais <- dados_aves_pais[which(dados_aves_pais$decimalLongitude > extent(shape_RS)[1] & dados_aves_pais$decimalLongitude < extent(shape_RS)[2]),]

### gerando points dataframe com os dados
coordenadas <- dados_aves_pais [,c("decimalLongitude", "decimalLatitude")]
coordinates (coordenadas) <- ~ decimalLongitude + decimalLatitude
crs(coordenadas)<- crs(shape_RS)

## sobrepor no shape do RS
over_coord_RS <- over (coordenadas,shape_RS)
dados_aves_pais <-dados_aves_pais [which(is.na(over_coord_RS$NM_MUNICIP)==F),]
 
## dados de Rhea americana ( ou de qualeur outra spp que vcs quiserem)
dados_rhea <- dados_aves_pais [which(dados_aves_pais$species == "Rhea americana"),]

### gerando spatialpoints dataframe com os dados de EMA, para sobrepor com o shape do RS
coordenadas_ema <- dados_rhea [,c("decimalLongitude", "decimalLatitude")]
coordinates (coordenadas_ema) <- ~ decimalLongitude + decimalLatitude
crs(coordenadas_ema)<- crs(shape_RS) # coordinate reference system

# obter municipios com registro de ema
mun_ema <- over (coordenadas_ema,shape_RS)

# obter um mapa destes municipios com  deteccao da ema
shape_mun_ema <- shape_RS [which (shape_RS@data$NM_MUNICIP %in% mun_ema$NM_MUNICIP),]

#plot(shape_RS,border="gray60",main="Dados GBIF")
plot(shape_RS,
     col = rgb(red = ifelse(shape_RS$NM_MUNICIP %in% mun_ema$NM_MUNICIP,1,0), 
                             green = 0, blue = 1, alpha = 0.5),
     border = "black",
     lwd=1,
     main="GBIF Data")

legend ('bottomleft',legend =c("Detection of Rhea",
                               "No detection of Rhea"),
        col = c("purple","blue"),pch=15,bty="n")

# plot de todos os pontos de amostragem, e dos pontos com Rhea 

points(dados_aves_pais$decimalLongitude,
       dados_aves_pais$decimalLatitude,
       pch=1,cex=0.7,col="gray40")

# pontos com RHea
points(dados_rhea$decimalLongitude,dados_rhea$decimalLatitude,
       col="red",pch=1,cex=0.7)## so pra dar uma engrossada na espessura da linha
points(dados_rhea$decimalLongitude,dados_rhea$decimalLatitude,
       col="red",pch=1,cex=0.7)

### obter um dataframe, 
#        com o cod do mun na linha, 
#        e a deteccao \ nao det de ema na coluna

det_ema <- data.frame (cod_mun = shape_RS@data$CD_GEOCMU,
                       det = ifelse (shape_RS@data$CD_GEOCMU %in% shape_mun_ema@data$CD_GEOCMU,1,0),
                       mun = shape_RS@data$NM_MUNICIP)

## abaixo, comando para conferir de seu certo
# det_ema$mun[which(det_ema$det == 1)] == shape_mun_ema@data$NM_MUNICIP

## agora, obter um dataframe com o numero total de especies de aves registradas por municipio do RS
## para isto, criar um spatialpoint dataframe para as coordenadas de cada especie
especies_aves_RS <- unique (dados_aves_pais$species) ## pegar a identidade de cada spp no banco de dados filtrado

## entao fazer um conjunto dos dados completos, para cada spp, e pegar os pontos com o registro 
subconjunto_aves_RS <- lapply (especies_aves_RS, function (especie)
                    
                                    dados_aves_pais [which (dados_aves_pais$species == especie),]
                              
                              )

## subset do conjunto de coordenadas de cada sp. 
coordenadas_especies <- lapply (subconjunto_aves_RS, function (especie)
    
                                    especie [,c("decimalLongitude", "decimalLatitude")]
                                )

## transformar cada objeto da lista de coordenadas de cada especie
# em spatialpoints dataframe

# PS: LAPPLY COM CHAVES

coordenadas_especies <- lapply (coordenadas_especies, function (especie) { # 
    
    coordinates (especie) <- ~ decimalLongitude + decimalLatitude # transformar em spdf
    crs(especie)<- crs(shape_RS) ## aproveitar para definir as projecao de referencia (do shapefile RS)
    ; 
    especie # retornar o spdf de cada sp.

 }
)

# obter municipios com registro de cada especie
mun_det_especies <- lapply (coordenadas_especies,function (especie)
                    
                                over (especie,shape_RS)
                
                            )

# obter um mapa destes municipios com det de cada especie
shape_mun_especies <- lapply (mun_det_especies, function (especie)
        
                                shape_RS [which (shape_RS@data$NM_MUNICIP %in% especie$NM_MUNICIP),]
                              
                              )

#  finalmente, obter um dataframe, com o cod do mun na linha, 
# e a deteccao \ nao det de cada especie na coluna

det_especies <- lapply (shape_mun_especies, function (especie)
        
                    data.frame (cod_mun = shape_RS@data$CD_GEOCMU,
                               det = ifelse (shape_RS@data$CD_GEOCMU %in% especie@data$CD_GEOCMU,1,0),
                                mun = shape_RS@data$NM_MUNICIP)
                    )

## colocar o nome de cada especie no respectivo dataframe com sua deteccao
names(det_especies) <- especies_aves_RS

## dissolver a lista com a deteccao de cada especie, e somar as linhas 
## para ter o numero de especies por municipio
df_todas_sp <- do.call (cbind, 
                        
                        lapply (det_especies, function (especie) 
        
                            especie$det)
                        )

## uso os codigos da primeira tabela da lista = det_especies[[1]]$cod_mun
# tudo esta na mesma ordem
rownames(df_todas_sp) <- det_especies[[1]]$mun

## obter a riqueza por municipio
riqueza_aves_mun_RS <- rowSums (df_todas_sp)

# plot 

f.mun<-fortify(shape_RS, region="NM_MUNICIP") # fortify mapa do RS, comum a todos os mapas
# modelo sem variaveis de ocupacao
cores_SR <- data.frame (cores= riqueza_aves_mun_RS,
                          NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_SR<- cbind (f.mun, 
                    Nespecies= cores_SR [match (f.mun$id, cores_SR$NM_MUNICIP),]$cores)

## inserir estimativas
c_SR <-    ggplot() + geom_polygon(data=f.mun_SR, aes(x=long, y=lat, group=group, 
                                                  color=Nespecies, 
                                             fill=Nespecies), colour = NA, size=1) + 
    labs (title= "Riqueza de aves do RS")+
    scale_fill_gradient2 (low='white', high='darkred', midpoint= 50,na.value = "white",
                          limits=c(0,max(cores_SR$cores)), 
                          breaks=seq(0,max(cores_SR$cores,na.rm=T),by=50),
                          name="Riqueza") ## para continuo
c_SR

########################
## os codigos dos municipios da riqueza fecham  com os cod das deteccoes de ema
# names (riqueza_aves_mun_RS) ==  det_ema$mun
# colar a riqueza às deteccoes
dados_det_ema_gbif <- cbind (det_ema, 
                             riqueza_aves= riqueza_aves_mun_RS)
save (dados_det_ema_gbif, file= here("data",
                                     "organized_data",
                                     "input_GBIF.RData"))
rm(list=ls())
S
