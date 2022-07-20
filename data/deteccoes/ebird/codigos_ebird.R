# --------------------------------------------------------------------------- #
#                     Dados eBird
# Dados ja processados pela Viviane Zulian, para ter  somente dados do RS
# --------------------------------------------------------------------------- #

# carregar os pacotes
source ('R/packagesR.R')

# carregar dados do eBird
ebirdDataRS <- read.csv (here("data","deteccoes","ebird","ebirdDataRS.csv"), h=T,
                         encoding = "UTF-8")

# carregar shapefile dos municipios do RS
shape_RS <- readOGR(dsn=here ("data","shape_munRS"), layer="43MUE250GC_SIR",
                    encoding = "UTF-8",use_iconv = T)
shape_RS <- shape_RS [-c(96,250),]## remover os lagos

## extrair as coordenadas dos registros do eBird (PS: eBird nao trabalha com municipios)
all_coord <-  ebirdDataRS[,c("LONGITUDE", "LATITUDE")]
coordinates (all_coord) <- ~ LONGITUDE + LATITUDE # transformar em spatialpoints DF
crs(all_coord) <- crs(shape_RS) # atribuir um sistema de referencia geografica

## sobrepor as coordenadas no shape do RS
over_coord <- over (all_coord,shape_RS) ## ver a qual municipio cada registro do eBird pertence
ebirdDataRS$NM_MUNICIP <- over_coord$NM_MUNICIP ## adicionar uma coluna de municipio a tabela do eBird

# remover os registros que não caem em nenhum municipio
ebirdDataRS <- ebirdDataRS [which(is.na(ebirdDataRS$NM_MUNICIP)!=T),]

## um DF para padronizar as id dos municipios
ID_MUN <- data.frame (NM_MUNICIP = shape_RS$NM_MUNICIP,
                      ID_MUN = seq (1, length(shape_RS$NM_MUNICIP)))

## ver a ID de cada mun da tabela do ebird
ID_MUN_ebird <- lapply (ebirdDataRS$NM_MUNICIP, function (i)
                                  
                                    ID_MUN$ID_MUN [which (ID_MUN$NM_MUNICIP == i )]
                        )
ID_MUN_ebird <- unlist (ID_MUN_ebird) # derreter a lista
ebirdDataRS$ID_MUN <- ID_MUN_ebird ## colar os IDs na tabela do ebird

# vamos criar uma lista de municipios com dados, para formar uma tabela de deteccoes por municipio
municipios <- unique (ebirdDataRS$NM_MUNICIP) # munucipios com dados

# subconjunto da tabela completa, extraindo somente os descritores que interessam (IDs e esforco) 
ebirdDataRS_mun <- lapply (municipios, function (i)
  ebirdDataRS [which(ebirdDataRS$NM_MUNICIP == i),c("GLOBAL.UNIQUE.IDENTIFIER", # identificador da lista de spp
									"DURATION.MINUTES",  # duracao da amostragem, em min,
                  "EFFORT.DISTANCE.KM",  # distancia percorrida para completar uma lista
									"NUMBER.OBSERVERS", # numero de observadores que produziram a lista
									"R_AMERIC", # deteccao de  Rhea americana
									"NSPECIES", # numero de spp
									"NM_MUNICIP")] # identidade do municipio da lista
  )

## criar uma seq que vai de 1 ate o numero de listas (checklists) que um municipio tem

ebirdDataRS_mun <- lapply (ebirdDataRS_mun, function (i)
 		cbind (i, cont_lista= seq (1,nrow(i))))


################################
####### DURATION.MINUTES #######
################################

## construir a tabela (mun x checklist)
mun_list <- lapply (ebirdDataRS_mun, function (i) 
			
                      cast (NM_MUNICIP ~ cont_lista,  # tabela dinamica do excel
				
                            data= i,value="DURATION.MINUTES",na.rm=T, FUN=sum) # valor que desejo ter nas celulas da tabela
                    
                    )

## inserir dados das listas inexistentes (inserir listas que faltam, inserindo NA)
mun_list <- lapply (mun_list, function (i)
	
  cbind (i,
	
	       matrix (rep(NA,length(seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))), # add tantas colunas quanto as listas inexistentes (em relacao ao maximo de checklists)
		
	               byrow=F,nrow=1, # 1 municipio
		
	               dimnames=list(NULL, 
			
	                             seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))
	               )
	       )
	)

## add municipios sem amostragem
para_add <- data.frame (
  
                matrix (rep (NA, length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]) * ncol(mun_list[[1]])), # quais municipios nao estao nos checklists
		
                        nrow=length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]),
		
                        byrow=T)
                )

para_add [,1] <- shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)] # id dos municipios fora das checklists
colnames (para_add) <- c ("NM_MUNICIP",seq (1,ncol (para_add)-1))

## df com os dados
df_mun_duration <- do.call (rbind,mun_list)# derrete a lista para termos um DF
df_mun_duration  <- rbind (df_mun_duration, para_add) # add municipios faltantes
df_mun_duration <- df_mun_duration [match(shape_RS$NM_MUNICIP, df_mun_duration$NM_MUNICIP),] # ordem (IMPORTANTE)
rownames (df_mun_duration) <- df_mun_duration$NM_MUNICIP
df_mun_duration <- df_mun_duration [,-1] # rm coluna 1 (nome mun)

##################################
####### EFFORT.DISTANCE.KM #######
##################################
# (mesmos procedimentos caracterizados anteriormente)
## construir a tabela 
mun_list <- lapply (ebirdDataRS_mun, function (i) 
			              cast (NM_MUNICIP ~ cont_lista,# tabela dinamica
				                  data= i,value= "EFFORT.DISTANCE.KM",# valor das celulas
			                  	na.rm=T, FUN=sum)
			              )

## inserir dados das listas inexistentes (add listas faltantes)
mun_list <- lapply (mun_list, function (i)
	cbind (i,
	  matrix (rep(NA,length(seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))),
		      byrow=F,nrow=1,
		      dimnames=list(NULL, 
			              seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))
		      )
	  )
	)

## add municipios sem amostragem
para_add <- data.frame (matrix (rep (NA, length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]) * ncol(mun_list[[1]])),
		nrow=length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]),
		byrow=T))
para_add [,1] <- shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)] # id do municipio
colnames (para_add) <- c ("NM_MUNICIP",seq (1,ncol (para_add)-1))

## df com os dados
df_mun_distance <- do.call (rbind,mun_list)
df_mun_distance <- rbind (df_mun_distance, para_add)
df_mun_distance <- df_mun_distance [match(shape_RS$NM_MUNICIP, df_mun_distance$NM_MUNICIP),]
rownames (df_mun_distance) <- df_mun_distance$NM_MUNICIP
df_mun_distance <- df_mun_distance [,-1]

##################################
####### NUMBER.OBSERVERS   #######
##################################

## construir a tabela
mun_list <- lapply (ebirdDataRS_mun, function (i) 
			                cast (NM_MUNICIP ~ cont_lista,
				                    data= i,value= "NUMBER.OBSERVERS",
				                    na.rm=T, FUN=sum)
                    )

## inserir dados das listas inexistentes
mun_list <- lapply (mun_list, function (i)
	                cbind (i,
                    	matrix (rep(NA,length(seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))),
	                          	byrow=F,nrow=1,
                            	dimnames=list(NULL, 
                              			      seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))
	                          	)
                    	)
	)

## add municipios sem amostragem
para_add <- data.frame (matrix (rep (NA, length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]) * ncol(mun_list[[1]])),
		nrow=length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]),
		byrow=T))
para_add [,1] <- shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]
colnames (para_add) <- c ("NM_MUNICIP",seq (1,ncol (para_add)-1))

## df com os dados
df_mun_observers <- do.call (rbind,mun_list)
df_mun_observers <- rbind (df_mun_observers, para_add)
df_mun_observers <- df_mun_observers [match(shape_RS$NM_MUNICIP, df_mun_observers$NM_MUNICIP),]
rownames (df_mun_observers) <- df_mun_observers$NM_MUNICIP
df_mun_observers <- df_mun_observers [,-1]

################################
####### NSPECIES ###############
################################

## construir a tabela
mun_list <- lapply (ebirdDataRS_mun, function (i) 
                        cast (NM_MUNICIP ~ cont_lista,
                              data= i,value="NSPECIES",na.rm=T, FUN=sum)
                    )

## inserir dados das listas inexistentes
mun_list <- lapply (mun_list, function (i)
                        cbind (i,
                                matrix (rep(NA,length(seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))), # add tantas colunas quanto as listas inexistentes
                                        byrow=F,nrow=1,
                                        dimnames=list(NULL, 
                                                      seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))
                                        )
                               )
                )

## add municipios sem amostragem
para_add <- data.frame (matrix (rep (NA, length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]) * ncol(mun_list[[1]])),
                                nrow=length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]),
                                byrow=T)
                        )
para_add [,1] <- shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]
colnames (para_add) <- c ("NM_MUNICIP",seq (1,ncol (para_add)-1))

## df com os dados
df_mun_species <- do.call (rbind,mun_list)
df_mun_species  <- rbind (df_mun_species, para_add)
df_mun_species <- df_mun_species [match(shape_RS$NM_MUNICIP, df_mun_species$NM_MUNICIP),]
rownames (df_mun_species) <- df_mun_species$NM_MUNICIP
df_mun_species <- df_mun_species [,-1]

##################################################
######    deteccoes de RHEA AMERICANA    #########
##################################################

## construir a tabela
mun_list <- lapply (ebirdDataRS_mun, function (i) 
              			cast (NM_MUNICIP ~ cont_lista,
				                  data= i,value= "R_AMERIC",
				                  na.rm=T, FUN=sum))

## inserir dados das listas inexistentes
mun_list <- lapply (mun_list, function (i)
	                        cbind (i,
                                	matrix (rep(NA,length(seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))),
                                      		byrow=F,nrow=1,
                                      		dimnames=list(NULL, 
                                                  			seq (ncol(i), max(unlist(lapply (mun_list,ncol)))))
                                      		)
	                               )
                        )

## add municipios sem amostragem
para_add <- data.frame (matrix (rep (NA, length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]) * ncol(mun_list[[1]])),
		nrow=length(shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]),
		byrow=T))

para_add [,1] <- shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% unlist(lapply (mun_list, function (i) i$NM_MUNICIP ))==F)]
colnames (para_add) <- c ("NM_MUNICIP",seq (1,ncol (para_add)-1))

## df com os dados
df_mun_rhea <- do.call (rbind,mun_list)
df_mun_rhea <- rbind (df_mun_rhea, para_add)
df_mun_rhea <- df_mun_rhea [match(shape_RS$NM_MUNICIP, df_mun_rhea$NM_MUNICIP),]
rownames (df_mun_rhea) <- df_mun_rhea$NM_MUNICIP
df_mun_rhea <- df_mun_rhea [,-1]

# formatar como uma matrix
y.ebird <- unname (as.matrix (df_mun_rhea)) 

## formatar as tabelas de covariaveis de observacao
# x <- 0 : se NA, significa que nao teve amostragem em um dado municipio, portanto == 0
dist.ebird <- as.matrix(df_mun_distance) 
dist.ebird [is.na(dist.ebird)] <-0
dist.ebird <- unname(dist.ebird)
dist.duration <- as.matrix(df_mun_duration)
dist.duration [is.na(dist.duration)] <-0
dist.duration <- unname(dist.duration)
dist.observers <- as.matrix(df_mun_observers)
dist.observers [is.na(dist.observers)] <-0  
dist.observers <- unname(dist.observers)

## plot

plot (shape_RS, main = "eBird", 
      col = rgb (rowSums(df_mun_rhea,na.rm=T)>0,0,1,alpha=0.5),
      border = "black")
#legend ('bottomleft',legend =c("Detection of Rhea",
#                               "No detection of Rhea"),
#        col = c("purple","blue"),pch=15,bty="n")
#

# salvar
save (y.ebird, 
      dist.ebird, 
      dist.duration, 
      dist.observers, 
      ID_MUN,
      file=here ("data","organized_data","input_ebird.RData"))


rm(list=ls())


