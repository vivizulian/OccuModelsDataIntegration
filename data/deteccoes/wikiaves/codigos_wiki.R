# ------------------------------------------------------------------ #
#              dados do wikiAves 
#              Obtidos por V Zulian
#       Dados ja vem em um formato adequado para integracao de dados
# ------------------------------------------------------------------- #

# carregar os pacotes
source ('R/packagesR.R')

# carregar dados
dados <- read.csv (here("data","deteccoes","wikiaves","DataWiki2008_2018.csv"),
                   h=T,sep=",")

# carregar shapefile RS
shape_RS <- readOGR (dsn=here("data","shape_munRS"), 
                     layer = "43MUE250GC_SIR",
                     encoding="UTF-8", use_iconv=T)
shape_RS <- shape_RS [-c (96,250),] # rm lagos

## BOTAR O COD_IBGE NOS MUNICIPIOS QUE NAO TEM NOS DADOS DO WIKIAVES
# identificados por busca manual
dados [is.na(dados$id),"id"] <- c(4301800, 4305603,4313334)

## # corrigir nomes
MUN_nomes_corretos <- unlist(lapply (dados$id, function (i)
  
      shape_RS$NM_MUNICIP [which (shape_RS$CD_GEOCMU %in% i)])
      
      )

## substituir os nomes antigos pelos corretos (com acentos e zaz)
dados$MUNI <- MUN_nomes_corretos

## colar na tabela dados os municipios sem dados (faltam tabela wiki)

dados <- rbind (dados,
                data.frame (X=NA,
                    id= shape_RS$CD_GEOCMU [which (shape_RS$NM_MUNICIP %in% dados$MUNI ==F)],
                    MUNI=shape_RS$NM_MUNICIP [which (shape_RS$NM_MUNICIP %in% dados$MUNI ==F)],
                    SITE="RS",
                    STATE = "RIO GRANDE DO SUL",
                    NSPECIES = 0,
                    NPIC = 0,
                    NSONG=0,
                    RHAMERICANA=0))

# matching names
dados_wikiaves <- dados [match(shape_RS$NM_MUNICIP,dados$MUNI),]

# plot
plot(shape_RS, col = rgb(dados_wikiaves$RHAMERICANA,0,1,alpha=0.5),
     main = "Wikiaves")


# save data
save (dados_wikiaves,file=here ("data","organized_data","input_wikiaves.RData"))

rm(list=ls())

                     


