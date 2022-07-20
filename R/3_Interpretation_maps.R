
# -----------------------------------------------------------------------
#                   Pacotes e dados
# -----------------------------------------------------------------------

# carregar os pacotes
source ('R/packagesR.R')

# abrir resultados dos modelos
load (here("output", "out_model1_glm.RData")) #     1 - modelos sem integracao, e sem consideracao de esforco amostral (eBird) - GLM
load (here("output", "out_model2_unmarked.RData")) #     2 - modelos sem integracao, e consideracao de esforco amostral (eBird) - HM 'unmarked'
load (here("output", "out_model2_bugs.RData")) #     2 - modelos sem integracao, e consideracao de esforco amostral (eBird) - HM BUGS
load (here("output", "model3_glm.RData")) #     3 - modelos integracao (agregacao), e sem esforco amostral (eBird) - GLM
load (here("output", "model3_bugs.RData")) #     3 - modelos integracao (agregacao), e sem esforco amostral (eBird) - BUGS
load (here("output", "model4_bugs.RData")) #      4 - modelos com integracao, e consideracao de esforco amostral (eBird, GBIF, WikiAves) - BUGS

# abrir dados das covariaveis
load (here("output", "covariaveis.RData")) #      4 - modelos com integracao, e consideracao de esforco amostral (eBird, GBIF, WikiAves) - BUGS

# -------------------------------------------------- #
#         ANALISE DOS COEFICIENTES
# modelo1
m1<-cbind (coef(model1),
       confint (model1),
       model = "NoIntegration,NoDetection")
# model 2 unmarked - coeficientes de ocupacao
m2<-cbind( coef(model2,"state"),
       confint(model2, type='state', method = 'normal'),
       model = "NoIntegration,Detection"
)
# modelo 2 BUGS
m2B<-cbind(
  out_model2$summary[grep("BETA",rownames(out_model2$summary)),c("mean", "2.5%", "97.5%")],
  model="NoIntegration,DetectionBUGS")
# model 3 
m3<-cbind (coef(model3),
       confint (model3),
       model = "Integration(Aggregation),NoDetection")
# m3 bugs
m3B<-cbind(
  model3_bugs$summary[grep("BETA",rownames(model3_bugs$summary)),c("mean", "2.5%", "97.5%")],
  model="Integration(Aggregation),NoDetectionBUGS")
# modelo4
m4 <- cbind(
  model4$summary[grep("BETA",rownames(model4$summary)),c("mean", "2.5%", "97.5%")],
  model="Integration,DetectionBUGS")

# colocar o mesmo nome de coluna em todos
colnames(m1)<- colnames(m2)<-colnames(m2B)<-colnames(m3)<-colnames(m3B)<-colnames(m4)

# rbind 
df_coeficientes <- rbind(m1,m2,m2B,m3,m3B,m4)
df_coeficientes <- data.frame(df_coeficientes,
                         coeficientes=rownames(df_coeficientes))
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="campo")]<- "CAMPO"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="BETA.CAMPO")]<- "CAMPO"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="agri")]<- "AGRI"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="BETA.AGRI")]<- "AGRI"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="psi(agri)")]<- "AGRI"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="psi(campo)")]<- "CAMPO"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="psi(Int)")]<- "INTERCEPTO"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="(Intercept)")]<- "INTERCEPTO"
df_coeficientes$coeficientes [which(df_coeficientes$coeficientes=="BETA0")]<- "INTERCEPTO"
# ajustar nomes das colunas
colnames(df_coeficientes) <- c("mean","lower","upper","model","Coeficientes")
# transformar em numero
df_coeficientes$mean <- as.numeric(df_coeficientes$mean)
df_coeficientes$lower <- as.numeric(df_coeficientes$lower)
df_coeficientes$upper <- as.numeric(df_coeficientes$upper)

# ordenar modelos para aparecere corretamente no grafico
df_coeficientes$model<-factor(df_coeficientes$model,
                              levels= c("NoIntegration,NoDetection",
                                        "Integration(Aggregation),NoDetection",
                                        "Integration(Aggregation),NoDetectionBUGS",
                                        "NoIntegration,Detection",
                                        "NoIntegration,DetectionBUGS",
                                        "Integration,DetectionBUGS"))

# plotar para ver os coeficientes
pd=position_dodge(0.5) # jitter coeff
p1 <- ggplot (df_coeficientes,
              aes (y=as.numeric(mean), x=model, fill=Coeficientes,
                   group = NULL,
                   colour=Coeficientes))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.5,position=pd,size=1) +
  geom_line(position=pd) +
  geom_point(position=pd, aes(colour=Coeficientes),size=2.5)  

p1 <-p1 + scale_color_manual(values=c("#FB9300", "#206A5D","#81B214"))

r <- p1 + scale_y_continuous(limits=c(-4,3),breaks=-2:5.5) + 
  geom_hline(yintercept=0, color="gray10", size=1,alpha=0.4)+
  #coord_cartesian(ylim=c(-1.1,1.1)) +
  theme(plot.margin=unit(c(0,0.5,0,0),"lines"))
r

s_coef <-r +  theme_classic()  + #coord_flip(ylim=c(-1,5))+
  scale_x_discrete()+
  theme(axis.text.x = element_text(angle=65,
                                   hjust = 1,
                                   size=8),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=10),
        legend.position="top") +
  xlab("Modelo") + ylab ("Tamanho de efeito padronizado")
s_coef

# ----------------------------------------------------------
#           Dados para mapeamento
# ----------------------------------------------------------

# -----------------------------------------------
#                modelo 1
# predizer ocorrencia para todos os municipios
predictions_m1 <-plogis( predict.glm(model1,
                              newdata = data.frame(campo=habitat_campo,
                                                   agri=habitat_lavoura)
))
range(
  (predictions_m1)
)

# -----------------------------------------------
#                modelo 2

# get predictions for the whole dataset
predictions_m2 <- predict(model2,type="state",
                          newdata = data.frame (campo = habitat_campo,
                                                agri=habitat_lavoura))
range(predictions_m2$Predicted)

# -----------------------------------------------
#           modelo 2 bugs

predictions_m2B <- out_model2$mean$z

# -----------------------------------------------
#           modelo 3 


predictions_m3 <- plogis(predict.glm(model3,
                              newdata = data.frame(campo=habitat_campo,
                                                   agri=habitat_lavoura)
))

range((predictions_m3))

# -------------------------------
#         model 3 bugs
predictions_m3B <- model3_bugs$mean$psi

# -------------------------------
#         model 4 bugs
predictions_m4B <- model4$mean$z


# histograma comparando predicoes
hist((predictions_m1),xlim=c(0,1),
     col=rgb(0,1,1,alpha=0.4),main="Histograma das predicoes")
hist(predictions_m2$Predicted,add=T,
     col=rgb(1,1,0,alpha=0.4))
hist(predictions_m2B,add=T,
     col=rgb(1,0.5,0,alpha=0.4))
hist(predictions_m3,add=T,
     col=rgb(0.5,0.5,0,alpha=0.4))
hist(predictions_m3B,add=T,
     col=rgb(0.1,0.5,0,alpha=0.4))
hist(predictions_m4B,add=T,
     col=rgb(1,0.5,1,alpha=0.4))

legend ("topright", c("m1","m2","m2B","m3","m3B", "m4B"),
        pch=19,
        col=c(rgb(0,1,1,alpha=0.4),
              rgb(1,1,0,alpha=0.4),
              rgb(1,0.5,0,alpha=0.4),
              rgb(0.5,0.5,0,alpha=0.4),
              rgb(0.1,0.5,0,alpha=0.4),
              rgb(1,0.5,1,alpha=0.4)
              
              ))

# ------------------------------
# deteccao
# unmarked
cbind(
  coef(model2,"det"),
  confint(model2, type='det', method = 'normal')
)

# bugs
out_model2$summary[grep("ALPHA",rownames(out_model2$summary)),]

# N municipios ocupados
sum(predictions_m1)
sum(predictions_m2)
(out_model2$mean$fs.z)*497
sum(predictions_m3)
sum(model3_bugs$mean$psi)
sum(model4$mean$fs.z)*497

# ----------------
# Mapeamento

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
b
# ----------------------------------------
# modelo 1
cores_m1 <- data.frame (cores= predictions_m1,
                          NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m1<- cbind (f.mun, 
                    Nespecies= cores_m1 [match (f.mun$id, cores_m1$NM_MUNICIP),]$cores)

## inserir estimativas
c_m1 <-   b + geom_polygon(data=f.mun_m1, aes(x=long, y=lat, group=group, 
                                                  color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 1")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m1$cores)), 
                        breaks=seq(0,max(cores_m1$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

c_m1 <- c_m1 + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))

# ------------
# modelo 2

cores_m2 <- data.frame (cores= predictions_m2$Predicted,
                        NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m2<- cbind (f.mun, 
                  Nespecies= cores_m2 [match (f.mun$id, cores_m2$NM_MUNICIP),]$cores)

## inserir estimativas
c_m2 <-   b + geom_polygon(data=f.mun_m2, aes(x=long, y=lat, group=group, 
                                              color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 2")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m2$cores)), 
                        breaks=seq(0,max(cores_m2$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

c_m2 <- c_m2 + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))


# --------------------------
# modelo2 bugs

cores_m2B <- data.frame (cores= predictions_m2B,
                        NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m2B<- cbind (f.mun, 
                  Nespecies= cores_m2B [match (f.mun$id, cores_m2B$NM_MUNICIP),]$cores)

## inserir estimativas
c_m2B <-   b + geom_polygon(data=f.mun_m2B, aes(x=long, y=lat, group=group, 
                                              color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 2B")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m2B$cores)), 
                        breaks=seq(0,max(cores_m2B$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

c_m2B <- c_m2B + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))

# -------------------
# modelo 3

cores_m3 <- data.frame (cores= predictions_m3,
                         NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m3<- cbind (f.mun, 
                   Nespecies= cores_m3 [match (f.mun$id, cores_m3$NM_MUNICIP),]$cores)

## inserir estimativas
c_m3 <-   b + geom_polygon(data=f.mun_m3, aes(x=long, y=lat, group=group, 
                                                color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 3")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m3$cores)), 
                        breaks=seq(0,max(cores_m3$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

c_m3 <- c_m3 + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))


# -------------------
# modelo 3 B

cores_m3B <- data.frame (cores= predictions_m3B,
                        NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m3B<- cbind (f.mun, 
                  Nespecies= cores_m3B [match (f.mun$id, cores_m3B$NM_MUNICIP),]$cores)

## inserir estimativas
c_m3B <-   b + geom_polygon(data=f.mun_m3B, aes(x=long, y=lat, group=group, 
                                              color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 3B")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m3B$cores)), 
                        breaks=seq(0,max(cores_m3B$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

c_m3B <- c_m3B + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))

# --------------------------
# MODELO 4
cores_m4B <- data.frame (cores= predictions_m4B,
                         NM_MUNICIP=shape_RS$NM_MUNICIP)
f.mun_m4B<- cbind (f.mun, 
                   Nespecies= cores_m4B [match (f.mun$id, cores_m4B$NM_MUNICIP),]$cores)

## inserir estimativas
c_m4B <-   b + geom_polygon(data=f.mun_m4B, aes(x=long, y=lat, group=group, 
                                                color=Nespecies, fill=Nespecies), colour = NA, size=1) + 
  labs (title= "Modelo 4B")+
  scale_fill_gradient2 (low='white', high='darkred', midpoint= 0.20,na.value = "white",
                        limits=c(0,max(cores_m4B$cores)), 
                        breaks=seq(0,max(cores_m4B$cores,na.rm=T),by=0.2),
                        name="Probabilidade") ## para continuo

## plot para retirar a legenda

c_m4B_legend <- c_m4B +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightcyan", 
                                        colour = "lightcyan", 
                                        size = 0.5, 
                                        linetype = "solid"),
        
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
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
        plot.margin = unit(c(-0.1, -0.1,-0.1, -0.2), "lines")) 
#        legend.margin = margin (0,0,0,0),
#       legend.box.margin = margin(1,0,0,0)) 

c_m4B_legend

## funcao para capturar legenda
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legenda_comum <- get_legend(c_m4B_legend) ## tem que primeiro gerar um mapa com legenda

## plot para o painel (sem a legenda desta vez)  
f_m4B <- c_m4B + 
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
        plot.margin = unit(c(0.2,-2, -3.5, -1.0), "lines"))

# top, right, bottom, and left margins
f_m4B <- f_m4B + annotate(geom="text", x=-56, y=-32, label="URUGUAY",color="black",size=2.1) +
  annotate(geom="text", x=-56.5, y=-27.8, label="ARGENTINA",color="black",size=2.1)+
  annotate(geom="text", x=-56.5, y=-27, label="PARAGUAY",color="black",size=2.1)+
  annotate(geom="text", x=-51.8, y=-26.8, label="Santa Catarina",color="black",size=2.1)+
  annotate(geom="text", x=-50.5, y=-32.5, label="ATLANTIC OCEAN",color="black",size=2.1)
# escala

g_m4B <- f_m4B + ggsn::north(f.mun_m4B, symbol=1,scale = 0.2,location = "bottomleft") +
  theme(legend.text=element_text(size=7),
        legend.title=element_text(size=8))


## arranjo dos mapas

pdf(file=here ("output","maps_estimates_models.pdf"),width = 7,height = 5,family="serif")

grid.arrange(c_m1,
             c_m2,
             c_m2B,
             c_m3,
             c_m3B,
             g_m4B,
             legenda_comum,
             ncol=3,nrow=7,
             layout_matrix = t (data.frame(rep(7,3),
                                    seq (1,3),
                                    seq(1,3),
                                    seq(1,3),
                                    seq(4,6),
                                    seq(4,6),
                                    seq(4,6))))


dev.off()

