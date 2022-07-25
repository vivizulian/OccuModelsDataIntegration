# ================================================================================= #

# funcao de mapeamento das variaveis
## PS: variavel em raiz quadrada
funcao_mapa_variaveis <- function (shape, variavel, titulo, cor.min, cor.max,cor.na) {
  
  # formato longo (adequado para ggplot)
  cores_hab <- data.frame (cores= sqrt(variavel),
                           NM_MUNICIP=shape$NM_MUNICIP)
  
  # fortify
  f.mun<-fortify(shape, region="NM_MUNICIP")
  f.mun_hab<- cbind (f.mun, 
                     Nespecies = cores_hab [match (f.mun$id, 
                                                   cores_hab$NM_MUNICIP),]$cores)
  
  # grassland cover
  # inserir vals
  c_hab <-   ggplot() + geom_polygon(data=f.mun_hab, aes(x=long, y=lat, group=group, 
                                                         color=Nespecies, 
                                                         fill=Nespecies), 
                                     colour = "gray", size=0.1) + 
    labs (title= titulo) +
    scale_fill_gradient2 (low=cor.min, high=cor.max, na.value = cor.na,
                          limits=c(0,max(cores_hab$cores)), 
                          breaks=seq(0,max(cores_hab$cores,na.rm=T),by=500),
                          name=expression(sqrt("Area (km2)"))) ## para continuo
  
  (c_hab <- c_hab + theme_classic() + 
      theme (axis.text = element_text(size=6),
             axis.title = element_text(size=8),
             legend.text = element_text(size=8),
             legend.title = element_text(size=9))+
      xlab("Longitude") + 
      ylab("Latitude")) 
  # retornar
  return (c_hab)
}




## funcao para capturar legenda
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}





# ================================================================================= #



# Mapeamento das observacoes



# funcao do mapa das observacoes
funcao_mapa_observacoes<- function (shape, observacoes, titulo) {
  
  ## mapas das deteccoes
  cores_data <- data.frame (cores= observacoes,
                            NM_MUNICIP=shape$NM_MUNICIP)
  # escala discreta
  cores_data$cores <- factor(cores_data$cores)
  levels(cores_data$cores) [which(levels (cores_data$cores) == 0)] <- "Not detected"
  levels(cores_data$cores) [which(levels (cores_data$cores) == 1)] <- "Detected"
  levels(cores_data$cores) [which(levels (cores_data$cores) == -Inf)] <- "Not sampled"
  levels(cores_data$cores) [which(levels (cores_data$cores) == "NA")] <- "Not sampled"
  
  # fortify RS map
  f.mun<-fortify(shape, region="NM_MUNICIP") # fortify mapa do RS, comum a todos os mapas
  # bind
  f.mun_data<- cbind (f.mun, 
                      Nespecies= cores_data [match (f.mun$id, 
                                                    cores_data$NM_MUNICIP),]$cores)
  
  ## inserir deteccoes do eBird
  c_data <-   b + geom_polygon(data=f.mun_data, aes(x=long, 
                                                    y=lat, 
                                                    group=group, 
                                                    color=Nespecies, 
                                                    fill=Nespecies), 
                               colour=NA,size=1) + 
    labs (title= titulo) +
    scale_fill_manual("Observation data",
                      values = c("Not detected" = "white",
                                 "Detected" = "darkred",
                                 "Not sampled" = "gray85"))
  
  # anotar
  e_data <- c_data + ggsn::scalebar(f.mun_data, 
                                    dist = 100, 
                                    st.dist=0.03,
                                    st.size=2.2, 
                                    height=0.02, 
                                    transform = TRUE, 
                                    dist_unit = "km",
                                    model = 'WGS84', 
                                    location = "bottomright")
  
  ## plot para extrair a legenda
  f_data_legend <- e_data +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightcyan", 
                                          colour = "lightcyan", 
                                          size = 0.5, 
                                          linetype = "solid"),
          
          legend.title = element_text(size=7),
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
          plot.title = element_text(size=10),
          plot.margin = unit(c(0.1, -0.1,-0.1, 0.2), "lines")) 
  
  # return
  return (f_data_legend)
  
}



