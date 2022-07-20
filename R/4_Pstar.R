# esfor?o
e <- seq (1,100) 

# deteccao
p <- 0.5

# probabilidade de nao detectar em funcao do esforco
par(mfrow=c(2,2),mar=c(4,4,2,4))
plot (sqrt(e),
      (1-p)^sqrt(e), type="l",
      xlab="Esforco amostral",
      ylab = "Probabilidade de deteccao",
      lwd = 3,
      main="(1-p)^E")


# calculo do Pstar (P*) - probabilidade de detecctar em ao menos um sitio em funcao do esfor?o
plot (sqrt(e),
      (1-(1-p)^sqrt(e)), 
      type="l",
      xlab="Esforco amostral",
      ylab = "",
      lwd = 3,
      main="P*=1-(1-p)^E")
