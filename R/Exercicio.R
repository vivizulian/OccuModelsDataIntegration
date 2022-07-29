
#### EXERCICIO ####

## Pesquisadores estavam interessados em estimar a probabilidade de ocupacao de 
## uma especie endemica de anfibio. Para isso, realizaram amostragem acustica 
## replicada em 120 pontos, visitando cada ponto por 4 vezes. 
## Para entender a relacao da especie com o microhabitat, coletaram informacao 
## referente a cada sitio amostrado (umidade). 
## Considerando que a atividade dos anfibios eh afetada pela temperatura, 
## em cada visita, coletaram dados de temperatura media.

# numero de sitios = 120
# numero de visitas = 4
# covariavel de sitio = umidade
# covariavel de amostragem = temperatura

# O objetivo eh modelar a probabilidade de ocupacao da especie em cada sitio, 
# de acordo com as covariaveis em uma abordagem bayesiana. 

# 1. Ler os dados
load(here("data","Exercicio", "dados_exercicio.RData"))

# 2. Crie o objeto para o modelo
str(data <- list(y = observ, nvisit = nvisit, nsitios = nsites, ...))

# 3. Especifique o modelo em linguagem BUGS
sink("exercicio.txt")
cat("
model {
  # Priors
  mean.p ~ dunif(0, 1) 
  alpha0 <- logit(mean.p) # Intercepto para deteccao
  mean.psi ~ dunif(0, 1) 
  beta0 <- logit(mean.psi) # Intercepto para ocupacao
  
  # Likelihood
  # Processo biologico parcialmente observado
  for (i in 1:nsitios) {
    z[i] ~ dbern(psi[i]) # status de cada sitio (ocupado ou nao)
    logit(psi[i]) <- beta0 + 
    
    # Processo amostral para os dados observados
    for (j in 1:nvisit) {
      y[i,j] ~ dbern(z[i] * p[i,j]) # Deteccao nao-deteccao
      logit(p[i,j]) <- alpha0 + 
    }
  }
  
  # Derived quantities
  N.occ <- sum(z[]) # Numero de sitios ocupados
  psi.fs <- N.occ/nsitios # Proporcao de sitios ocupados dentre os sitios amostrados
  
}
",fill = TRUE)
sink()

# 4. De Valores iniciais para o modelo
zst <- apply(observ, 1, max) 

inits <- function(){list(z = zst, mean.p = runif(1), alpha1 = runif(1), mean.psi
                         = runif(1), beta1 = runif(1))}

# 5. Defina os parametros a serem monitorados
params <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs", "p", "z", 
            "psi")

# 6. Defina MCMC 
ni <- 2500 ; nt <- 10 ; nb <- 2000 ; nc <- 3 #teste
#ni <- 25000 ; nt <- 10 ; nb <- 2000 ; nc <- 3

# 7. Chamar o jags
out <- jags(data, inits, params, "exercicio.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb)

# 8. Printe e interprete os resultados
print(out, dig = 3)

