library(tidyverse)
library(network)
library(mixer)
library(RColorBrewer)
library(randnet)
library(irlba)

options(OutDec = ",")
load("matriz_adj.RData")
heatmap(matriz, Rowv=NA,Colv=NA,scale="none", cexRow = 1.5, cexCol = 1.5)
rowSums(matriz) %>% hist(main = element_blank(),
                         ylab = "Frequência de departamentos",
                         xlab = "Número de ligações",
                         cex.lab = 1.5,
                         cex.axis = 1.5,
                         xlim = c(0,20))
set.seed(1)
rede <- network(matriz)
layout <- network.layout.fruchtermanreingold(rede,
                                             layout.par=NULL)
par(mfrow = c(1, 1), oma = c(.2,.2,0.2,0.2))
plot(rede, 
     label = rownames(matriz),
     mode = "fruchtermanreingold",
     coord = layout, label.pos = 0, label.pad = -0.7, 
     edge.col = "gray")
rowSums(matriz) %>% hist(main = element_blank(),
                         ylab = "Frequência de departamentos",
                         xlab = "Número de ligações")
set.seed(2)
#https://cran.r-project.org/src/contrib/Archive/mixer/
#install.packages("mixer_1.9.tar.gz",repos=NULL,source=TRUE)
#Candidatos:1a5clusteres
fitb <- mixer(matriz, qmin=1, qmax=10,method= "variational")
criterios <- purrr::map(1:length(fitb$output), function(i) data.frame(q = i,
                                                                      criterio = fitb$output[[i]]$criterion))
criterios <- list_rbind(criterios)
plot(criterios$q, criterios$criterio, ylab = "ICL", xlab = "Número de clusters",
     cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, pch = 16)
criterios %>% arrange(-criterio) %>% head(3)

sbm4 <- mixer(matriz, qmin = 4, qmax = 4,method= "variational")
mod4<-getModel(sbm4)
mod4$Pis %>% apply(1, which.max) # comunidades

sbm5 <- mixer(matriz, qmin = 5, qmax = 5,method= "variational")
mod5<-getModel(sbm5)
mod5$Pis %>% apply(1, which.max)

sbm7 <- mixer(matriz, qmin = 7, qmax = 7,method= "variational")
mod7<-getModel(sbm7)
mod7$Pis %>% apply(1, which.max)

cores <- brewer.pal(3, "Dark2")

z <- t(mod4$Taus)
(1 - apply(z, 1, max)) %>% sort %>% plot(ylim = c(0,.6), col = cores[1],
                                         pch = 15, cex = 1.5,
                                         cex.lab = 1.5,
                                         cex.axis = 1.5,
                                         xlab = "Índice do vértice",
                                         ylab = "Incerteza a posteriori")
z <- t(mod5$Taus)
(1 - apply(z, 1, max)) %>% sort %>% points(ylim = c(0,.6), col = cores[2],
                                           pch = 16, cex = 1.5)
z <- t(mod7$Taus)
(1 - apply(z, 1, max)) %>% sort %>% points(ylim = c(0,.6), col = cores[3],
                                           pch = 18, cex = 2)
legend("topleft", col = cores, pch = c(15, 16, 18), legend = c(4,5,7),
       title = "Número de clusters", cex = 1.5)


sbm <- mixer(matriz, qmin = 7, qmax = 7,method= "variational")

sbm$edges
mod<-getModel(sbm)
mod$Pis %>% apply(1, which.max)
z <- t(mod$Taus)
clusters <- mclust::map(z)

plota_rede_clusterizada <- function(clusters,...){
  
  names(clusters) <- rownames(matriz)
  names(clusters)[clusters == 5]
  color <- brewer.pal(max(clusters), "Accent")
  ccolor <- color[clusters] # Etiquetas de acordo com cluster
  plot(rede,
       label = names(clusters),
       coord = layout, vertex.col = ccolor,
       label.pos = 0, label.pad = -0.7,
       ...)
  legend("topright", col = color, pch = 20, legend = 1:max(clusters),
         title = "Agrupamento", cex = 1.5, pt.cex = 1.5)
}
par(mar = c(0,0,0.2,7))
plota_rede_clusterizada(mclust::map(z), edge.col = "gray")
dev.off()
(1 - apply(z, 1, max)) %>% sort %>% plot()

table(clusters)
mod$Pis %>% round(4) %>% xtable::xtable(digits = 2) 


set.seed(5)
rlang::env_unlock(env = asNamespace('randnet'))
rlang::env_binding_unlock(env = asNamespace('randnet'))
reg.SSP <- function (A, K, tau = 1, lap = FALSE, nstart = 30, iter.max = 100) 
{
  avg.d <- mean(colSums(A))
  A.tau <- A + tau * avg.d/nrow(A)
  if (!lap) {
    SVD <- irlba(A.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  else {
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau)) %*% A.tau %*% diag(1/sqrt(d.tau))
    SVD <- irlba(L.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  km <- kmeans(V.normalized, centers = K, nstart = nstart, 
               iter.max = iter.max)
  return(list(cluster = km$cluster, loss = km$betweenss))
}
assign('reg.SSP', reg.SSP, envir = asNamespace('randnet'))
rlang::env_binding_lock(env = asNamespace('randnet'))
rlang::env_lock(asNamespace('randnet'))
loss_ssp <- purrr::map(2:20, function(K){
  regssp <- randnet::reg.SSP(A = matriz, K = K, tau = 1, lap = T)
  regssp$loss
})
plot(2:20, loss_ssp, xlab = "Número de agrupamentos",
     pch = 16, cex = 1.5,
     ylab = "Soma de quadrados entre agrupamentos", 
     cex.axis = 1.5, cex.lab = 1.5)
regssp <- reg.SSP(A = matriz, K = 7, tau = 1, lap = T)

regssp$cluster
par(mar = c(0,0,0.2,7))
plota_rede_clusterizada(regssp$cluster, edge.col = "gray")
table(regssp$cluster)
clusters2 <- regssp$cluster

ligacoes <- matrix(0, nrow = 7, ncol = 7)
for(i in 1:7){
  for(j in 1:7){
    ligacoes[i,j] <- matriz[clusters2 == i, clusters2 == j] %>% mean
  }
}
ligacoes %>% round(2) %>% xtable::xtable(digits = 2)
