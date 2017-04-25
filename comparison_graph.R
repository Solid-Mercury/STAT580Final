source("fpfunctions.R")

lambda = 100:1
teste <- rep(0,20000)
traine <- rep(0,20000)
rankzhat <- rep(0,20000)
method <- rep(c(1,2,3,4),each = 5000)
p <- 0.5
truerank <- 10
SNR <- 1
#for(p in c(0.3,0.5,0.8)){
#  for(SNR in c(1,10)){
#    for(truerank in c(5,10,20)){
##########################This is method 1 #######################
      for(i in 1:50){
        gen.matrix <- incomp.sim(100,100,truerank,SNR,p)
        uv <- (gen.matrix)$true
        z <- (gen.matrix)$true_wth_noise
        omega <- (gen.matrix)$omega
        incomp <- (gen.matrix)$incomp
        sim.matrix <- mysoftimpute1(incomp,lambda,100)
        for (j in 1:length(lambda)){
          teste[i*100 - 100 + j] <- testerror(uv,sim.matrix[[j]],omega)
          traine[i*100 - 100 + j] <- trainingerror(z,sim.matrix[[j]], omega)
          rankzhat[i*100 - 100 + j] <- qr(sim.matrix[[j]])$rank
        }
      }

################This is method 2- svd package#######################

for(i in 1:50){
  gen.matrix <- incomp.sim(100,100,truerank,SNR,p)
  uv <- (gen.matrix)$true
  z <- (gen.matrix)$true_wth_noise
  omega <- (gen.matrix)$omega
  incomp <- (gen.matrix)$incomp
  sim.matrix <- mysoftimpute2(incomp,lambda,truerank)
  for (j in 1:length(lambda)){
    teste[i*100 + 100*49 + j] <-  testerror(uv,sim.matrix[[j]],omega)
    traine[i*100 + 100*49 + j] <- trainingerror(z,sim.matrix[[j]], omega)
    rankzhat[i*100 + 100*49 + j] <-  qr(sim.matrix[[j]])$rank
  }
}
################This is method 3- irlba#######################
for(i in 1:50){
  gen.matrix <- incomp.sim(100,100,truerank,SNR,p)
  uv <- (gen.matrix)$true
  z <- (gen.matrix)$true_wth_noise
  omega <- (gen.matrix)$omega
  incomp <- (gen.matrix)$incomp
  sim.matrix <- mysoftimpute3(incomp,lambda,30)
  for (j in 1:length(lambda)){
    teste[i*100 + 100*99 + j] <- testerror(uv,sim.matrix[[j]],omega)
    traine[i*100 + 100*99 + j] <- trainingerror(z,sim.matrix[[j]], omega)
    rankzhat[i*100 + 100*99 + j] <- qr(sim.matrix[[j]])$rank
  }
}
################This is method 4- RCPP#######################
for(i in 1:50){
  gen.matrix <- incomp.sim(100,100,truerank,SNR,p)
  uv <- (gen.matrix)$true
  z <- (gen.matrix)$true_wth_noise
  omega <- (gen.matrix)$omega
  incomp <- (gen.matrix)$incomp
  sim.matrix <- mysoftimpute4(incomp,lambda)
  for (j in 1:length(lambda)){
    teste[i*100 + 100*149 + j] <- testerror(uv,sim.matrix[[j]],omega)
    traine[i*100 + 100*149 + j] <- trainingerror(z,sim.matrix[[j]], omega)
    rankzhat[i*100 + 100*149 + j] <- qr(sim.matrix[[j]])$rank
  }
}
df <- data.frame(testerror = teste, trainingerror = traine, rank = rankzhat, method = method)
save(df, file = "graphpart.Rda" )
#    }
#  }
#}



load("graphpart.Rda")
require(ggplot2)
df$method[df$method == 1] <- "original svd"
df$method[df$method == 2] <- "svd package"
df$method[df$method == 3] <- "irlba package"
df$method[df$method == 4] <- "rcpparmadillo"
df$method <- factor(df$method)
ggplot(df, aes(x = rank,y = testerror, color = method)) + geom_smooth() + theme_bw()





