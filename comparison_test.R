source("fpfunctions2.R")

time_save <- c()
gen.matrix <- list()
uv <- list()
z <- list()
omega <- list()
incomp <- list()
teste_save1 <- rep(0,18)
time_save1 <- rep(0,18)
teste_save2 <- rep(0,18)
time_save2 <- rep(0,18)
teste_save3 <- rep(0,18)
time_save3 <- rep(0,18)
teste_save4 <- rep(0,18)
time_save4 <- rep(0,18)
teste_temp <- rep(0,50)
lambda <- 500:1


for(k in 1:10){
i = 1
#######first generate 18 matrices
for(p in c(0.3,0.5,0.8)){
  for(SNR in c(1,10)){
    for(truerank in c(5,10,20)){
      gen.matrix[[i]] <- incomp.sim(100,100,truerank,SNR,p)
      uv[[i]] <- (gen.matrix[[i]])$true
      z[[i]] <- (gen.matrix[[i]])$true_wth_noise
      omega[[i]] <- (gen.matrix[[i]])$omega
      incomp[[i]] <- (gen.matrix[[i]])$incomp
##########################This is method 1 #######################
  time_start <- proc.time()
  sim.matrix <- mysoftimpute1(incomp[[i]],lambda,100)
  time_end <- proc.time()
  for (j in 1:length(lambda)){
    teste_temp[j] <- testerror(uv[[i]],sim.matrix[[j]], omega[[i]])
  }
  teste_save1[i] <- teste_save1[i] + min(teste_temp)
  time_save1[i] <- time_save1[i] + time_end[1] - time_start[1]
  ##########################This is method 2 #######################
  time_start <- proc.time()
  sim.matrix <- mysoftimpute2(incomp[[i]],lambda,50)
  time_end <- proc.time()
  for (j in 1:length(lambda)){
    teste_temp[j] <- testerror(uv[[i]],sim.matrix[[j]], omega[[i]])
  }
  teste_save2[i] <- teste_save2[i] + min(teste_temp)  
  time_save2[i] <- time_save2[i] + time_end[1] - time_start[1]
  ##########################This is method 3 #######################
  time_start <- proc.time()
  sim.matrix <- mysoftimpute3(incomp[[i]],lambda,30)
  time_end <- proc.time()
  for (j in 1:length(lambda)){
    teste_temp[j] <- testerror(uv[[i]],sim.matrix[[j]], omega[[i]])
  }
  teste_save3[i] <- teste_save3[i] + min(teste_temp)  
  time_save3[i] <- teste_save3[i] + time_end[1] - time_start[1]
  ##########################This is method 4 #######################
  time_start <- proc.time()
  sim.matrix <- mysoftimpute4(incomp[[i]],lambda)
  time_end <- proc.time()
  for (j in 1:length(lambda)){
    teste_temp[j] <- testerror(uv[[i]],sim.matrix[[j]], omega[[i]])
  }
  teste_save4[i] <- teste_save4[i] + min(teste_temp)  
  time_save4[i] <- time_save4[i] + time_end[1] - time_start[1]
  i <- i+1
    }
  }
}
}
rank_save <- rep(c(5,10,20),6)
SNR_save <- rep(c(1,1,1,10,10,10),3)
p_save <- rep(c(0.3,0.5,0.8),each = 6)
df2 <- data.frame(time1 = time_save1/10,time2 = time_save2/10,time3 = time_save3/10,time4 = time_save4/10,error1 = teste_save1/10,error2 = teste_save2/10,
                  error3 = teste_save3/10,error4 = teste_save4/10,truerank = rank_save, SNR = SNR_save, p = p_save)
save(df2,file = "comparison_testpart.Rda")



load("comparison_testpart.Rda")
df2
