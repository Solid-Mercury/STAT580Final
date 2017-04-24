########################################################################
#########             functions to do softimpute             ###########
#Do SVD and soft-threshold using method1
mysoft1 <- function(x,zold,lambda,rank.max){
  res = matrix(rep(0, nrow(x)*ncol(x)), nrow = nrow(x))
  res[!is.na(x)] = res[!is.na(x)] + x[!is.na(x)]
  res[is.na(x)]= res[is.na(x)] + zold[is.na(x)]
  s <- svd(res,nv = min(rank.max,qr(res)$rank), nu = min(rank.max,qr(res)$rank))
  D <- s$d[1:min(rank.max,qr(res)$rank)]
  D <- D - lambda
  D[D < 0] <- 0
  U <- (s$u)
  V <- (s$v)
  D <- diag(D)
  res <- (U%*%D)%*%(t(V))
  return(res)
}

#Do SVD and soft-threshold using method2
require(svd)
mysoft2 <- function(x,zold,lambda,truerank){
  res = matrix(rep(0, nrow(x)*ncol(x)), nrow = nrow(x))
  res[!is.na(x)] = res[!is.na(x)] + x[!is.na(x)]
  res[is.na(x)]= res[is.na(x)] + zold[is.na(x)]
  s <- propack.svd(res,neig = truerank)
  D <- s$d
  D <- D - lambda
  D[D < 0] <- 0
  U <- (s$u)
  V <- (s$v)
  D <- diag(D)
  res <- (U%*%D)%*%(t(V))
  return(res)
}

#Do SVD and soft-threshold using method3
require(irlba)
mysoft3 <- function(x,zold,lambda,rank.max){
  res = matrix(rep(0, nrow(x)*ncol(x)), nrow = nrow(x))
  res[!is.na(x)] = res[!is.na(x)] + x[!is.na(x)]
  res[is.na(x)]= res[is.na(x)] + zold[is.na(x)]
  s <- irlba(res,nv = rank.max)
  D <- s$d[1:min(rank.max,qr(res)$rank)]
  D <- D - lambda
  D[D < 0] <- 0
  U <- (s$u)
  V <- (s$v)
  D <- diag(D)
  res <- (U%*%D)%*%(t(V))
  return(res)
}
#Do SVD and soft-threshold using method4 rcpp

library(Rcpp)
library(RcppArmadillo)
sourceCpp("armadillo_svd.cpp")
mysoft4 <- function(x,zold,lambda){
  res = matrix(rep(0, nrow(x)*ncol(x)), nrow = nrow(x))
  res[!is.na(x)] = res[!is.na(x)] + x[!is.na(x)]
  res[is.na(x)]= res[is.na(x)] + zold[is.na(x)]
  s <- armadillo_svd(res)
  D <- s$s
  D <- D - lambda
  D[D < 0] <- 0
  U <- (s$U)
  V <- (s$V)
  D <- diag(c(D))
  res <- (U%*%D)%*%(t(V))
  return(res)
}


#the soft impute algorithm
mysoftimpute1 <- function (x,lambda,rank.max){
  zhat <- list()
  zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
  znew = mysoft1(x,zold,lambda[1],rank.max)
  zold = znew
  for (i in 1:length(lambda)){
    zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
    while(1){
      znew = mysoft1(x,zold,lambda[i],rank.max)
      diff=(sum((znew-zold)^2))/(sum(zold^2))
      if (is.na(diff)){
        break
      }
      if(diff < 10^-5){
        break
      }
      zold <- znew
    }
    zhat[[i]] <- znew
  }
  return (zhat)
}
#method 2
mysoftimpute2 <- function (x,lambda,truerank){
  zhat <- list()
  zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
  znew = mysoft2(x,zold,lambda[1],truerank)
  zold = znew
  for (i in 1:length(lambda)){
    zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
    while(1){
      znew = mysoft2(x,zold,lambda[i],truerank)
      diff=(sum((znew-zold)^2))/(sum(zold^2))
      if (is.na(diff)){
        break
      }
      if(diff < 10^-5){
        break
      }
      zold <- znew
    }
    zhat[[i]] <- znew
  }
  return (zhat)
}
#method 3
mysoftimpute3 <- function (x,lambda,rank.max){
  zhat <- list()
  zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
  znew = mysoft3(x,zold,lambda[1],rank.max)
  zold = znew
  for (i in 1:length(lambda)){
    zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
    while(1){
      znew = mysoft3(x,zold,lambda[i],rank.max)
      diff=(sum((znew-zold)^2))/(sum(zold^2))
      if (is.na(diff)){
        break
      }
      if(diff < 10^-5){
        break
      }
      zold <- znew
    }
    zhat[[i]] <- znew
  }
  return (zhat)
}
#method 4 rcpp
mysoftimpute4 <- function (x,lambda){
  zhat <- list()
  zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
  znew = mysoft4(x,zold,lambda[1])
  zold = znew
  for (i in 1:length(lambda)){
    zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
    while(1){
      znew = mysoft4(x,zold,lambda[i])
      diff=(sum((znew-zold)^2))/(sum(zold^2))
      if (is.na(diff)){
        break
      }
      if(diff < 10^-5){
        break
      }
      zold <- znew
    }
    zhat[[i]] <- znew
  }
  return (zhat)
}









##################################################
#####     functions to generate matrix     #######

# function to generate low rank matrix. m and n are number of rows and columns of the generated matrix. Rank is the desired rank of generated matrix.
lrm.gen <- function(m, n, rank){
  U <- matrix(rnorm(m*rank), nrow = m)
  Vt <- matrix(rnorm(rank*n), nrow = rank)
  return(U%*%Vt)
}

# function to add noise to a matrix. X is the true matrix and SNR is the desired signal-to-noise ratio. Large SNR preferred for imputing.
add_noise <- function(X, SNR){
  noise <- rnorm(length(X))
  
  # scale noise to satisfy SNR
  noise <- sqrt(var(X[1:length(X)])/var(noise))/SNR * noise
  
  # shape the noise into the shape of true matrix
  dim(noise) <- dim(X)
  return(X + noise)
}

# function to generate p*100% missing postions in matrix X. Input p is in (0,1). Return a matrix of 0 and 1. (observed is 1).
obs.gen <- function(X, p){
  obs_ind <- sample(length(X), p*length(X))
  omega <- array(1, dim = dim(X))
  omega[obs_ind] <- 0
  return(omega)
}

# function to set unobserved postions NA
set_NA <- function(X, omega){
  X[omega == 0] <- NA
  return(X)
}

# generate low rank matrix, add noise and randomly remove
incomp.sim <- function(m,n,rank,SNR,p){
  # generate
  lrm <- lrm.gen(m,n,rank)
  
  # add noise
  lrm_with_noise <- add_noise(lrm, SNR)
  
  # omega
  omega <- obs.gen(lrm, p)
  
  #remove
  incomp <- set_NA(lrm_with_noise, omega)
  
  res <- list(incomp = incomp, true = lrm, true_wth_noise = lrm_with_noise, omega = omega)
  return(res)
}


##########################################################
#####        functions to do comparison      #############
#function to find test error
testerror <- function(uv,zhat,omega){
  resnum <- (norm((uv-zhat)*(1-omega), type = "F"))^2
  resden <- (norm(uv*(1-omega), type = "F"))^2
  res <- resnum/resden
  return(res)
}
#function to find training error
trainingerror <- function(z,zhat,omega){
  resnum <- (norm((z-zhat)*omega, type = "F"))^2
  resden <- (norm(z*omega, type = "F"))^2
  res <- resnum/resden
  return(res)
}









