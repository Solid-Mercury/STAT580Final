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
  incomp <- add_noise(lrm, SNR)
  
  #remove
  incomp <- set_NA(incomp, obs.gen(incomp, p))
  
  res <- list(incomp <- incomp, true = lrm)
  return(res)
}


