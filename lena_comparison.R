##########################################################
#####        functions to do comparison      #############
##########################################################
#####        functions to do comparison      #############
#function to find test error. Compare over the rest 40% of image, which is not taken as data.
testerror <- function(uv,zhat,omega){
  resnum <- (norm((uv-zhat)*(1-omega), type = "F"))^2
  resden <- (norm(uv*(1-omega), type = "F"))^2
  res <- resnum/resden
  return(res)
}
#function to find training error. Compare over validation set.
trainingerror <- function(z,zhat,omega){
  resnum <- (norm((z-zhat)*omega, type = "F"))^2
  resden <- (norm(z*omega, type = "F"))^2
  res <- resnum/resden
  return(res)
}

############################################################
############ original svd in R##############################
#a) R interval svd function
z_new_r_svd<-function(z,z_old,lambda){
  r_zero<-row(z)[which(z!=0)]
  c_zero<-col(z)[which(z!=0)]
  z_old[cbind(r_zero,c_zero)]<-0
  p<-z+z_old
  s<-svd(p)
  D<-diag(s$d)
  U<-(s$u)
  V<-(s$v)
  D<-D-lambda*diag(length(s$d))
  D[D<0]<-0
  z_new<-U%*%D%*%t(V)
  return (z_new)
}
r_svd<-function(z,lambda){
  z_hat<-array(rep(0,nrow(z)*ncol(z)*length(lambda)),dim=c(nrow(z),ncol(z),length(lambda)))
  z_old<-matrix(rep(0,nrow(z)*ncol(z)),nrow=nrow(z))
  for (i in 1:length(lambda)){
    z_new<-z_new_r_svd(z,z_old,lambda[i])
    diff=1
    while(diff>10^-5){
      z_old<-z_new
      z_new<-z_new_r_svd(z,z_old,lambda[i])
      diff<-sum((z_new-z_old)^2)/(sum(z_old^2))	}
    z_hat[,,i]<-z_new}
  return (z_hat)
}




############################################################
################ read in data###############################
lena_test<-read.table("~/Documents/isu/STAT 580/lena/test.txt")
lena_training<-read.table("~/Documents/isu/STAT 580/lena/training.txt")
lena_validation<-read.table("~/Documents/isu/STAT 580/lena/validating.txt")
lena_test_omega<-read.table("~/Documents/isu/STAT 580/lena/test_omega.txt")
lena_validation_omega<-read.table("~/Documents/isu/STAT 580/lena/validating_omega.txt")
lena_training_omega<-read.table("~/Documents/isu/STAT 580/lena/trainng_omega.txt")
lena<-matrix(scan("~/Documents/isu/STAT 580/lena/lena256",skip=1),nrow=256)

lena_test[is.na(lena_test)]<-0

lena_training[is.na(lena_training)]<-0

lena_validation[is.na(lena_validation)]<-0
lena_validation_omega[is.na(lena_validation_omega)]<-0
lena_training<-as.matrix(lena_training)
lena_test<-as.matrix(lena_test)
lena_validation<-as.matrix(lena_validation)
lena_validation_omega<-as.matrix(lena_validation_omega)
lena_test_omega<-as.matrix(lena_test_omega)
lena_training_omega<-as.matrix(lena_training_omega)
lena<-as.matrix(lena)

lambda=seq(500,10,-10)

lena_training_hat<-array(rep(0,nrow(lena_training)*ncol(lena_training)*length(lambda)),dim=c(nrow(lena_training),ncol(lena_training),length(lambda)))
lena_training_error<-rep(0,length(lambda))

for (i in 1:length(lambda)){
  lena_training_hat[,,i]<-r_svd(lena_training,lambda[i])
  lena_training_error[i]<-trainingerror(lena_validation,lena_training_hat[,,i], lena_validation_omega)
  
}
best_lambda<-lambda[which.min(lena_training_error)]
best_lambda
#190
min(lena_training_error)
#0.04585135
lena_test_hat<-r_svd(lena_test,best_lambda)[,,1]

lena_test_error<-testerror(lena,lena_test_hat,lena_test_omega)
lena_test_error
#0.03104175
library(Matrix)
rankMatrix(lena)
rankMatrix(lena_test_hat)
#original image
image(lena[256:1,256:1],col = grey(seq(0, 1, length = 256)))

#before pic amendment
image(lena_test[256:1,256:1],col = grey(seq(0, 1, length = 256)))

#after pic amendment
image(lena_test_hat[256:1,256:1],col = grey(seq(0, 1, length = 256)))


