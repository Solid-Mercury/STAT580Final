##########################################################
#####        functions to do comparison      #############
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

############################################################
############ Rcpp Armadillo##############################
library(Rcpp)
library(RcppArmadillo)
sourceCpp("armadillo_svd.cpp")
z_new_armadillo_svd<-function(z,z_old,lambda){
	r_zero<-row(z)[which(z!=0)]
	c_zero<-col(z)[which(z!=0)]
	z_old[cbind(r_zero,c_zero)]<-0
	p<-z+z_old
	p<-armadillo_svd(p)
	D<-diag(as.vector(p$s))
	U<-(p$U)
	V<-(p$V)
	D<-D-lambda*diag(length(p$s))
	D[D<0]<-0
	Da<-matrix(rep(0,dim(U)[2]*dim(V)[2]),nrow=dim(U)[2])
	Da[1:dim(D)[1],1:dim(D)[2]]<-Da[1:dim(D)[1],1:dim(D)[2]]+D
	D<-Da
	z_new<-U%*%D%*%t(V)
	return (z_new)
}
armadillo_cpp_svd<-function(z,lambda){
 z_hat<-array(rep(0,nrow(z)*ncol(z)*length(lambda)),dim=c(nrow(z),ncol(z),length(lambda)))
 z_old<-matrix(rep(0,nrow(z)*ncol(z)),nrow=nrow(z))
for (i in 1:length(lambda)){
		z_new<-z_new_armadillo_svd(z,z_old,lambda[i])
		diff=1
	while(diff>0.00001){
		z_old<-z_new
		z_new<-z_new_armadillo_svd(z,z_old,lambda[i])
		diff<-sum((z_new-z_old)^2)/(sum(z_old^2))
        if (is.na(diff)){
        break
      }
	}
z_hat[,,i]<-z_new}
return (z_hat)
}

############################################################
################ read in data###############################
m_test<-read.table("test.txt")
m_training<-read.table("trainng.txt")
m_validation<-read.table("validating.txt")
m_test_omega<-read.table("test_omega.txt")
m_validation_omega<-read.table("validating.txt")
m_training_omega<-read.table("trainng_omega.txt")
m_validation2<-read.table("validating_2.txt")
m_validation2_omega<-read.table("validating_2_omega.txt")

m_test[is.na(m_test)]<-0

m_training[is.na(m_training)]<-0

m_validation[is.na(m_validation)]<-0
m_validation2[is.na(m_validation2)]<-0
m_validation_omega[is.na(m_validation_omega)]<-0
m_validation2_omega[is.na(m_validation2_omega)]<-0


m_training<-as.matrix(m_training)
m_test<-as.matrix(m_test)
m_validation<-as.matrix(m_validation)
m_validation2<-as.matrix(m_validation2)
m_validation_omega<-as.matrix(m_validation_omega)
m_validation2_omega<-as.matrix(m_validation2_omega)
m_test_omega<-as.matrix(m_test_omega)
m_training_omega<-as.matrix(m_training_omega)


lambda=seq(100,10,-10)




############comparing lambda##################
lambda_com<-function(td,vd,lambda,td_omega){
td_hat<-array(rep(0,nrow(td)*ncol(td)*length(lambda)),dim=c(nrow(td),ncol(td),length(lambda)))
td_error<-rep(0,length(lambda))
for (i in 1:length(lambda)){
	td_hat[,,i]<-armadillo_cpp_svd(td,lambda[i])
	td_error[i]<-trainingerror(vd,td_hat[,,i], td_omega)	
}
best_lambda<-lambda[which.min(td_error)]
min_td_error<-td_error
return (list(best_lambda,min_td_error))
}

#ts<-proc.time()
l<-lambda_com(m_training,m_validation,lambda,m_validation_omega)
t<-system.time(lambda_com(m_training,m_validation,lambda,m_validation_omega))
#t<-proc.time()-ts


best_lambda<-unlist(l[1])

min_td_error<-unlist(l[2])

o<-list(t,best_lambda,min_td_error)
save(o, file = "p_cpp_100_10.Rda" )










