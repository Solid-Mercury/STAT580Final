z<-scan(file="~/Documents/isu/STAT 580/final_proj/lena256",skip=1)
r_index<-seq(1,length(z),by=1)
r_index<-sample(r_index)
r_index_zero<-r_index[1:(length(z)*0.4)]
z[r_index_zero]<-0
z<-matrix(z,nrow=sqrt(length(z)),byrow=T)
image(z)
lambda<-seq(1,0.2,by=-0.2)