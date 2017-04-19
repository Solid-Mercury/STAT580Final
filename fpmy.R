lena <- scan("C:\\Users\\Lynn\\Desktop\\Courses\\stat 510\\lena256")
x <- lena
for(i in 1:256^2){
  if (runif(1) > 0.6){
    x[i] = 0
  }
}
x <- matrix(x, nrow = 256)
x <- x[256:1, 256:1]
lena <- matrix(lena, nrow = 256)
lena <- lena[256:1, 256:1]
mysoft <- function(x,zold,lambda,rank.max){
  res = matrix(rep(0, nrow(x)*ncol(x)), nrow = nrow(x))
  res[x!=0] = res[x!=0] + x[x!=0]
  res[x == 0]= res[x==0] + zold[x == 0]
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





lambda = 90
sbchengxu <- function (x,lambda,rank.max){
zhat <- array(rep(0, nrow(x) * ncol(x)*length(lambda)), dim = c(nrow(x), ncol(x),length(lambda)))
zold = matrix(rep(0, nrow(x) * ncol(x)), nrow = nrow(x))
znew = mysoft(x,zold,lambda[1],rank.max)
zold = znew
for (i in 1:length(lambda)){
  while(1){
  znew = mysoft(x,zold,lambda[i],rank.max)
  diff=(sum((znew-zold)^2))/(sum(zold^2))
  print(diff)
  if(diff < 10^-5){
    break
  }
  zold <- znew
  }
  zhat[,,i] <- znew
}
return (zhat)
}


out <- sbchengxu(x,lambda,255)
image(out[,,1], col = gray(0:255/255))
