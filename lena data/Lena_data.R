
lena <- scan("/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/lena256", skip = 1)
x <- rep(1:256, 256)
y <- rep(1:256, each = 256)
data <- as.data.frame(cbind(x, y, lena))

x <- 1:(256*256)
x1 <- sample(x, round(256*256*0.6), replace = F)
x2 <- sample(x1, round(256*256*0.6*0.7), replace = F)
x3 <- x1[!x1 %in% x2]

training <- matrix(nrow = 256, ncol = 256)
validating <- matrix(nrow = 256, ncol = 256)
test <- matrix(nrow = 256, ncol = 256)
training_omega <- matrix(0, nrow = 256, ncol = 256)
validating_omega <- matrix(0, nrow = 256, ncol = 256)
test_omega <- matrix(0, nrow = 256, ncol = 256)

for (i in 1:length(x1)){
  u <- data[x1[i],1]
  v <- data[x1[i],2]
  w <- data[x1[i],3]
  test[u,v] <- w
  test_omega[u,v] <- 1
}

for (i in 1:length(x2)){
  u <- data[x2[i],1]
  v <- data[x2[i],2]
  w <- data[x2[i],3]
  training[u,v] <- w
  training_omega[u,v] <- 1
}

for (i in 1:length(x3)){
  u <- data[x3[i],1]
  v <- data[x3[i],2]
  w <- data[x3[i],3]
  validating[u,v] <- w
  validating_omega[u,v] <- 1
}

write.table(training, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/trainng.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/validating.txt", sep = "\t", col.names = F, row.names = F)
write.table(test, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/test.txt", sep = "\t", col.names = F, row.names = F)
write.table(training_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/trainng_omega.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/validating_omega.txt", sep = "\t", col.names = F, row.names = F)
write.table(test_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/lena/test_omega.txt", sep = "\t", col.names = F, row.names = F)


image(test[256:1, 256:1])
image(training[256:1, 256:1])
image(validating[256:1, 256:1])


