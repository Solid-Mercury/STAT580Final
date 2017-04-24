data <- read.table(file = "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/ml-100k/u.data")
training <- matrix(nrow = 943, ncol = 1682)
validating <- matrix(nrow = 943, ncol = 1682)
test <- matrix(nrow = 943, ncol = 1682)
training_omega <- matrix(0, nrow = 943, ncol = 1682)
validating_omega <- matrix(0, nrow = 943, ncol = 1682)
test_omega <- matrix(0, nrow = 943, ncol = 1682)
validating_2 <- matrix(0, nrow = 943, ncol = 1682)
validating_2_omega <- matrix(0, nrow = 943, ncol = 1682)

x <- 1:100000
x1 <- sample(x, 70000, replace = F)
x2 <- sample(x[-x1], 15000, replace = F)
x3 <- x[-c(x1,x2)]

for (i in 1:70000){
  custom_id <- data[x1[i],1]
  movie_id <- data[x1[i],2]
  score <- data[x1[i],3]
  training[custom_id, movie_id] <- score
  training_omega[custom_id, movie_id] <- 1
}

for (i in 1:15000){
  custom_id <- data[x2[i],1]
  movie_id <- data[x2[i],2]
  score <- data[x2[i],3]
  validating[custom_id, movie_id] <- score
  validating_omega[custom_id, movie_id] <- 1
}

test <- ifelse(is.na(training),ifelse(is.na(validating),NA,validating), ifelse(is.na(validating), training, training+validating)) 
test_omega <- training_omega + validating_omega

for (i in 1:15000){
  custom_id <- data[x3[i],1]
  movie_id <- data[x3[i],2]
  score <- data[x3[i],3]
  validating_2[custom_id, movie_id] <- score
  validating_2_omega[custom_id, movie_id] <- 1
}



write.table(training, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/trainng.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/validating.txt", sep = "\t", col.names = F, row.names = F)
write.table(test, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/test.txt", sep = "\t", col.names = F, row.names = F)
write.table(training_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/trainng_omega.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/validating_omega.txt", sep = "\t", col.names = F, row.names = F)
write.table(test_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/test_omega.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating_2, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/validating_2.txt", sep = "\t", col.names = F, row.names = F)
write.table(validating_2_omega, "/Users/ganghan/Documents/Study/Statistics/Stat 580/Final project/100k/validating_2_omega.txt", sep = "\t", col.names = F, row.names = F)

