load("m_test_hat_irlba.Rda")
m_test_hat_irlba <- m_test_hat
load("m_test_hat_cpp.Rda")
m_test_hat_cpp <- m_test_hat
test_set <- read.table("./Movie data/validating_2.txt")
test_set <- as.matrix(test_set)
test_set[is.na(test_set)] <- 0
test_set_omega <- read.table("./Movie data/validating_2_omega.txt")
test_set_omega <- as.matrix(test_set_omega)


# result from irlba
test_error_irlba <- (norm(m_test_hat_irlba*test_set_omega - test_set, "F"))^2/(norm(test_set, "F"))^2

RMSE_irlba <- sqrt((norm(m_test_hat_irlba*test_set_omega - test_set, "F"))^2/(sum(test_set_omega)))

# result from Rcpp
test_error_cpp <- (norm(m_test_hat_cpp*test_set_omega - test_set, "F"))^2/(norm(test_set, "F"))^2

RMSE_cpp <- sqrt((norm(m_test_hat_cpp*test_set_omega - test_set, "F"))^2/(sum(test_set_omega)))


