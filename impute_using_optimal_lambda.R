best_lambda <- 13

m_test_hat_rcpp<-armadillo_cpp_svd(m_test,best_lambda)[,,1]

save(m_test_hat,file="m_test_hat_rcpp.Rda")

m_test_hat_irlba<-irlba_svd(m_test,best_lambda)[,,1]

save(m_test_hat_irlba,file="m_test_hat_irlbaRda")