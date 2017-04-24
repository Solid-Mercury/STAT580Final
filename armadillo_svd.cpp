#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List armadillo_svd(arma::mat &p){
   arma::mat U; 
   arma::vec s;
   arma::mat V;
   
   svd(U,s,V,p);

   List ret;
   ret["U"]=U;
   ret["s"]=s;
   ret["V"]=V;
   

   return (ret);
}