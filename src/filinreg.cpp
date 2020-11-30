#include <RcppEigen.h>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
double test(double x){
  return x + 1;
}
