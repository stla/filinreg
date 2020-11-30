#include <RcppEigen.h>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
double test(double x){
  return x + 1;
}

Eigen::VectorXd qnorm(const Eigen::VectorXd& p){
  boost::math::normal gaussian;
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++){
    out(i) = quantile(gaussian, p.coeff(i));
  }
  return out;
}

double logdnorm(const Eigen::VectorXd& x){
  return -(x.size() * log(2 * M_PI) + x.cwiseProduct(x).sum()) / 2;
}

// [[Rcpp::export]]
Rcpp::List f_normal(
  const Eigen::MatrixXd& centers,
  const Eigen::MatrixXd& XI,
  const Eigen::MatrixXd& XmI,
  const Eigen::VectorXd& yI,
  const Eigen::VectorXd& ymI,
  const size_t M,
  const size_t n
){
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  for(size_t m = 0; m < ncenters; m++){
    Eigen::MatrixXd H(q,q);
    H << XI, qnorm(centers.col(m));
    Eigen::MatrixXd Ht = H.transpose();
    Eigen::VectorXd theta = (Ht * H).inverse() * Ht * yI;
    double sigma = theta.coeff(q-1);
    if(sigma > 0){
      Eigen::VectorXd beta = theta.topRows(q-1);
      Eigen::VectorXd v = ymI - XmI * beta;
      J(counter) = logdnorm(v/sigma) - (n-q) * log(sigma);
      Theta.col(counter) = theta;
      counter++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Theta") = Theta.transpose(),
                            Rcpp::Named("logWeights") = J);
}