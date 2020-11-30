#include <RcppEigen.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::normal;
using boost::math::students_t;
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
//  return -(x.size() * log(2 * M_PI) + x.cwiseProduct(x).sum()) / 2;
  return -x.cwiseProduct(x).sum() / 2.0;
}

Eigen::VectorXd qt(const Eigen::VectorXd& p, const double nu){
  boost::math::students_t t(nu);
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++){
    out(i) = quantile(t, p.coeff(i));
  }
  return out;
}

double C_student(const double nu){ // useless !!!
  return boost::math::lgamma((nu+1.0)/2.0) -
    boost::math::lgamma(nu/2.0) - log(nu*M_PI)/2.0;
}

double logdt(const Eigen::VectorXd& x, const double nu){
/*  return x.size() *
    (boost::math::lgamma((nu+1.0)/2.0) -
    boost::math::lgamma(nu/2.0) - log(nu*M_PI)/2.0) -
    (nu+1.0)/2.0 * log1p(x.cwiseProduct(x).array() / nu).sum();
*/
  return -(nu+1.0)/2.0 * log1p(x.cwiseProduct(x).array() / nu).sum();
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
      Eigen::VectorXd v = ymI - XmI * theta.topRows(q-1);
      J(counter) = logdnorm(v/sigma) - (n-q) * log(sigma);
      Theta.col(counter) = theta;
      counter++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Theta") = Theta.transpose(),
                            Rcpp::Named("logWeights") = J);
}

// [[Rcpp::export]]
Rcpp::List f_student(
    const Eigen::MatrixXd& centers,
    const Eigen::MatrixXd& XI,
    const Eigen::MatrixXd& XmI,
    const Eigen::VectorXd& yI,
    const Eigen::VectorXd& ymI,
    const size_t M,
    const size_t n,
    const double nu
){
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  //double C = C_student(nu);
  for(size_t m = 0; m < ncenters; m++){
    Eigen::MatrixXd H(q,q);
    H << XI, qt(centers.col(m), nu);
    Eigen::MatrixXd Ht = H.transpose();
    Eigen::VectorXd theta = (Ht * H).inverse() * Ht * yI;
    double sigma = theta.coeff(q-1);
    if(sigma > 0){
      Eigen::VectorXd v = ymI - XmI * theta.topRows(q-1);
      J(counter) = logdt(v/sigma, nu) - (n-q) * log(sigma);
      Theta.col(counter) = theta;
      counter++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Theta") = Theta.transpose(),
                            Rcpp::Named("logWeights") = J);
}
