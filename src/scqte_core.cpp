# include <RcppArmadillo.h>
# include <cmath>
//# include "basicOp.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double mad(const arma::vec& x) {
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx1, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) * sx1(i);
  }
  return X;
}

// Asymmetric huber regression adjusted to quantile tau for initialization 
// [[Rcpp::export]]
void updateHuber(const arma::mat& Z, const arma::vec& res, const double tau, arma::vec& der, arma::vec& grad, const int n, const double rob, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur > rob) {
      der(i) = -2 * tau * rob;
    } else if (cur > 0) {
      der(i) = -2 * tau * cur;
    } else if (cur > -rob) {
      der(i) = 2 * (tau - 1) * cur;
    } else {
      der(i) = 2 * (1 - tau) * rob;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void updateGauss(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const double tau, const double n1, const double h1) {
  der = arma::normcdf(-res * h1) - tau;
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void hessGauss(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::mat& hess, const double tau, const double n1, const double h1) {
  Environment base("package:base");
  Function sweep = base["sweep"];
  der = arma::normpdf(res * h1) * h1;
  arma::mat WZ = as<arma::mat>(sweep(Z, 1, der,"*"));
  hess = n1 * Z.t() * WZ;
}

// [[Rcpp::export]]
arma::vec Loggrad(arma::mat X, arma::vec y, arma::vec beta){
  int n = X.n_rows;
  arma::vec Xb = X*beta;
  arma::vec mu = exp(Xb)/(exp(Xb)+1);
  arma::vec grad = X.t()*(y-mu)/n;
  return(grad);
}

// [[Rcpp::export]]
arma::mat Loghess(arma::mat X, arma::vec beta){
  Environment base("package:base");
  Function sweep = base["sweep"];
  int n = X.n_rows;
  arma::vec Xb = X*beta;
  arma::vec mu = exp(Xb)/(exp(Xb)+1);
  arma::vec sigma = mu%(1-mu);
  arma::mat XW = as<arma::mat>(sweep(X, 1, sigma, "*"));
  arma::mat hess = X.t()*XW/n;
  return(hess);
}

// Functions with an upper bound for the GD step size.
// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& Z, const arma::vec& Y, const double tau, arma::vec& der, arma::vec& gradOld, arma::vec& gradNew, const int n, const int p, 
                   const double n1, const double tol = 0.0001, const double constTau = 1.345, const int iteMax = 5000, const double stepMax = 100.0) {
  double rob = constTau * mad(Y);
  updateHuber(Z, Y, tau, der, gradOld, n, rob, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  rob = constTau * mad(res);
  updateHuber(Z, res, tau, der, gradNew, n, rob, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), stepMax);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    rob = constTau * mad(res);
    updateHuber(Z, res, tau, der, gradNew, n, rob, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// Low-dimensional conquer: estimation
// [[Rcpp::export]]
Rcpp::List SQRcpp(const arma::mat& X, arma::vec Y, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                     const double tol = 0.0001, const int iteMax = 5000, const double stepMax = 100.0) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec beta = huberReg(Z, Y, tau, der, gradOld, gradNew, n, p, n1, tol, constTau, iteMax, stepMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, tau, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, tau, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), stepMax);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, tau, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) %= sx1;
  beta(0) += my - arma::as_scalar(mx * beta.rows(1, p));
  return Rcpp::List::create(Rcpp::Named("coeff") = beta, Rcpp::Named("ite") = ite, Rcpp::Named("residual") = res, Rcpp::Named("bandwidth") = h);
}

// [[Rcpp::export]]
Rcpp::List seqSQRcpp(const arma::mat& X, arma::vec Y, arma::vec sizeM, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                     const double tol = 0.0001, const int iteMax = 5000, const double stepMax = 100.0) {
  const int n = X.n_rows, p = X.n_cols, M = sizeM.size();
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double h1 = 1.0 / h;
  double nm1 = 1.0 / sizeM(0), record = sizeM(0);
  arma::mat betaiter(p+1, M), betarecord(p+1, M), hessiter(p+1, p+1), hess(p+1, p+1);
  arma::vec gradOld(p+1), gradNew(p+1), betare(p+1), der1(sizeM(0)), iterecord(M-1);
  arma::mat X1 = X.rows(0, sizeM(0)-1);
  arma::vec Y1 = Y.subvec(0, sizeM(0)-1);
  arma::vec beta = SQRcpp(X1, Y1, tau, h)["coeff"];
  betarecord.col(0) = beta;
  arma::rowvec mx1 = arma::mean(X1, 0), mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X1, 0, 0).t(), sx = 1.0 / arma::stddev(X, 0, 0).t();
  double my1 = arma::mean(Y1), my = arma::mean(Y);
  arma::mat Z1 = arma::join_rows(arma::ones(sizeM(0)), standardize(X1, mx1, sx1, p));
  Y1 -= my1;
  beta(0) -= my1 - arma::as_scalar(mx1 * beta.rows(1, p));
  beta.rows(1, p) %= arma::stddev(X1, 0, 0).t();
  betaiter.col(0) = beta;
  arma::vec res1 = Y1 - Z1 * beta;
  hessGauss(Z1, res1, der1, hess, tau, nm1, h1);
  hessiter = hess;
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  Y -= my;
  for(int m = 1; m < M; m++){
    arma::mat Zm = Z.rows(record, record+sizeM(m)-1);
    arma::vec Ym = Y.subvec(record, record+sizeM(m)-1);
    arma::vec resm = Ym - Zm * beta;
    arma::vec derm(sizeM(m));
    nm1 = 1.0 / sizeM(m);
    updateGauss(Zm, resm, derm, gradOld, tau, nm1, h1);
    beta -= gradOld;
    arma::vec betaDiff = -gradOld;
    resm -= Zm * betaDiff;
    updateGauss(Zm, resm, derm, gradNew, tau, nm1, h1);
    gradNew += hessiter*(beta - betaiter.col(m-1));
    arma::vec gradDiff = gradNew - gradOld;
    int ite = 1;
    while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
      double alpha = 1.0;
      double cross = arma::as_scalar(betaDiff.t() * gradDiff);
      if (cross > 0) {
        double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
        double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
        alpha = std::min(std::min(a1, a2), stepMax);
      }
      gradOld = gradNew;
      betaDiff = -alpha * gradNew;
      beta += betaDiff;
      resm -= Zm * betaDiff;
      updateGauss(Zm, resm, derm, gradNew, tau, nm1, h1);
      gradNew += hessiter*(beta- betaiter.col(m-1));
      gradDiff = gradNew - gradOld;
      ite++;
    }
    hessGauss(Zm, resm, derm, hess, tau, nm1, h1);
    hessiter += hess;
    record += sizeM(m);
    betaiter.col(m) = beta;
    betare.rows(1, p) = beta.rows(1, p)%sx;
    betare(0) = beta(0)+my-arma::as_scalar(mx * betare.rows(1, p));
    betarecord.col(m) = betare;
    iterecord(m-1) = ite;
  }
  return Rcpp::List::create(Rcpp::Named("coeff") = betarecord.col(M-1), Rcpp::Named("record") = betarecord, Rcpp::Named("ite") = iterecord, Rcpp::Named("bandwidth") = h);
}

// [[Rcpp::export]]
arma::vec Newtoncpp(arma::mat X, arma::vec y, arma::vec betaini, double eps = 0.000001, int maxstep = 100){
  double dis = 2*eps;
  int iteration = 0;
  int p = betaini.size();
  arma::vec beta = arma::zeros(p);
  while((dis > eps)&&(iteration < maxstep)){
    beta = betaini+inv(Loghess(X,betaini))*Loggrad(X,y,betaini);
    dis = sum(abs(beta-betaini));
    betaini = beta;
    iteration = iteration+1;
  }
  return(beta);
}

// [[Rcpp::export]]
arma::vec seqNewtoncpp(arma::mat X, arma::vec y, arma::vec betaini, arma::vec sizeM, double eps = 0.000001, int maxstep = 100){
  
  int p = betaini.size();
  int M = sizeM.size();
  int record = 0;
  arma::mat betarecord = arma::zeros(p, M);
  arma::mat X1 = X.rows(0, sizeM(0)-1);
  arma::vec y1 = y.subvec(0, sizeM(0)-1);
  betarecord.col(0) = Newtoncpp(X1, y1, betaini, eps, maxstep);  
  record = record+sizeM(0);
  arma::mat hessrecord = Loghess(X1, betarecord.col(0));
  for(int m = 1; m < M; m++){
    arma::mat Xm = X.rows(record, record+sizeM(m)-1);
    arma::vec ym = y.subvec(record, record+sizeM(m)-1);
    arma::vec betamini = betarecord.col(m-1);
    arma::mat Finv = inv(hessrecord+Loghess(Xm,betamini));
    double dism = 2*eps;
    int iterm = 0;
    arma::vec betam = betamini;
    while((dism>eps)&(iterm<maxstep)){
      betam = betamini+Finv*(Loggrad(Xm,ym,betamini)+hessrecord*(betarecord.col(m-1)-betamini));
      dism = sum(abs(betam-betamini));
      betamini = betam;
      iterm = iterm+1;
    }
    betarecord.col(m) = betam;
    hessrecord = hessrecord+Loghess(Xm, betam);
    record = record + sizeM(m);
  }
  return(betarecord.col(M-1));
}

// [[Rcpp::export]]
Rcpp::List SQRcppW(const arma::mat& Z, arma::vec Y, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                       const double tol = 0.0001, const int iteMax = 5000, const double stepMax = 100.0) {
  const int n = Z.n_rows, p = Z.n_cols;
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double h1 = 1.0 / h, n1 = 1.0 / n;
  arma::vec gradOld(p), gradNew(p), der(n);
  arma::vec beta = huberReg(Z, Y, tau, der, gradOld, gradNew, n, p, n1, tol, constTau, iteMax, stepMax);
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, tau, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, tau, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), stepMax);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, tau, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return Rcpp::List::create(Rcpp::Named("coeff") = beta, Rcpp::Named("ite") = ite, Rcpp::Named("residual") = res, Rcpp::Named("bandwidth") = h);
}

// [[Rcpp::export]]
Rcpp::List seqSQRcppW(const arma::mat& Z, arma::vec Y, arma::vec sizeM, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                     const double tol = 0.0001, const int iteMax = 5000, const double stepMax = 100.0) {
  const int n = Z.n_rows, p = Z.n_cols, M = sizeM.size();
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n/M) + p) / (n/M), 0.4), 0.05);
  }
  const double h1 = 1.0 / h;
  double nm1 = 1.0 / sizeM(0), record = sizeM(0);
  arma::mat betaiter(p, M), hessiter(p, p), hess(p,p);
  arma::vec gradOld(p), gradNew(p), der1(sizeM(0)), iterecord(M-1);
  arma::mat Z1 = Z.rows(0, sizeM(0)-1);
  arma::vec Y1 = Y.subvec(0, sizeM(0)-1);
  arma::vec beta = SQRcppW(Z1, Y1, tau, h)["coeff"];
  betaiter.col(0) = beta;
  arma::vec res1 = Y1 - Z1 * beta;
  hessGauss(Z1, res1, der1, hess, tau, nm1, h1);
  hessiter = hess;
  for(int m = 1; m < M; m++){
    arma::mat Zm = Z.rows(record, record+sizeM(m)-1);
    arma::vec Ym = Y.subvec(record, record+sizeM(m)-1);
    arma::vec resm = Ym - Zm * beta;
    arma::vec derm(sizeM(m));
    nm1 = 1.0 / sizeM(m);
    updateGauss(Zm, resm, derm, gradOld, tau, nm1, h1);
    beta -= gradOld;
    arma::vec betaDiff = -gradOld;
    resm -= Zm * betaDiff;
    updateGauss(Zm, resm, derm, gradNew, tau, nm1, h1);
    gradNew += hessiter*(beta - betaiter.col(m-1));
    arma::vec gradDiff = gradNew - gradOld;
    int ite = 1;
    while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
      double alpha = 1.0;
      double cross = arma::as_scalar(betaDiff.t() * gradDiff);
      if (cross > 0) {
        double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
        double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
        alpha = std::min(std::min(a1, a2), stepMax);
      }
      gradOld = gradNew;
      betaDiff = -alpha * gradNew;
      beta += betaDiff;
      resm -= Zm * betaDiff;
      updateGauss(Zm, resm, derm, gradNew, tau, nm1, h1);
      gradNew += hessiter*(beta- betaiter.col(m-1));
      gradDiff = gradNew - gradOld;
      ite++;
    }
    betaiter.col(m) = beta;
    hessGauss(Zm, resm, derm, hess, tau, nm1, h1);
    hessiter += hess;
    record += sizeM(m);
    iterecord(m-1) = ite;
  }
  return Rcpp::List::create(Rcpp::Named("coeff") = beta, Rcpp::Named("record") = betaiter, Rcpp::Named("ite") = iterecord, Rcpp::Named("bandwidth") = h);
}

// [[Rcpp::export]]
Rcpp::List seqSQRcppW1R(const arma::mat& Z, arma::vec A, arma::vec Y, arma::vec sizeM, arma::vec betaAini, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                        const double tol = 0.0001, const double tolN = 0.000001, const int iteMax = 5000, const int iteMaxN = 100, const double stepMax = 100.0) {
  const int n = Z.n_rows, p = Z.n_cols, M = sizeM.size();
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n/M) + p) / (n/M), 0.4), 0.05);
  }
  const double h1 = 1.0 / h;
  double nm1 = 1.0 / sizeM(0);
  int record = sizeM(0);
  arma::mat betaAiter(p, M), hessAiter(p, p), betaiter(2, M), hessiter(2, 2), hess(2, 2);
  arma::vec gradOld(2), gradNew(2), der1(sizeM(0)), iterecord(M-1);
  arma::mat Z1 = Z.rows(0, sizeM(0)-1);
  arma::vec A1 = A.subvec(0, sizeM(0)-1);
  arma::vec betaA0 = Newtoncpp(Z1, A1, betaAini, tolN, iteMaxN);  
  betaAiter.col(0) = betaA0;
  hessAiter = Loghess(Z1, betaA0);
  arma::vec muAseq = exp(Z1*betaA0)/(1+exp(Z1*betaA0));
  arma::vec Wseq = A1/muAseq+(1-A1)/(1-muAseq);
  arma::vec Yw1 = Wseq%Y.subvec(0, sizeM(0)-1);
  arma::mat Aw1 = arma::join_rows(Wseq, Wseq%A1);
  arma::vec beta = SQRcppW(Aw1, Yw1, tau, h)["coeff"];
  betaiter.col(0) = beta;
  arma::vec res1 = Yw1 - Aw1 * beta;
  hessGauss(Aw1, res1, der1, hess, tau, nm1, h1);
  hessiter = hess;
  for(int m = 1; m < M; m++){
    arma::mat Zm = Z.rows(record, record+sizeM(m)-1);
    arma::vec Am = A.subvec(record, record+sizeM(m)-1);
    arma::vec betaAmini = betaAiter.col(m-1);
    arma::mat Finv = inv(hessAiter+Loghess(Zm,betaAmini));
    double dism = 2*tolN;
    int item = 0;
    arma::vec betaAm = betaAmini;
    while((dism>tolN)&(item<iteMaxN)){
      betaAm = betaAmini+Finv*(Loggrad(Zm,Am,betaAmini)+hessAiter*(betaAiter.col(m-1)-betaAmini));
      dism = sum(abs(betaAm-betaAmini));
      betaAmini = betaAm;
      item += 1;
    }
    betaAiter.col(m) = betaAm;
    hessAiter += Loghess(Zm, betaAm);
    arma::vec muAseqm = exp(Zm*betaAm)/(1+exp(Zm*betaAm));
    arma::vec Wseqm = Am/muAseqm+(1-Am)/(1-muAseqm);
    arma::vec Ywm = Wseqm%Y.subvec(record, record+sizeM(m)-1);
    arma::mat Awm = arma::join_rows(Wseqm, Wseqm%Am); 
    arma::vec resm = Ywm - Awm * beta;
    arma::vec derm(sizeM(m));
    nm1 = 1.0 / sizeM(m);
    updateGauss(Awm, resm, derm, gradOld, tau, nm1, h1);
    beta -= gradOld;
    arma::vec betaDiff = -gradOld;
    resm -= Awm * betaDiff;
    updateGauss(Awm, resm, derm, gradNew, tau, nm1, h1);
    gradNew += hessiter*(beta - betaiter.col(m-1));
    arma::vec gradDiff = gradNew - gradOld;
    int ite = 1;
    while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
      double alpha = 1.0;
      double cross = arma::as_scalar(betaDiff.t() * gradDiff);
      if (cross > 0) {
        double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
        double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
        alpha = std::min(std::min(a1, a2), stepMax);
      }
      gradOld = gradNew;
      betaDiff = -alpha * gradNew;
      beta += betaDiff;
      resm -= Awm * betaDiff;
      updateGauss(Awm, resm, derm, gradNew, tau, nm1, h1);
      gradNew += hessiter*(beta- betaiter.col(m-1));
      gradDiff = gradNew - gradOld;
      ite++;
    }
    betaiter.col(m) = beta;
    hessGauss(Awm, resm, derm, hess, tau, nm1, h1);
    hessiter += hess;
    record += sizeM(m);
    iterecord(m-1) = ite;
  }
  return Rcpp::List::create(Rcpp::Named("coeff") = beta, Rcpp::Named("record") = betaiter, Rcpp::Named("ite") = iterecord, Rcpp::Named("bandwidth") = h);
}


// [[Rcpp::export]]
Rcpp::List seqSQRcppW1Ren(const arma::mat& Z, arma::vec D, arma::vec A, arma::vec Y, arma::vec sizeM, arma::vec betaDini, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                        const double tol = 0.0001, const double tolN = 0.000001, const int iteMax = 5000, const int iteMaxN = 100, const double stepMax = 100.0) {

  const int n = Z.n_rows, p = Z.n_cols, M = sizeM.size();
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n/M) + p) / (n/M), 0.4), 0.05);
  }
  const double h1 = 1.0 / h;
  double nm1 = 1.0 / sizeM(0);
  int record = sizeM(0);
  arma::mat betaDiter(p, M), hessDiter(p, p), betaiter(2, M), hessiter(2, 2), hess(2, 2);
  arma::vec gradOld(2), gradNew(2), der1(sizeM(0)), iterecord(M-1);
  arma::mat Z1 = Z.rows(0, sizeM(0)-1);
  arma::vec D1 = D.subvec(0, sizeM(0)-1);
  arma::vec A1 = A.subvec(0, sizeM(0)-1);
  arma::vec betaD0 = Newtoncpp(Z1, D1, betaDini, tolN, iteMaxN);  
  betaDiter.col(0) = betaD0;
  hessDiter = Loghess(Z1, betaD0);
  arma::vec muDseq = exp(Z1*betaD0)/(1+exp(Z1*betaD0));
  arma::vec Wseq = (D1-muDseq)%(2*A1-1)/(muDseq%(1-muDseq));
  arma::vec Yw1 = Wseq%Y.subvec(0, sizeM(0)-1);
  arma::mat Aw1 = arma::join_rows(Wseq, Wseq%Z1.col(1),Wseq%Z1.col(2), Wseq%A1);
  arma::vec beta = SQRcppW(Aw1, Yw1, tau, h)["coeff"];
  betaiter.col(0) = beta;
  arma::vec res1 = Yw1 - Aw1 * beta;
  hessGauss(Aw1, res1, der1, hess, tau, nm1, h1);
  hessiter = hess;
  for(int m = 1; m < M; m++){
    arma::mat Zm = Z.rows(record, record+sizeM(m)-1);
    arma::vec Am = A.subvec(record, record+sizeM(m)-1);
    arma::vec Dm = D.subvec(record, record+sizeM(m)-1);
    arma::vec betaDmini = betaDiter.col(m-1);
    arma::mat Finv = inv(hessDiter+Loghess(Zm,betaDmini));
    double dism = 2*tolN;
    int item = 0;
    arma::vec betaDm = betaDmini;
    while((dism>tolN)&(item<iteMaxN)){
      betaDm = betaDmini+Finv*(Loggrad(Zm,Dm,betaDmini)+hessDiter*(betaDiter.col(m-1)-betaDmini));
      dism = sum(abs(betaDm-betaDmini));
      betaDmini = betaDm;
      item += 1;
    }
    betaDiter.col(m) = betaDm;
    hessDiter += Loghess(Zm, betaDm);
    arma::vec muDseqm = exp(Zm*betaDm)/(1+exp(Zm*betaDm));
    arma::vec Wseqm = (Dm-muDseqm)%(2*Am-1)/(muDseqm%(1-muDseqm));
    arma::vec Ywm = Wseqm%Y.subvec(record, record+sizeM(m)-1);
    arma::mat Awm = arma::join_rows(Wseqm, Wseqm%Zm.col(1), Wseqm%Zm.col(2), Wseqm%Am); 
    arma::vec resm = Ywm - Awm * beta;
    arma::vec derm(sizeM(m));
    nm1 = 1.0 / sizeM(m);
    updateGauss(Awm, resm, derm, gradOld, tau, nm1, h1);
    beta -= gradOld;
    arma::vec betaDiff = -gradOld;
    resm -= Awm * betaDiff;
    updateGauss(Awm, resm, derm, gradNew, tau, nm1, h1);
    gradNew += hessiter*(beta - betaiter.col(m-1));
    arma::vec gradDiff = gradNew - gradOld;
    int ite = 1;
    while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
      double alpha = 1.0;
      double cross = arma::as_scalar(betaDiff.t() * gradDiff);
      if (cross > 0) {
        double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
        double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
        alpha = std::min(std::min(a1, a2), stepMax);
      }
      gradOld = gradNew;
      betaDiff = -alpha * gradNew;
      beta += betaDiff;
      resm -= Awm * betaDiff;
      updateGauss(Awm, resm, derm, gradNew, tau, nm1, h1);
      gradNew += hessiter*(beta- betaiter.col(m-1));
      gradDiff = gradNew - gradOld;
      ite++;
    }
    betaiter.col(m) = beta;
    hessGauss(Awm, resm, derm, hess, tau, nm1, h1);
    hessiter += hess;
    record += sizeM(m);
    iterecord(m-1) = ite;
  }
  
  return Rcpp::List::create(Rcpp::Named("coeff") = beta, Rcpp::Named("record") = betaiter, Rcpp::Named("ite") = iterecord, Rcpp::Named("bandwidth") = h);
}



// [[Rcpp::export]]
arma::vec SQRcppWIni(const arma::mat& Z, arma::vec Y, arma::vec betaHat, const double tau = 0.5, double h = 0.05, const double constTau = 1.345, 
                      const double tol = 0.0001, const int iteMax = 5000, const double stepMax = 100.0) {
  const int n = Z.n_rows, p = Z.n_cols;
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double h1 = 1.0 / h, n1 = 1.0 / n;
  arma::vec gradOld(p), gradNew(p), der(n);
  arma::vec beta = betaHat;
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, tau, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, tau, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), stepMax);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, tau, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::mat SQRWINF(arma::mat X, arma::vec Y, arma::vec betaHat, int n, int p, double h = 0.05, double tau = 0.5, 
                       int B = 1000, double tol = 0.0001, int iteMax = 5000, double stepMax = 100.0){
  arma::mat rst(p, B);
  if (h <= 0.0) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  for (int b = 0; b < B; b++) {
    arma::uvec idx = arma::find(arma::randi(n, arma::distr_param(0, 1)) == 1);
    arma::mat mbX = X.rows(idx);
    arma::mat mbY = Y.rows(idx);
    rst.col(b) = SQRcppWIni(mbX, mbY, betaHat, tau, h);
  }
  return rst;
}


// [[Rcpp::export]]
arma::mat asymptoticCI(const arma::mat& X, const arma::vec& res, const arma::vec& coeff, const double tau, 
                       const int n, const double h, const double z) {
  arma::mat Z = X;
  arma::vec Wh = arma::square(arma::normcdf(-res / h) - tau);
  arma::vec Wh2 = arma::normpdf(res / h) / h;
  arma::mat temp = Z;
  for (int i = 0; i < n; i++) {
    temp.row(i) = Wh.row(i) * Z.row(i);
  }
  arma::mat Stau = temp.t() * Z;
  for (int i = 0; i < n; i++) {
    temp.row(i) = Wh2.row(i) * Z.row(i);
  }
  arma::mat Dh = temp.t() * Z;
  arma::mat Dhinv = Dh.i();
  arma::vec tm = z * arma::sqrt(arma::diagvec(Dhinv * Stau * Dhinv));
  return arma::join_rows(coeff - tm, coeff + tm);
}











