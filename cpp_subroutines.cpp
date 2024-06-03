#include <Rcpp.h>
#include <cmath>
#include <iostream>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double epi(Rcpp::NumericVector q_all, double rho, double psi, bool obs, int N, int w, int M, double tol){
  double out = 0;

  std::vector<double> q;
  std::vector<double> qsmall;

  int nsmall = 0;
  int nbig = 0;

  for(int j = 0; j < q_all.length(); ++j){
    if(q_all[j] <= log(tol)){
      qsmall.push_back(q_all[j]);
      nsmall = nsmall + 1;
    }else{
      q.push_back(exp(q_all[j]));
      nbig = nbig + 1;
    }
  }

  std::vector<double> qtmp(nbig, 1.0);

  double lcoef = rho*log(psi) + log(1-psi) + log(rho);
  double lterm = lcoef;

  for(int j = 0; j < nbig; ++j){
    qtmp[j] = qtmp[j] * (1 - q[j]);
    lterm = lterm + log(1 - qtmp[j]);
  }
  out = lterm;
  double ratio = 0;

  for(int i = 2; i <= M; ++i) {
    lcoef = lcoef + log(1-psi) + log(i + rho - 1) - log(i);
    lterm = lcoef;
    for(int j = 0; j < nbig; ++j){
      qtmp[j] = qtmp[j] * (1 - q[j]);
      lterm = lterm + log(1 - qtmp[j]);
    }
    lterm = lterm + nsmall*log(i);

    ratio = log(1 + (exp(lterm - out)));
    out = out + ratio;
  }

  for(int j = 0; j < nsmall; ++j){
    out = out + qsmall[j];
  }

  if(obs){
    out = out - log(N);
  }

  if(w > 0){
    out = epi(Rcpp::NumericVector (1, out), rho, psi, FALSE, N, w-1, M, tol);
  }

  return out;
}

// [[Rcpp::export]]
std::vector<double> get_qs(int n, int N, int G, double rho, double psi){
  std::vector<double> q;
  double qtmp = log(1 - (double)n/N);
  q.push_back(qtmp);

  for(int i = 1; i < G; ++i){
    q.insert(q.begin(), (
      qtmp + rho * (log(psi) - log(1 - (1-psi) * exp(q[0])))
    ));
  }
  q.push_back(0);
  return q;
}

// [[Rcpp::export]]
double e_lik_node(int d, int w, std::vector<double> qs, double rho, double psi){
  double out = 0;
  for(int i=0; i < qs.size() - 1; ++i){
    out = out + log(1-psi) + rho*log(psi) - (rho+1)*log(1 + exp(qs[i])*(psi - 1)) + log(rho);
  }
  out = out + d*log(1-psi) + rho*log(psi) - (rho+d)*log(1 + exp(qs.back())*(psi - 1)) + std::lgamma(rho + d) - std::lgamma(rho) - std::lgamma(d + 1);

  return out;
}

// [[Rcpp::export]]
std::vector<int> gen(std::vector<int> h, std::vector<int> w){
  int n = h.size();
  std::vector<int> out(n, 0);
  std::deque<int> queue = {1};
  while(queue.size() > 0){
    if(queue[0] != 1){
      out[queue[0] - 1] = w[queue[0] - 1] + out[h[queue[0] - 1] - 1] + 1;
    }
    for(int j=0; j<n; ++j){
      if(h[j] == queue[0]){
        queue.push_back(j+1);
      }
    }
    queue.pop_front();
  }
  return(out);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
