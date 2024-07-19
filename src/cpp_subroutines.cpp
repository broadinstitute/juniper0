#include <Rcpp.h>
using namespace Rcpp;

// Code taken from Didelot et al. (2017)

// [[Rcpp::export]]
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t)
{
  int n = std::round((dateT-tinf)/delta_t);
  NumericVector grid(n);
  for(int i=0; i<n; ++i) // use the left point of each subinterval
    grid[i] = dateT-n*delta_t+i*delta_t;

  NumericVector pi2 = pi*pgamma(dateT-grid, shSam, scSam);
  NumericVector F = 1-pgamma(dateT-grid, shGen, scGen);

  NumericVector w(n), out(n);

  IntegerVector seq = seq_len(n);
  NumericVector gam = dgamma(as<NumericVector>(seq)*delta_t,shGen,scGen);
  double sumPrev = 0.5 * gam[0];
  out[n-1]=std::min(1.0,F[n-1]+sumPrev*delta_t);
  for(int i=n-1; i>0; --i){
    w[i] = (1-pi2[i]) * pow((1-pOff)/(1-pOff*out[i]), rOff);

    sumPrev = 0.0;
    for(int j=0; j<n-i; ++j)
      sumPrev += gam[j]*w[i+j];
    sumPrev += 0.5 * gam[n-i];
    out[i-1] = std::min(1.0,F[i-1] + sumPrev*delta_t);
  }
  return log(out);
}
