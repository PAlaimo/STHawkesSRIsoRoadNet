#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// // [[Rcpp::export]]
// NumericVector kernelWithEdgeCorr_trend(NumericVector x, NumericVector mus, double sigma, double endt){
//   int n_len = x.size();
//   NumericVector oo(n_len);
//   for(int i = 0; i < n_len; i++){
//     double num = R::dnorm4(x[i], mus[i], sigma, 0) + R::dnorm4(x[i], -mus[i], sigma, 0) + R::dnorm4(x[i], 2*endt-mus[i], sigma, 0);
//     double den = R::pnorm5(2*endt, mus[i], sigma, 0, 0) -  R::pnorm5(-endt, mus[i], sigma, 0, 0);
//     oo[i] = num/den;
//   }
//   return oo;
// }

// [[Rcpp::export]]
NumericMatrix kernelWithEdgeCorr_trend(NumericVector x, NumericVector mus, double sigma, double startt, double endt){
  int n_len = x.size();
  int mu_len = mus.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < n_len; i++){
    for(int j = 0; j < mu_len; j++){
      double num = R::dnorm4(x[i], mus[j], sigma, 0);
      double den = R::pnorm5(endt, mus[j], sigma, 1, 0) -  R::pnorm5(startt, mus[j], sigma, 1, 0);
      oo(j,i) = num/den; 
    }
  }
  return oo;
}

// // [[Rcpp::export]]
// NumericVector kernelWithEdgeCorr_daily(NumericVector x, NumericVector mus, double sigma){
//   int n_len = x.size();
//   NumericVector oo(n_len);
//   for(int i = 0; i < n_len; i++){
//     oo[i] = R::dnorm4(x[i], mus[i], sigma, 0) + R::dnorm4(x[i], -1+mus[i], sigma, 0) + R::dnorm4(x[i], 1+mus[i], sigma, 0);
//   }
//   return oo;
// }
// 
// // [[Rcpp::export]]
// NumericVector kernelWithEdgeCorr_weekly(NumericVector x, NumericVector mus, double sigma){
//   int n_len = x.size();
//   NumericVector oo(n_len);
//   for(int i = 0; i < n_len; i++){
//     oo[i] = R::dnorm4(x[i], mus[i], sigma, 0) + R::dnorm4(x[i], -1+mus[i], sigma, 0) + R::dnorm4(x[i], 1+mus[i], sigma, 0);
//   }
//   return oo;
// }

// [[Rcpp::export]]
NumericMatrix kernel_periodic(NumericVector x, NumericVector mus, double sigma, int c){
  int n_len = x.size();
  int mu_len = mus.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < n_len; i++){
    for(int j = 0; j < mu_len; j++){
      oo(j,i) = R::dnorm4(x[i], mus[j], sigma, 0) + R::dnorm4(x[i], -c+mus[j], sigma, 0) + R::dnorm4(x[i], c+mus[j], sigma, 0);
    }
  }
  return oo;
}

// // [[Rcpp::export]]
// NumericVector kernelWithEdgeCorr_background(NumericVector x, NumericVector y, NumericVector musx, NumericVector musy, NumericVector sigmas){
//   int n_len = x.size();
//   NumericVector oo(n_len);
//   for(int i = 0; i < n_len; i++){
//     oo[i] = R::dnorm4(x[i], musx[i], sigmas[i], 0) * R::dnorm4(y[i], musy[i], sigmas[i], 0);
//   }
//   return oo;
// }

// [[Rcpp::export]]
NumericMatrix kernelWithEdgeCorr_background(NumericVector x, NumericVector y, NumericVector musx, NumericVector musy, NumericVector sigmas, NumericVector integ){
  int x_len = x.size();
  // int y_len = y.size();
  int mux_len = musx.size();
  // int muy_len = musy.size();
  NumericMatrix oo(mux_len, x_len);
  for(int i = 0; i < mux_len; i++){
    for(int j = 0; j < x_len; j++){
      oo(i,j) = (R::dnorm4(x[j], musx[i], sigmas[i], 0) * R::dnorm4(y[j], musy[i], sigmas[i], 0))/integ[i];
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericMatrix kernelWithEdgeCorr_backgroundx(NumericVector x, NumericVector musx, NumericVector sigmas, NumericVector integ){
  int x_len = x.size();
  // int y_len = y.size();
  int mux_len = musx.size();
  // int muy_len = musy.size();
  NumericMatrix oo(mux_len, x_len);
  for(int i = 0; i < mux_len; i++){
    for(int j = 0; j < x_len; j++){
      oo(i,j) = (R::dnorm4(x[j], musx[i], sigmas[i], 0)/sqrt(integ[i]));
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericMatrix kernelWithEdgeCorr_backgroundy(NumericVector y, NumericVector musy, NumericVector sigmas, NumericVector integ){
  // int x_len = x.size();
  int y_len = y.size();
  // int mux_len = musx.size();
  int muy_len = musy.size();
  NumericMatrix oo(muy_len, y_len);
  for(int i = 0; i < muy_len; i++){
    for(int j = 0; j < y_len; j++){
      oo(i,j) = (R::dnorm4(y[j], musy[i], sigmas[i], 0)/sqrt(integ[i]));
    }
  }
  return oo;
}


// // [[Rcpp::export]]
// NumericVector smooth_background(NumericMatrix kernx, NumericMatrix kerny, NumericVector phi){
//   int x_len = kernx.ncol();
//   int y_len = kerny.ncol();
//   int p_len = phi.size();
//   // int muy_len = musy.size();
//   NumericVector oo(x_len*y_len);
//   for(int i = 0; i < x_len; i++){
//     for(int j = 0; j < y_len; j++){
//       for(int k = 0; k < p_len; k++){
//         oo[i+j*x_len] += phi[k]*kernx(k, i)*kerny(k, j);
//       }
//     }
//   }
//   return oo;
// }

// [[Rcpp::export]]
NumericVector smooth_background(NumericMatrix kernx, NumericMatrix kerny, NumericVector phi, 
                                NumericMatrix idgrid, int gridsize){
  // int x_len = kernx.ncol();
  // int y_len = kerny.ncol();
  int p_len = phi.size();
  // int muy_len = musy.size();
  NumericVector oo(gridsize);
  for(int i = 0; i < gridsize; i++){
    int xid = idgrid(i, 0);
    int yid = idgrid(i, 1);
    for(int k = 0; k < p_len; k++){
      oo[i] += phi[k]*kernx(k, xid)*kerny(k, yid);
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericVector smooth_background_sparse(const sp_mat& kernx, const sp_mat& kerny, NumericVector phi, 
                                NumericMatrix idgrid, int gridsize){
  // int x_len = kernx.ncol();
  // int y_len = kerny.ncol();
  int p_len = phi.size();
  // int muy_len = musy.size();
  NumericVector oo(gridsize);
  for(int i = 0; i < gridsize; i++){
    int xid = idgrid(i, 0);
    int yid = idgrid(i, 1);
    arma::sp_mat kernxi(kernx.col(xid));
    arma::sp_mat kernyi(kerny.col(yid));
    // for(int k = 0; k < p_len; k++){
    //   oo[i] += phi[k]*kernxi(k)*kernyi(k);
    // }
    for(arma::sp_mat::const_iterator k = kernxi.begin(); k != kernxi.end(); ++k){
      oo[i] += *k * phi[k.row()]*kernyi(k.row());
    }
  }
  return oo;
}


// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_tempRes(NumericVector x, NumericVector musx, double sigma, double endt){
  int n_len = x.size();
  int mu_len = musx.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < mu_len; i++){
    for(int j = 0; j < n_len; j++){
      double num = R::dnorm4(x[j], musx[i], sigma, 0) + R::dnorm4(x[j], -musx[i], sigma, 0);
      double den = R::pnorm5(endt, musx[i], sigma, 1, 0) - R::pnorm5(-endt, musx[i], sigma, 1, 0);
      oo(i,j) = num/den; 
    }
  }
  return oo;
}


// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_spatRes(NumericVector x, NumericVector y, NumericVector musx, NumericVector musy, double sigma){
  int n_len = x.size();
  int mu_len = musx.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < mu_len; i++){
    for(int j = 0; j < n_len; j++){
      double num = R::dnorm4(x[j], musx[i], sigma, 0) * R::dnorm4(y[j], musy[i], sigma, 0);
      double den = (R::pnorm5(1, musx[i], sigma, 1, 0) - R::pnorm5(-1, musx[i], sigma, 1, 0))*(R::pnorm5(1, musy[i], sigma, 1, 0) - R::pnorm5(-1, musy[i], sigma, 1, 0));
      oo(i,j) = num/den; 
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_spatResISO(NumericVector x, NumericVector musx, double sigma, double ends){
  int n_len = x.size();
  int mu_len = musx.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < mu_len; i++){
    for(int j = 0; j < n_len; j++){
      double num = (R::dnorm4(x[j], musx[i], sigma, 0) - R::dnorm4(x[j], -musx[i], sigma, 0))/(2*M_PI*x[j]);
      double den = R::pnorm5(ends, musx[i], sigma, 1, 0) - R::pnorm5(-ends, musx[i], sigma, 1, 0);
      oo(i,j) = num/den; 
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_spatResISOrepNeg(NumericVector x, NumericVector musx, double sigma, double ends){
  int n_len = x.size();
  int mu_len = musx.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < mu_len; i++){
    for(int j = 0; j < n_len; j++){
      double num = (R::dnorm4(x[j], musx[i], sigma, 0) - R::dnorm4(x[j], -musx[i], sigma, 0));
      double den = R::pnorm5(ends, musx[i], sigma, 1, 0) - R::pnorm5(-ends, musx[i], sigma, 1, 0);
      oo(i,j) = num/den; 
    }
  }
  return oo;
}

// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_spatResISOrep(NumericVector x, NumericVector musx, double sigma, double ends){
  int n_len = x.size();
  int mu_len = musx.size();
  NumericMatrix oo(mu_len, n_len);
  for(int i = 0; i < mu_len; i++){
    for(int j = 0; j < n_len; j++){
      double num = (R::dnorm4(x[j], musx[i], sigma, 0));
      double den = R::pnorm5(ends, musx[i], sigma, 1, 0) - R::pnorm5(0, musx[i], sigma, 1, 0);
      oo(i,j) = num/den; 
    }
  }
  return oo;
}


// [[Rcpp::export]]
NumericVector bwCalc_cpp(NumericVector x, NumericVector y, int nmin){
  int n_len = x.size();
  NumericVector oo(n_len);
  NumericVector d_nmin(n_len);
  for(int i = 0; i < n_len; i++){
    for(int j = 0; j < n_len; j++){
      oo[j] = sqrt((x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]));
    } 
    d_nmin[i] = oo.sort()[nmin-1];
  }
  
  return d_nmin;
}

// [[Rcpp::export]]
NumericMatrix kerneledgeCorr_spatResx(NumericVector x, NumericVector musx, double sigma, double end){
  int x_len = x.size();
  int mux_len = musx.size();
  NumericMatrix oo(mux_len, x_len);
  for(int i = 0; i < mux_len; i++){
    double den = (R::pnorm5(end, musx[i], sigma, 1, 0) - R::pnorm5(-end, musx[i], sigma, 1, 0));
    for(int j = 0; j < x_len; j++){
      double num = R::dnorm4(x[j], musx[i], sigma, 0);
      oo(i,j) = num/den; 
    }
  }
  return oo;
}


