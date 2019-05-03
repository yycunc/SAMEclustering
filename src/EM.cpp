#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

//Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector rdirichletC (int n){
  arma::vec v = arma::randg<arma::vec>(n);
  double sum = arma::sum(v);
  arma::vec out = v/sum;
  return Rcpp::as<Rcpp::NumericVector>(wrap(out));
}

// [[Rcpp::export]]
std::map<int, int> tableC (Rcpp::NumericVector x) {
  std::map<int, int> counts;
  
  int n = x.size();
  for (int i = 0; i < n; i++){
    if(!Rcpp::NumericVector::is_na(x[i])){
      counts[x[i]]++;
    }
  }
  
  return counts;
}

// [[Rcpp::export]]
Rcpp::List EM (const Rcpp::NumericMatrix & data, const int M){
  int N = data.nrow();  //number of subjects
  int H = data.ncol();  //number of clustering methods
  Rcpp::IntegerVector K(H);
  for(int j=0; j<H; j++){
    K[j] = max(na_omit(data(_,j)));
  }
  
  //start with random weight
  Rcpp::NumericVector alpha = rdirichletC(M);
  //creat M lists of nu
  Rcpp::List nu(M);
  for(int m=0; m<M; m++){
    Rcpp::List L(H);
    for(int j=0; j<H; j++){
      std::map<int,int> count = tableC(data(_,j));
      Rcpp::NumericVector nuj(count.size());
      nuj = rdirichletC(count.size());
      L[j] = nuj;
    }
    nu[m] = L;
  }

  //initialize Z vector
  Rcpp::NumericMatrix Z(N,M);
  std::fill(Z.begin(), Z.end(), 0);
  
  //start EM algorithm
  Rcpp::NumericVector loglik = Rcpp::NumericVector::create(1.0); 
  double diff = 1.0;
  
  while(diff > 0.0001){
    for(int i=0; i<N; i++){
      Rcpp::NumericVector prod(M);
      for(int m=0; m<M; m++){
       double prod_m = 1.0;
        Rcpp::List num = nu[m];
        for(int j=0; j<H; j++){
          if(!Rcpp::NumericVector::is_na(data(i,j))){
            prod_m = prod_m * Rcpp::as<Rcpp::NumericVector>(num[j])[data(i,j)-1];
          }
        }
        Z(i,m) = alpha[m] * prod_m;
        prod[m] = Z(i,m);
      }
      Z(i,_) = Z(i,_)/sum(prod);
    }
  
  //maximization step Equation 17
  //maximize alpha
  for(int m=0; m<M; m++){
    alpha[m] = sum(Z(_,m))/sum(Z);
  }
  
  //maximize nus Equation 19
  for(int m=0; m<M; m++){
    Rcpp::List updatenu = nu[m];
    Rcpp::List L(H);
      for(int j=0; j<H; j++){
        Rcpp::NumericVector nuj(H);
        for(int k=0; k<K[j]; k++){
          Rcpp::NumericVector numerator(K[j]);
          std::fill(numerator.begin(), numerator.end(), 0);
          for(int i=0; i<N; i++){
            if(!Rcpp::NumericVector::is_na(data(i,j))){
              if(k==data(i,j)-1){
                numerator[k] += Z(i,m);
              }
            } else if(Rcpp::NumericVector::is_na(data(i,j))){
              //calculate E[z_im,data_i^mis|data_i^obs, theta] for missing data
              Rcpp::List num = nu[m];
              numerator[k] += Rcpp::as<Rcpp::NumericVector>(num[j])[k] * Z(i,m);
            }
          }
          Rcpp::as<Rcpp::NumericVector>(updatenu[j])[k] = numerator[k];
        }
        nuj = Rcpp::as<Rcpp::NumericVector>(updatenu[j])/sum(Rcpp::as<Rcpp::NumericVector>(updatenu[j]));
        L[j] = nuj;
      }
      nu[m] = L;
    }
    
    //likelihood
    double loglikelihood = 0;
    for(int i=0; i<N; i++){
      Rcpp::NumericVector prod(M);
      for(int m=0; m<M; m++){
        double prod_m = 1.0;
        Rcpp::List num = nu[m];
        for(int j=0; j<H; j++){
          if(!Rcpp::NumericVector::is_na(data(i,j))){
            prod_m *= Rcpp::as<Rcpp::NumericVector>(num[j])[data(i,j)-1];
          }
        }
        prod[m] = alpha[m] * prod_m;
      }
      loglikelihood += log(sum(prod));
    }
    loglik.push_back(loglikelihood);
    if(loglik.size() == 1){
      diff = loglik[0];
    } else{
      diff = std::abs(loglik[loglik.size()-1]-loglik[loglik.size()-2]);
    }
  }

  Rcpp::NumericVector clusterresults(N);
  for(int i=0; i<N; i++){
    double j = which_max(Z(i,_));
    clusterresults[i] = j + 1.0;
  }
  std::map<int, int> ensemblecluster = tableC(clusterresults);
  Rcpp::NumericVector clusters = unique(clusterresults);
  int k = clusters.size();
  
  //number of parameters
  std::size_t nu0_size = Rcpp::as<Rcpp::List>(nu[0]).size();
  std::size_t total_length = 0;
  for(std::size_t i=0; i<nu0_size; i++){
    total_length += Rf_length(Rcpp::as<Rcpp::List>(nu[0])[i]);
  }
  double p = alpha.size() -1 + (total_length - H) * M;
  
  //AIC and BIC
  double AIC = -2 * loglik(loglik.size() - 1) + 2 * p;
  double BIC = -2 * loglik(loglik.size() - 1) + log(N) * p;
  double LL = loglik(loglik.size() - 1);

  return Rcpp::List::create(Rcpp::Named("Number of clusters") = k,
                            Rcpp::Named("Cluster results") = clusterresults,
                            Rcpp::Named("Ensemble cluster") = ensemblecluster,
                            Rcpp::Named("AIC") = AIC,
                            Rcpp::Named("BIC") = BIC,
                            Rcpp::Named("Number of parameter") = p,
                            Rcpp::Named("Log Likelihood") = LL,
                            Rcpp::Named("Alpha") = alpha);

}

