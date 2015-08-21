#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericMatrix calcPairwiseDist(NumericVector x,
                               NumericVector y) {

  int n = x.size(); // Number of individuals
  double dij; // Distance between individuals i and j
  NumericMatrix Mdij(n, n); // Distance matrix of individuals

  for(int i = 0; i < (n - 1) ; i++) for(int j = i+1; j <= (n - 1); j++){ // Loops through all pairs of individuals

   // Rcpp:Rcout << i << " - " << j << "\n";

    dij = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])); // Calculate distance

    Mdij(i, j) = dij; // Assign pariwise distance into matrix

  } // End loop through pairs

  return(Mdij); // Return distance matrix
}


