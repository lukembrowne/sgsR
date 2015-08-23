#include <Rcpp.h>
using namespace Rcpp;

// ***************************
// ***************************
// ** Calculate pairwise distance between all individuals
// ***************************
// ***************************

// [[Rcpp::export]]
NumericMatrix calcPairwiseDist(NumericVector x,
                               NumericVector y,
                               int Nind) {

  double dij; // Distance between individuals i and j
  NumericMatrix Mdij(Nind, Nind); // Distance matrix of individuals

  for(int i = 0; i < (Nind - 1); i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

   // Rcpp:Rcout << i << " - " << j << "\n";
    dij = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])); // Calculate distance
    Mdij(i, j) = dij; // Assign pariwise distance into matrix

  } // End loop through pairs

  return(Mdij); // Return distance matrix
}





// ***************************
// ***************************
// ** Place pairwise comparisions into distance classes
// ***************************
// ***************************

// Need to add ability to deal with distance classes bigger than max dist int given
// Looping through the pairwise matrix for each distance interval is likely not the most efficient way to do this

// [[Rcpp::export]]
NumericMatrix findDIs(NumericMatrix Mdij,
                      NumericVector DistInts,
                      int Nind) {

    NumericMatrix Mcij(Nind, Nind); // Matrix saying which distance class each pairwise combo is in

    for(int di = 0; di < DistInts.size(); di++){
      for(int i = 0; i < (Nind - 1); i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

        if(di == 0){ // Special loop for first interval because there's no lower interval
          if(Mdij(i, j) < DistInts[di]) Mcij(i, j) = di;
        }
        if(di > 0){ // Find interval where the pairwise distance fits inside
          if((Mdij(i, j) > DistInts[di - 1]) & (Mdij(i, j) < DistInts[di])) Mcij(i, j) = di;
        }
      } // End loop through pairs
    } // End loop through distance interval
  return(Mcij); // Return distance matrix
}




// ***************************
// ***************************
// ** Place distance intervals with equal pairwise comparisons
// ***************************
// ***************************

// Need to add ability to deal with distance classes bigger than max dist int given
// Looping through the pairwise matrix for each distance interval is likely not the most efficient way to do this

// [[Rcpp::export]]
NumericVector findEqualDIs(NumericMatrix Mdij,
                      NumericVector DistInts,
                      int Nind) {
  int ctr = 0; // Counter
  int npairs_per_c = 0; // Number of pairs per class

    // Find out how many pairs there are to flatten out distance matrix
      for(int i = Nind; i > 0; i--){
        ctr += i;
      }
      ctr = ctr - Nind; // Subtract diagonal comparisons
      NumericVector flatdij(ctr); // Initialize flattened Mdij

      int k = 0;
      for(int i = 0; i < (Nind - 1); i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals
       flatdij[k] = Mdij(i, j); // Flatten Mdij into a vector
        k += 1;
       } // End loop through pairs

      std::sort(flatdij.begin(), flatdij.end()); // Sort from smallest to greatest distance

    // Find distance intervals to cut

    // Calculate number of pairs per interval
    int ndis = -DistInts[0]; // Take negative of negative

    npairs_per_c = flatdij.size() / ndis;

    // Save new distance intervals
    NumericVector DistIntsEqual(ndis); // Reinitialize distance interval vector

    int index = npairs_per_c - 1; // Minus 1 for zero indexing, or else it messes up with factors
    for(int i = 0; i < ndis; i++){
      DistIntsEqual[i] = flatdij[index] * 1.001; // Make sure distance is a little bit higher by multiplying by 1.001
      index += npairs_per_c;
    }
    // Make sure max distance is included
    DistIntsEqual[ndis-1] = max(flatdij) * 1.001;

 return(DistIntsEqual);
}


// ***************************
// ***************************
// ** Summarize information for distance intervals
// ***************************
// ***************************

// Add CV and & participation information

// [[Rcpp::export]]
NumericMatrix summarizeDIs(NumericMatrix Mdij,
                           NumericMatrix Mcij,
                           NumericVector DistInts,
                           int Nind){

  int ndis = DistInts.size();
  int npairs[ndis];
  float tot[ndis]; // Total distance
  NumericMatrix summary(4, ndis); // nrows determined by how many stats to summarize

  // Initalize pairs and total variables
  for(int i = 0; i < ndis; i++){
    npairs[i] = 0;
    tot[i] = 0;
  }

  for(int i = 0; i < (Nind - 1); i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

    int di = Mcij(i, j); // Distance interval for this pair
    tot[di] += Mdij(i, j); // Add distance to total distance
    npairs[di] += 1;      // Increment number of pairs by one
    } // End loop through pairs

  /*
  Fill in summary information
  Row order is ...
  1. Name of distance interval
  2. Max distance in distance interval
  3. Average distance in distance itnerval
  4. Number of pairs
  */

  for(int i = 0; i < ndis; i++){

    summary(0, i) = i; // Name of distance interval
    summary(1, i) = DistInts[i]; // Max distance in distance interval
    summary(2, i) = tot[i] / npairs[i]; // Average distance
    summary(3, i) = npairs[i]; // Number of pairs
  }

  return(summary);

} // End function





