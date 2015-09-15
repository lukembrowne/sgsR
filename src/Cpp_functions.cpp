#include <Rcpp.h>
using namespace Rcpp;


// Declare functions
std::vector<float> calcFijPairwiseCpp(List ref_gen,
                         NumericMatrix alfreq1,
                         NumericMatrix alfreq2,
                         int Nloci, NumericVector Nallele,
                         int n);

NumericVector calcAlleleFreqPop(NumericVector alleles_1,
                                NumericVector alleles_2,
                                int Nallele);


std::vector<float> calcAlleleFreqCppInd(int alleles_1,
                                        int alleles_2,
                                        int Nallele);

NumericMatrix findDIs(NumericMatrix Mdij,
                      NumericVector distance_intervals,
                      int Nind);

NumericMatrix calcPairwiseDist(NumericVector x,
                               NumericVector y,
                               int Nind);



// ***************************************************
// ***************************************************
// *** CALCULATE PAIRWISE FIJ BETWEEN A PAIR OF INDIVIDUALS
// ***************************************************
// ***************************************************

// Returns a vector with Fij estimate at each locus

// [[Rcpp::export]]
std::vector<float> calcFijPairwiseCpp(List ref_gen,
                         NumericMatrix alfreq1,
                         NumericMatrix alfreq2,
                         int Nloci, NumericVector Nallele,
                         int n){

  float denom[Nloci];
  float numer[Nloci];
  std::vector<float> fij(Nloci);
  NumericVector ref_gen_loc;

  // Initialize numerator and denominator
  for(int locus = 0 ; locus < Nloci; ++locus){
    denom[locus] = 0;
    numer[locus] = 0;
  }

  for(int locus = 0;  locus < Nloci; ++locus){ // Loop through loci

    for(int allele = 0; allele < Nallele[locus]; ++allele){ // Loop through alleles

    ref_gen_loc = ref_gen[locus]; // Subset ref gen list to just row we're interested in

      // Calculate Numerator and Denominator of Loiselle et al. 1995
      numer[locus]  +=  (alfreq1(locus, allele) - ref_gen_loc[allele]) *
                        (alfreq2(locus, allele) - ref_gen_loc[allele]) +
        (ref_gen_loc[allele]*(1 - ref_gen_loc[allele])) / (n - 1);

      denom[locus] += ref_gen_loc[allele] * (1 - ref_gen_loc[allele]);

    } // End allele loop

    fij[locus] = numer[locus] / denom[locus]; // Calculate Fij per locus

  } // End loci loop

  return(fij);
}



// ***************************************************
// ***************************************************
// ****CALCULATE PAIRWISE FIJ AMONG ALL INDIVIDUALS IN A POPULATION
// ***************************************************
// ***************************************************

// Returns a matrix with first section as Fij estimates for each locus individiually, then averaged across loci
// Second section is average estimates after permutation of spatial location among individuals
// Third and fourth section are the 2.5 % and 97.5 % quantiles of the permuted data

// [[Rcpp::export]]
NumericMatrix calcFijPopCpp(NumericMatrix Mcij,
                            NumericMatrix Mdij,
                            NumericVector distance_intervals,
                            NumericMatrix genotype_data,
                            List ref_gen,
                            int Nloci,
                            NumericVector Nallele,
                            int Nind,
                            int Ngenecopies,
                            int Nperm,
                            NumericVector x_coord,
                            NumericVector y_coord){

  // Initialize variables
  int MNallele = max(Nallele);
  float ind_al_freq[Nloci][MNallele][Nind]; // 3d array that saves allele freq for each individual
  int ndis = distance_intervals.size(); // Number of distance intervals
  float perm_results[Nloci + 1][ndis][Nperm + 1]; // 3d array that saves permuted results
  int perm_results_npairs[Nloci + 1][ndis][Nperm + 1]; // 3d array that saves permuted results
  int row = 0; // Needed for loop through 3d array
  std::vector<float> fij(Nloci);
  NumericMatrix  alfreq1(Nloci, MNallele);
  NumericMatrix  alfreq2(Nloci, MNallele);
  NumericMatrix fijsummary((Nloci + 1) * 4, ndis); // Stores summary information of pairwise Fij
  NumericVector quant975(Nperm);
  NumericVector quant025(Nperm);


  // Initialize matrices - fill them with 0s, or could cause problems later on
  for(int locus = 0; locus < (Nloci + 1); locus++) for(int di = 0; di < ndis; di++) {

    fijsummary(locus, di) = 0;
    fijsummary(locus + Nloci, di) = 0;

     for(int perm = 0; perm < (Nperm + 1); perm++){
      perm_results[locus][di][perm] = 0;
      perm_results_npairs[locus][di][perm] = 0;
     }
   }

  // Loop through individuals and calculate allele frequencies for each individual
  // Frequency data looks something like: 0 0 0 .5 .5 0 0 (for a heterozygote)
  // Save results in a 3d array that is arraged [locus][allele][individual]
  for(int indi = 0; indi < Nind; ++indi){

    row = 0; // Row of individual allele frequency 3d array - corresponds to locus
    for(int locus = 0; locus < (Nloci * 2); locus += 2){ // Loop through loci


      std::vector<float> freq_table(Nallele[row]);

      freq_table = calcAlleleFreqCppInd(genotype_data(indi, locus),
                                                         genotype_data(indi, locus + 1),
                                                         Nallele[row]);
      /*
      // Code to print out individual frequency tables
      std::cout << "Indi: " << indi <<" | Locus - "<< locus << "|:";
      for (std::vector<float>::iterator i = freq_table.begin(); i != freq_table.end(); ++i)
      std::cout << *i << ' ';
      Rcpp::Rcout << "\n";
      */

      // Have to then loop through allele frequency table to fill in individual part
      // of the array
      int allele = 0;
      for (std::vector<float>::iterator i = freq_table.begin(); i != freq_table.end(); ++i){

        ind_al_freq[row][allele][indi] = *i;
        allele += 1;
      } // End allele loop


      row += 1; // Jump to next row in array, which is the next locus

    } // End loci loop
  } // End individuals loop
  // END CALCULATING INDIVIDUAL ALLELE FREQUENCIES

  // PERMUTATION AND CALCULATING FIJ

  // Initialize variables
  int spat_size = x_coord.size(); // Save length of spatial coordinates
  NumericVector loc_index(spat_size); // Create vector of indices that will be shuffled later
  for(int i = 0; i < spat_size; i++) loc_index[i] = i; // Sequence of integers from 0 to ...

  NumericVector x_coord_shuf(spat_size); // Initialize vectors that will hold shuffled locations
  NumericVector y_coord_shuf(spat_size);



  // Start permutations of individuals among locations
  // Must ensure that x and y coordinates are shuffled in unison so that completely new spatial locations are not created during the permutation process

  // Note that the OBSERVED data are held in perm_results[locus][distance_interval][0]

  for(int perm = 0; perm < (Nperm + 1); perm++){ // Loop through observed (perm = 0) and then total number of permutations

    // Print status update on how far along we are in the permutations..
    if(perm % 100 == 0 ){
      Rcpp::Rcout << "Working on permutation: " << perm << "... \n";
    }

      // Randomizing spatial locations among individuals
      if(perm > 0){ // Skip if it's the first permutation and use observed Mdij and Mcij

        // Otherwise, shuffle location indices and corresponding spatial coordinates
        // This function alters the object in place, so be careful!
        std::random_shuffle(loc_index.begin(), loc_index.end()); // Shuffle location indices

        for(int i = 0; i < spat_size; i++){
          x_coord_shuf[i] = x_coord[loc_index[i]];
          y_coord_shuf[i] = y_coord[loc_index[i]];
        }

        // Recalculate distance between individuals and distance interval classification
        Mdij = calcPairwiseDist(x_coord_shuf, y_coord_shuf, Nind);
        Mcij = findDIs(Mdij, distance_intervals, Nind);

      } // End shuffling spatial locations




      // CALCULATING PAIRWISE Fij

      // Calculate pairwise Fij between all pairs of individuals
      // Loop through all pairs of individuals
      for(int i = 0; i < (Nind - 1) ; i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

        int di = Mcij(i, j); // Save distance class as an index to save in perm_results

        // Loop through and save allele frequency separately for each individual
        for(int locus = 0; locus < Nloci; ++locus){
          for(int allele = 0; allele < Nallele[locus]; ++allele){
            alfreq1(locus, allele) = ind_al_freq[locus][allele][i];
            alfreq2(locus, allele) = ind_al_freq[locus][allele][j];
          }
        }

        // Should return an array / vector with Fij estimate for each locus
        fij  =   calcFijPairwiseCpp(ref_gen, alfreq1, alfreq2,
                                         Nloci, Nallele, Ngenecopies);

        // Save per locus Fij estimate
        int locus = 0; // Initialize before loop
        for (std::vector<float>::iterator k = fij.begin(); k != fij.end(); ++k){

          // First index of perm results is the OBSERVED
          perm_results[locus][di][perm] += *k; // *k accesses value in fij vector
          perm_results[Nloci][di][perm] += *k; // Save sum across all loci in last row

          perm_results_npairs[locus][di][perm] += 1; // Save number of pairwise combinations to use in calculate averages
          perm_results_npairs[Nloci][di][perm] += 1;

          locus += 1;

        } // End locus loop
      } // End loop through pairs of individuals

    } // End permutation loop


  //Calculate average Fij per locus by dividing sum by total pairs
  // Note that locus <= Nloci in loop to calculate for average across loci as well
   for(int locus = 0; locus <= Nloci; locus++) for(int di = 0; di < ndis; di++) for(int perm = 0; perm < (Nperm + 1); perm++){
      perm_results[locus][di][perm] = perm_results[locus][di][perm] / perm_results_npairs[locus][di][perm];
   }

   // Save out to summary matrix
   for(int locus = 0; locus <= Nloci; locus++) for(int di = 0; di < ndis; di++){

     // Observed values - saved in first index of perm_results
     fijsummary(locus, di) = perm_results[locus][di][0];

     for(int perm = 1; perm < (Nperm + 1); perm++){ // Start perm at 1 to not overwrite observed values

       fijsummary(locus + Nloci + 1, di) += perm_results[locus][di][perm];
       fijsummary(Nloci + Nloci + 1, di) += perm_results[Nloci][di][perm]; // Sum across loci

       // Save values to calculate 95% quantiles
       quant975[perm-1] = perm_results[locus][di][perm]; // -1 index bc perm is looping starting on perm = 1
       quant025[perm-1] = perm_results[locus][di][perm];

     } // End permutation loop

     // Calculate 2.5 % and 95 % quantiles
     std::sort(quant975.begin(), quant975.end()); // Sort vectors
     std::sort(quant025.begin(), quant025.end());
     float quant975_val = quant975[quant975.size()*(0.975 - 0.000000001)]; // Assign quantile value
     float quant025_val = quant025[quant025.size()*(0.025 - 0.000000001)];


     fijsummary(locus + Nloci* 2 + 2, di) = quant025_val; // Save in fijsummary after observed and perm values
     fijsummary(locus + Nloci* 3 + 3, di) = quant975_val;

   } // End locus and di loop

   for(int locus = 0; locus <= Nloci; locus++) for(int di = 0; di < ndis; di++){
     fijsummary(locus + Nloci + 1, di) = fijsummary(locus + Nloci + 1, di) / Nperm; // Take avg based on permutation
     fijsummary(Nloci + Nloci + 1, di) = fijsummary(Nloci + Nloci + 1, di) / Nperm;
   }

   // Save out summary information of permuted results
   return(fijsummary);

}



// ***************************************************
// ***************************************************
// **** CALCULATE ALLELE FREQUENCY OF A POPULATION
// ***************************************************
// ***************************************************

// [[Rcpp::export]]
NumericVector calcAlleleFreqPop(NumericVector alleles_1,
                                NumericVector alleles_2,
                                int Nallele){


  NumericVector freq_table(Nallele);
  float AlleleTotal = alleles_1.size() + alleles_2.size();

  // Loop through number of alleles and add them to frequency table
  for(int allele = 0; allele < alleles_1.size(); ++allele){
    freq_table[alleles_1[allele]] += 1;
    freq_table[alleles_2[allele]] += 1;
  }

  // Divide by total number of alleles to get frequency
  for(int allele = 0; allele < Nallele; ++allele){
    freq_table[allele] = freq_table[allele] / AlleleTotal;
  }

  return(freq_table);

}


// ***************************************************
// ***************************************************
// *** CALCULATE ALLELE FREQUENCY FOR AN INDIVIDUAL
// ***************************************************
// ***************************************************

// [[Rcpp::export]]
std::vector<float> calcAlleleFreqCppInd(int alleles_1,
                                        int alleles_2,
                                        int Nallele){


  std::vector<float> freq_table(Nallele);
  float AlleleTotal = 2; // FOr looking at individual alleles


  // Loop through number of alleles and add them to frequency table
  for(int allele = 0; allele < 1; ++allele){

    freq_table[alleles_1] += 1; // Minus 1 is because of zero indexing in Cpp
    freq_table[alleles_2] += 1;
  }

  // Divide by total number of alleles to get frequency
  for(int allele = 0; allele < Nallele; ++allele){
    freq_table[allele] = freq_table[allele] / AlleleTotal;
  }

  return(freq_table);

}






// ***************************
// ***************************
// ***************************  SPATIAL FUNCTIONS
// ***************************
// ***************************
// ***************************




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

// Looping through the pairwise matrix for each distance interval is likely not the most efficient way to do this

// [[Rcpp::export]]
NumericMatrix findDIs(NumericMatrix Mdij,
                      NumericVector distance_intervals,
                      int Nind) {

  NumericMatrix Mcij(Nind, Nind); // Matrix saying which distance class each pairwise combo is in

  for(int di = 0; di < distance_intervals.size(); di++){
    for(int i = 0; i < (Nind - 1); i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

      if(di == 0){ // Special loop for first interval because there's no lower interval
        if(Mdij(i, j) < distance_intervals[di]) Mcij(i, j) = di;
      }
      if(di > 0){ // Find interval where the pairwise distance fits inside
        if((Mdij(i, j) > distance_intervals[di - 1]) & (Mdij(i, j) < distance_intervals[di])) Mcij(i, j) = di;
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

// Looping through the pairwise matrix for each distance interval is likely not the most efficient way to do this

// [[Rcpp::export]]
NumericVector findEqualDIs(NumericMatrix Mdij,
                           NumericVector distance_intervals,
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
  int ndis = -distance_intervals[0]; // Take negative of negative

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
                           NumericVector distance_intervals,
                           int Nind){

  int ndis = distance_intervals.size();
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
    summary(1, i) = distance_intervals[i]; // Max distance in distance interval
    summary(2, i) = tot[i] / npairs[i]; // Average distance
    summary(3, i) = npairs[i]; // Number of pairs
  }

  return(summary);

} // End function












