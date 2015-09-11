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



// ***************************************************
// ***************************************************
// *** CALCULATE PAIRWISE FIJ BETWEEN A PAIR OF INDIVIDUALS
// ***************************************************
// ***************************************************

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

// [[Rcpp::export]]
NumericMatrix calcFijPopCpp(NumericMatrix Mcij,
                            NumericVector distance_intervals,
                            NumericMatrix genotype_data,
                            List ref_gen,
                            int Nloci,
                            NumericVector Nallele,
                            int Nind,
                            int Ngenecopies){

  // Initalize a 3d array that will save allele frequency for each individual
  int MNallele = max(Nallele);
  float ind_al_freq[Nloci][MNallele][Nind];
  int ndis = distance_intervals.size(); // Number of distance intervals
  int row = 0; // Needed for loop through 3d array
  std::vector<float> fij(Nloci);
  NumericMatrix  alfreq1(Nloci, MNallele);
  NumericMatrix  alfreq2(Nloci, MNallele);
  NumericMatrix fijsummary(Nloci + 1, ndis); // Stores summary information of pairwise Fij
  NumericMatrix npairs(Nloci + 1, ndis); // Stores number of comparisons

  // Initialize matrices
  for(int locus = 0; locus < (Nloci + 1); locus++) for(int di = 0; di < ndis; di++){
    fijsummary(locus, di) = 0;
    npairs(locus, di) = 0;
  }


  // Loop through individuals and calculate allele frequency
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


  // Calculate pairwise Fij between all pairs of individuals
  for(int i = 0; i < (Nind - 1) ; i++) for(int j = i+1; j <= (Nind - 1); j++){ // Loops through all pairs of individuals

    // Loop through and save allele frequency for each individual quickly
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

      fijsummary(locus, Mcij(i, j)) += *k;
      fijsummary(Nloci, Mcij(i, j)) += *k; // Save sum across all loci in last row

      npairs(locus, Mcij(i, j)) += 1;
      npairs(Nloci, Mcij(i, j)) += 1; // Save average in last row

      locus += 1;

    } // End locus loop
  } // End loop through pairs of individuals


  //Calculate average Fij per locus by dividing sum by total pairs
   for(int locus = 0; locus <= Nloci; locus++) for(int di = 0; di < ndis; di++){
      fijsummary(locus, di) = fijsummary(locus, di) / npairs(locus, di);
   }

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
// Needs integer inputs
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



