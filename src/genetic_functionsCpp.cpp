#include <Rcpp.h>
using namespace Rcpp;


// Declare functions
float calcFijPairwiseCpp(NumericMatrix ref_gen,
                         NumericMatrix alfreq1,
                         NumericMatrix alfreq2,
                         int Nloci, int Nallele,
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
float calcFijPairwiseCpp(NumericMatrix ref_gen,
                         NumericMatrix alfreq1,
                         NumericMatrix alfreq2,
                         int Nloci, int Nallele,
                         int n){

  float denom = 0;
  float numer = 0;
  float fij;

  for(int locus = 0;  locus < Nloci; ++locus){ // Loop through loci

    for(int allele = 0; allele < Nallele; ++allele){ // Loop through alleles

      numer  +=  (alfreq1(locus, allele) - ref_gen(locus, allele))
      *  (alfreq2(locus, allele) -
        ref_gen(locus, allele)) +
        (ref_gen(locus, allele)*(1 - ref_gen(locus, allele))) / (n - 1);


      denom += ref_gen(locus, allele) * (1 - ref_gen(locus, allele));

    }
  }

  fij = numer / denom;

  return(fij);
}



// ***************************************************
// ***************************************************
// ****CALCULATE PAIRWISE FIJ AMONG ALL INDIVIDUALS IN A POPULATION
// ***************************************************
// ***************************************************

// [[Rcpp::export]]
NumericMatrix calcFijPopCpp(NumericVector ids,
                            NumericMatrix genotype_data,
                            NumericMatrix ref_gen,
                            int Nloci,
                            int Nallele,
                            int Nind,
                            int Ngenecopies){

  // Initalize a 3d array that will save allele frequency for each individual
  float ind_al_freq[Nloci][Nallele][Nind];
  int row = 0; // Needed for loop through 3d array
  NumericMatrix fij(Nind, Nind);
  NumericMatrix  alfreq1(Nloci, Nallele);
  NumericMatrix  alfreq2(Nloci, Nallele);
  int id1;
  int id2;


  // Loop through individuals and calculate allele frequency
  for(int indi = 0; indi < Nind; ++indi){

    row = 0; // Row of individual allele frequency 3d array - corresponds to locus

    for(int locus = 0; locus < (Nloci * 2); locus += 2){ // Loop through loci

      std::vector<float> freq_table = calcAlleleFreqCppInd(genotype_data(indi, locus),
                                                           genotype_data(indi, locus + 1),
                                                           Nallele);
      // Code to print out individual frequency tables
      /*
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
      for(int allele = 0; allele < Nallele; ++allele){
        alfreq1(locus, allele) = ind_al_freq[locus][allele][i];
        alfreq2(locus, allele) = ind_al_freq[locus][allele][j];
      }
    }



    fij(i, j) =   calcFijPairwiseCpp(ref_gen,
                                     alfreq1,
                                     alfreq2,
                                     Nloci, Nallele, Ngenecopies);



  } // End Pairwise Fij

  return(fij);

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
    freq_table[alleles_1[allele] - 1] += 1; // Minus 1 is because of zero indexing in Cpp
    freq_table[alleles_2[allele] - 1] += 1;
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
std::vector<float> calcAlleleFreqCppInd(int alleles_1,
                                        int alleles_2,
                                        int Nallele){


  std::vector<float> freq_table(Nallele);
  float AlleleTotal = 2; // FOr looking at individual alleles


  // Loop through number of alleles and add them to frequency table
  for(int allele = 0; allele < 1; ++allele){

    freq_table[alleles_1 - 1] += 1; // Minus 1 is because of zero indexing in Cpp
    freq_table[alleles_2 - 1] += 1;
  }

  // Divide by total number of alleles to get frequency
  for(int allele = 0; allele < Nallele; ++allele){
    freq_table[allele] = freq_table[allele] / AlleleTotal;
  }

  return(freq_table);

}



