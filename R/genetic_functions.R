


#### Main SGS analysis function

sgs <- function(sgsObj,
                distance_intervals,
                nperm = 999){

  # Check to make sure data is sgs object
  if(!is(sgsObj, "sgsObj")){
    stop("Input data must be of class sgsObj... \n")
  }


  ## Calculate reference allele frequencies
    ## Could probably speed this up by converting loop to Cpp
  ref_gen <- matrix(data = 0, nrow = Nloci, ncol = Nallele)

  row = 1
  for(locus in seq(1, (Nloci*2), by = 2)){
    ref_gen[row, ] = calcAlleleFreqPop(pony$gen_data[ , locus],
                                       pony$gen_data[ , locus + 1],
                                       Nallele = pony$Nallele )
    row = row + 1
  }



  ## Calculate pairwise Fij across entire population
  fijs = calcFijPopCpp(ids = pony$ids,
                       genotype_data = as.matrix(pony$gen_data),
                       ref_gen = ref_gen,
                       Nloci = pony$Nloci,
                       Nallele = pony$Nallele,
                       Nind = pony$Nind,
                       Ngenecopies= pony$Ngenecopies)


  ## Output a summary table with distance classes, Fij estimates, permutation results, etc


  return(fijs)
}
