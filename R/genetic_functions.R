


#### Main SGS analysis function

sgs <- function(sgsObj,
                distance_intervals,
                nperm = 999){

  # Check to make sure data is sgs object
  if(!is(sgsObj, "sgsObj")){
    stop("Input data must be of class sgsObj... \n")
  }

  ## Set up output data structure

  sgsOut <- structure(list(), class = "sgsOut")


  #####
  ## DISTANCE INTERVALS
  ####

    ## Calculate distance intervals
      Mdij = calcPairwiseDist(sgsObj$x, sgsObj$y, sgsObj$Nind ) ## Distance matrix

    ## Add a new max distance interval if last interval isn't already large enough
      if(max(Mdij) > max(distance_intervals)){
        distance_intervals <- c(distance_intervals, max(Mdij) * 1.001)
        cat("Adding an aditional distance interval to encompass all pairwise distances.. \n")
      }

    # Assign each pairwise combination to a distance interval
      Mcij = findDIs(Mdij, distance_intervals, sgsObj$Nind) ## Class

    # Calculate summary of distance intervals
      DIsummary = summarizeDIs(Mdij, Mcij, distance_intervals, Nind = sgsObj$Nind)
      rownames(DIsummary) <- c("Distance class", "Max distance", "Average distance",
                               "Number of pairs")

      sgsOut$DIsummary <- DIsummary


  #####
  ## REFERENCE ALLELE FREQUENCY
  ####

    ## Calculate reference allele frequencies
      ## Could probably speed this up by converting loop to Cpp
    ref_gen <- matrix(data = 0, nrow = Nloci, ncol = Nallele)

    row = 1
    for(locus in seq(1, (Nloci*2), by = 2)){
      ref_gen[row, ] = calcAlleleFreqPop(sgsObj$gen_data[ , locus],
                                         sgsObj$gen_data[ , locus + 1],
                                         Nallele = sgsObj$Nallele )
      row = row + 1
    }

  #####
  ## PAIRWISE RELATEDNESS COEFFICIENT
  ####

  ## Calculate pairwise Fij across entire population
  fijsummary = calcFijPopCpp(ids = sgsObj$ids,
                       genotype_data = as.matrix(sgsObj$gen_data),
                       distance_intervals = distance_intervals,
                       Mcij = Mcij,
                       ref_gen = ref_gen,
                       Nloci = sgsObj$Nloci,
                       Nallele = sgsObj$Nallele,
                       Nind = sgsObj$Nind,
                       Ngenecopies= sgsObj$Ngenecopies)

    rownames(fijsummary) <- c(sgsObj$loci_names, "ALL LOCI")

    sgsOut$Fijsummary <- fijsummary


  ## Output a summary table with distance classes, Fij estimates, permutation results, etc


  return(sgsOut)
}
