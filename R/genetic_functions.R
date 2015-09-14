


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

      # If equalized distance interval option is chosen (by setting option to negative)
      # Find new distance intervals with approximately equal number of pairwise comparisons
    if(distance_intervals[1] < 0){
      cat("Finding --", -distance_intervals, "-- distance intervals with approximately equal pairwise comparisons...\n")
      distance_intervals =  findEqualDIs(Mdij, distance_intervals, sgsObj$Nind)
    }

      ## Add a new max distance interval if last interval isn't already large enough
      if(max(Mdij) > max(distance_intervals)){
        distance_intervals <- c(distance_intervals, max(Mdij) * 1.001)
        cat("Adding an aditional distance interval --", max(Mdij) * 1.001,
            "-- to encompass all pairwise distances.. \n")
      }

      Mcij = findDIs(Mdij, distance_intervals, sgsObj$Nind)

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
    ref_gen <- matrix(data = NA, nrow = sgsObj$Nloci, ncol = max(sgsObj$Nallele))

    row = 1
    for(locus in seq(1, (sgsObj$Nloci*2), by = 2)){
      ref_gen[row, 1:sgsObj$Nallele[row]] = calcAlleleFreqPop(sgsObj$gen_data_f[, locus],
                                         sgsObj$gen_data_f[ , locus + 1],
                                         Nallele = sgsObj$Nallele[row] )
      row = row + 1
    }

    ## reference genotypes as a list
    ref_gen <- list()
    row = 1
    for(locus in seq(1, (sgsObj$Nloci*2), by = 2)){
      ref_gen[[row]] = calcAlleleFreqPop(sgsObj$gen_data_f[, locus],
                                         sgsObj$gen_data_f[ , locus + 1],
                                         Nallele = sgsObj$Nallele[row])
      row = row + 1
    }


  #####
  ## PAIRWISE RELATEDNESS COEFFICIENT
  ####

  ## Calculate pairwise Fij across entire population

    ### Save x and y coordinates in their own vector so that random_shuffle doesn't overwrite
    xcopy = rep(NA, length(sgsObj$x))
    xcopy = replace(xcopy, 1:length(sgsObj$x), sgsObj$x )

    ycopy = rep(NA, length(sgsObj$y))
    ycopy = replace(ycopy, 1:length(sgsObj$y), sgsObj$y )

  fijsummary = calcFijPopCpp(genotype_data = as.matrix(sgsObj$gen_data_f),
                       distance_intervals = distance_intervals,
                       Mdij = Mdij,
                       Mcij = Mcij,
                       ref_gen = ref_gen,
                       Nloci = sgsObj$Nloci,
                       Nallele = sgsObj$Nallele,
                       Nind = sgsObj$Nind,
                       Ngenecopies= sgsObj$Ngenecopies,
                       Nperm = nperm,
                       x_coord = xcopy,
                       y_coord = ycopy)

    rownames(fijsummary) <- c(sgsObj$loci_names, "ALL LOCI",
                              paste(sgsObj$loci_names, "_perm", sep = ""), "ALL_LOCI_perm")

    sgsOut$Fijsummary <- fijsummary


  ## Output a summary table with distance classes, Fij estimates, permutation results, etc


  return(sgsOut)
}
