

###############################
#### Main SGS analysis function
###############################

## This is the meat of the package, the main function to do the main analysis.
# It performs several steps along the way (finding distance intervals, calculating SGS,
# conducting permutation tests). Calls several C++ functions to speed up the process.

### Make it clear what checkDi's is looking for and how to interpret the cautionary values.

# Returns an sgsOut object that holds summary information about the results of the analysis

#' Title
#'
#' @param sgsObj
#' @param distance_intervals
#' @param nperm
#'
#' @return 1
#' @export
#'
#' @examples 1
sgs <- function(sgsObj,
                distance_intervals,
                nperm = 999){

  # Check to make sure data is sgs object
  if(!is(sgsObj, "sgsObj")){
    stop("Input data must be of class sgsObj... \n")
  }

  ## Set up output data structure
  sgsOut <- structure(list(), class = "sgsOut")

  sgsOut$sgsObj <- sgsObj

  #####
  ## DISTANCE INTERVALS
  ####

    ## Calculate distance intervals
      Mdij = calcPairwiseDist(sgsObj$x, sgsObj$y, sgsObj$Nind ) ## Distance matrix - C++ func

      # If equalized distance interval option is chosen (by setting option to negative)
      # Find new distance intervals with approximately equal number of pairwise comparisons
    if(distance_intervals[1] < 0){
      cat("Finding --", -distance_intervals, "-- distance intervals with approximately equal pairwise comparisons...\n")
      distance_intervals =  findEqualDIs(Mdij, distance_intervals, sgsObj$Nind) # C++ func
    }

      ## Add a new max distance interval if last interval isn't already large enough
      if(max(Mdij) > max(distance_intervals)){
        distance_intervals <- c(distance_intervals, max(Mdij) * 1.001)
        cat("----------------------------\n")
        cat("Adding an aditional distance interval --", max(Mdij) * 1.001,
            "-- to encompass all pairwise distances.. \n\n")
      }

      ## Make sure that first distance interval include all pairs of neighbors. As suggested by
      ## Vekemans and Hardy 2004 for calculating Sp statistics

        # Make distance matrix symmetric
        dist_mat <- Mdij + t(Mdij)
        diag(dist_mat) <- NA # Set diagonal to NA

        nn_dist <- apply(dist_mat, 1, min, na.rm = TRUE)

        if(any(nn_dist > distance_intervals[1])){
          cat("----------------------------\n")
          cat("Caution - first distance interval does not include all pairs of neighbors, as suggested by Vekemans and Hardy 2004 for calculating Sp statistic. \n")
          cat("First distance interval would have to be at least:", max(nn_dist), "to include all pairs of neighbors. \n\n")
        }


       ### Find distance intervals
      Mcij = findDIs(Mdij, distance_intervals, sgsObj$Nind) # C++ func

    # Calculate summary of distance intervals
      DIsummary = summarizeDIs(Mdij, Mcij, distance_intervals, Nind = sgsObj$Nind) # C++ func
      rownames(DIsummary) <- c("Distance class", "Max distance", "Average distance",
                               "Number of pairs", "% participation", "CV participation")

      sgsOut$di <- round(DIsummary, 2) ## Save summary information about distance intervals

      ## Add column names to di output
      colnames(sgsOut$di) <- paste("D", 1:length(distance_intervals), sep = "")

      ## Check to make sure DIs follow rules of thumbs in Spagedi manual
      checkDIs(sgsOut$di)



  #####
  ## REFERENCE ALLELE FREQUENCY
  ####


    ## reference genotypes as a list
    ##
      ## Could probably speed up a bit by converting loop to C++
    ref_gen <- list()
    row = 1
    for(locus in seq(1, (sgsObj$Nloci*2), by = 2)){
      ref_gen[[row]] = calcAlleleFreqPop(sgsObj$gen_data_int[, locus],
                                         sgsObj$gen_data_int[ , locus + 1],
                                         Nallele = sgsObj$Nallele[row],
                                         Ngenecopies = sgsObj$Ngenecopies[row])
      row = row + 1
    }

    ## Check to make sure allele frequencies at each locus sum to one
    if(any(sapply(ref_gen, sum) != 1)){
      warning("Reference allele frequencies do not sum to 1... \n")
    }


  #####
  ## PAIRWISE RELATEDNESS COEFFICIENT
  ####

  ## Calculate pairwise Fij across entire population
  # This is a huge C++ function that perfroms the brunt of the SGSG analysis

  ## Returns a matrix with information on Fij estimates, permutation results, for each loci,
    # averaged across loci for each distance class
    ## number of columns is equal to length of distance_intervals

    ## Sets of data - each set is Nloci + 1 rows long
      # 1) Fij - observed by distance class
      # 2) Fij - Permutation average by distance class
      # 3) Fij - 2.5% permutation quantile by distance class
      # 4) Fij - 97.5% permutation quantile by distance class
      # 5) Sp - observed values - only first column
      # 6) Slope - observed values - only first clumn

  fijsummary = calcFijPopCpp(genotype_data = as.matrix(sgsObj$gen_data_int),
                       distance_intervals = distance_intervals,
                       Mdij = Mdij,
                       Mcij = Mcij,
                       ref_gen = ref_gen,
                       Nloci = sgsObj$Nloci,
                       Nallele = sgsObj$Nallele,
                       Nind = sgsObj$Nind,
                       Ngenecopies = sgsObj$Ngenecopies,
                       Nperm = nperm,
                       x_coord = sgsObj$x,
                       y_coord = sgsObj$y)



    rownames(fijsummary) <- c(sgsObj$loci_names, "ALL LOCI",
                        paste(sgsObj$loci_names, "_perm_avg", sep = ""), "ALL_LOCI_perm_avg",
                        paste(sgsObj$loci_names, "_perm_025", sep = ""), "ALL_LOCI_perm_025",
                        paste(sgsObj$loci_names, "_perm_975", sep = ""), "ALL_LOCI_perm_975",
                        paste(sgsObj$loci_names, "_sp", sep = ""), "ALL_LOCI_sp",
                        paste(sgsObj$loci_names, "_slope", sep = ""), "ALL_LOCI_slope")

    ## Add column names to fijsummary
    colnames(fijsummary) <- paste("D", 1:length(distance_intervals), sep = "")

    ## Subset out specific sections
    sgsOut$fij_obs <- fijsummary[1:(sgsObj$Nloci + 1), ]

    sgsOut$sp_obs <- fijsummary[(sgsObj$Nloci*4 + 5):(sgsObj$Nloci*5 + 5), 1, drop = FALSE ]
    colnames(sgsOut$sp_obs) <- c("Observed")
    sgsOut$slope_obs <- fijsummary[(sgsObj$Nloci*5 + 6):(sgsObj$Nloci*6 + 6), 1, drop = FALSE ]
    colnames(sgsOut$slope_obs) <- c("Observed")

    if(nperm > 0){ # If permutation selected, add it to the output
      sgsOut$fij_perm_avg <- fijsummary[(sgsObj$Nloci + 2):(sgsObj$Nloci*2 + 2), ]
      sgsOut$fij_perm_025 <- fijsummary[(sgsObj$Nloci*2 + 3):(sgsObj$Nloci*3 + 3), ]
      sgsOut$fij_perm_975<- fijsummary[(sgsObj$Nloci*3 + 4):(sgsObj$Nloci*4 + 4), ]
    }


  ## Output a summary table with distance classes, Fij estimates, permutation results, etc
   return(sgsOut)

}
