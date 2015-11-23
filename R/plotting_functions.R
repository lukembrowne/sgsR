
###############################
#### Plot results of SGS analysis
###############################

#' Plot spatial autocorrelation graph
#'
#'Plots a spatial autocorrelation graph with Kinship estimator on y axis and distance on x axis. Each point represents the average kinship estimate within that distance interval. Dashed lines represent 95% confidence intervals for the permuted kinship estimates (if permutations were performed in the \code{\link{sgs}} analysis).
#'
#' @param x An sgsOut objects that contains results of an \code{\link{sgs}} analysis
#' @param overlay If TRUE, overlay plot on current plot.
#' @param max_distance Maximum distance to display on x axis. If NULL, defaults to value of largest distance interval
#' @param lwd_avg Line width of line connecting kinship estimates in each distance interval
#' @param col_avg Color of points and line connecting kinship estimates in each distance interval
#' @param pch Point type for average value of kinsihp estimates in each distance interval
#' @param lwd_CI Line width of confidence interval lines
#' @param lty_CI Line type of confidence interval lines
#' @param col_CI Color of confidence interval lines
#' @param ylab Y axis label
#' @param xlab X axis label
#'
#' @param ... Arguments to pass to generic plotting functions
#'
#' @return
#' A spatial autocorrelation plot
#' @export
#'
#' @examples
#' 1
plot.sgsOut <- function(x, overlay = FALSE,
                        max_distance = NULL,
                        lwd_avg = 2,
                        col_avg = "black",
                        pch = 19,
                        lwd_CI = 1,
                        lty_CI = 4,
                        col_CI = "black",
                        ylab = "Kinship",
                        xlab = "Distance",
                        ...){

  sgsOut = x
  # Check to make sure input is the right class
  if(!is(sgsOut, "sgsOut")){
    stop("Data must be of class sgsOut.. \n")
  }

    ## Set values
    estimate <- sgsOut$fij_obs["ALL LOCI",] ## Save relatedness estimate
    dist <- sgsOut$di["Max distance", ] ## Save max distance intervals

    conf_hi = 0 # Placeholders
    conf_low = 0

    if(is.null(max_distance)) max_distance = max(dist) # If max distance not set, set it

    ## Saving permutation results
    if(!is.null(sgsOut$fij_perm_avg)){ # If permutation results found, continue..

      conf_hi <- sgsOut$fij_perm_975["ALL_LOCI_perm_975",]
      conf_low <- sgsOut$fij_perm_025["ALL_LOCI_perm_025",]

    }

    if(overlay){ # Option to overlay plot on top of another
      par(new = TRUE)
      points(dist, estimate, type = "b", pch = pch, xlab = "", ylab = "",
             yaxt = "n", xaxt = "n", col = col_avg, lwd = lwd_avg, lty = 1, ...)
      if(!is.null(sgsOut$fij_perm_avg)){
        lines(dist, conf_hi,  lty = lty_CI, col = col_CI, lwd = lwd_CI)
        lines(dist, conf_low, lty = lty_CI, col = col_CI, lwd = lwd_CI)
      }

    } else{

      ## Main plotting section
      plot(dist, estimate, type = "b", las = 1, pch = pch,
           ylab = ylab, xlab = xlab, lwd = lwd_avg,
           ylim = c(min(estimate, conf_low, na.rm = TRUE) * 1.1,
                    max(estimate, conf_hi, na.rm = TRUE) * 1.1),
           xlim = c(0, max_distance), col = col_avg, ...)
      abline(h = 0, lty = 1, col = "grey50") # Horizontal line where relatedness = 0

      if(!is.null(sgsOut$fij_perm_avg)){ ## Add permutation lines
         lines(dist, conf_hi,  lty = lty_CI, col = col_CI, lwd = lwd_CI)
         lines(dist, conf_low, lty = lty_CI, col = col_CI, lwd = lwd_CI)
       }
    }

  } ## End plotting function



###############################
#### Plot sgsObj data structure
###############################

#' Plot spatial locations of individuals in sgsObj data structure
#'
#' @param x An sgsObj data structure to be plotted
#' @param pch Type of plotting character to display
#' @param ... Any arguments to pass to the generic plot function
#'
#' @return
#' Plots the spatial location of individuals in the sgsObj data structure.
#' @export
#'
#' @examples
#'
#' ## Simulate genetic data
#' Nind = 100 # Number of individuals
#' Nloci = 5 # Number of loci
#' Nallele = 10 # Number of alleles per loci
#'
#' ## Set up data frame and generate random spatial locations
#' dat <- data.frame(id = 0:(Nind - 1))
#' dat$x = runif(Nind, 0, 100)
#' dat$y = runif(Nind, 0, 100)
#'
#' ## Simulate Random genetic data and assign loci names
#' for (loci in 1:Nloci) {
#'  loci_name_a = paste("Loc", loci, "_A", sep = "")
#'  loci_name_b = paste("Loc", loci, "_B", sep = "")
#'  dat[loci_name_a] <- sample.int(Nallele, Nind, replace = TRUE)
#'  dat[loci_name_b] <- sample.int(Nallele, Nind, replace = TRUE)
#'}
#'
#' ## Convert to sgsObj
#' sgsObj = createSgsObj(sample_ids = dat$id,
#'                      genotype_data = dat[, 4:(Nloci*2 + 3)],
#'                      ploidy = 2,
#'                      x_coords = dat$x,
#'                      y_coords = dat$y,
#'                      missing_val = -999)
#'
#'summary(sgsObj)
#'
#'plot(sgsObj)
#'
#'plot(sgsObj, pch = 21, col = "steelblue", cex = 1.5)
plot.sgsObj <- function(x, pch = 19, ...){

  sgsObj = x

  # Check to make sure input is the right class
  if(!is(sgsObj, "sgsObj")){
    stop("Data must be of class sgsObj.. \n")
  }

    ## Main plotting section
    plot(sgsObj$x, sgsObj$y, pch = pch, las = 1,
         ylab = "Y coordinates", xlab = "X coordinates",
         asp = 1, ...)
} ## End plotting function





