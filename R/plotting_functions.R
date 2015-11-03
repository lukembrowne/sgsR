
###############################
#### Plot results of SGS analysis
###############################

# This function plots the amount of spatial autocorrelation in relatedness among the samples,
# which is a typical graph made from these type of analyses

plot.sgsOut <- function(sgsOut, overlay = FALSE, color = "black",
                    max_distance = NULL, ...){

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
      points(dist, estimate, type = "b", pch = 19, xlab = "", ylab = "",
             yaxt = "n", xaxt = "n", col = color, lwd = 2, lty = 1)
      if(!is.null(conf_hi)){
        lines(dist, conf_hi, lty = 4, col = color)
        lines(dist, conf_low, lty = 4, col = color)
      }

    } else{

      ## Main plotting section
      plot(dist, estimate, type = "b", pch = 19, las = 1,
           ylab = "Kinship", xlab = "Distance (m)", lwd = 2,
           ylim = c(min(estimate, conf_low, na.rm = TRUE) * 1.1,
                    max(estimate, conf_hi, na.rm = TRUE) * 1.1),
           xlim = c(0, max_distance), col = color)
      abline(h = 0, lty = 1, col = "grey50") # Horizontal line where relatedness = 0

      if(!is.null(conf_hi)){ ## Add permutation lines
         lines(dist, conf_hi, lty = 4, col = color)
         lines(dist, conf_low, lty = 4, col = color)
       }
    }

  } ## End plotting function


