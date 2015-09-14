


plotSgs <- function(sgsOut, overlay = FALSE, color = "black",
                    max_distance = NULL){

  # Check to make sure input is the right class
  if(!is(sgsOut, "sgsOut")){
    stop("Data must be of class sgsOut.. \n")
  }

    symbols <- 19

    estimate <- sgsOut$Fijsummary["ALL LOCI",] ## Save relatedness estimate
    dist <- sgsOut$DIsummary["Max distance", ] ## Save max distance intervals

    conf_hi = 0
    conf_low = 0

    if(is.null(max_distance)) max_distance = max(dist)

    ## Plotting permutations
    if(!is.null(sgsOut$PermAvg)){

      #estimate <- sgsOut$PermAvg["ALL_LOCI_perm_avg",] # chop off last 3 columns and 1st col
      conf_hi <- sgsOut$Perm975["ALL_LOCI_perm_975",]
      conf_low <- sgsOut$Perm025["ALL_LOCI_perm_025",]

    }

    if(overlay){
      par(new = TRUE)
      points(dist, estimate, type = "b", pch = symbols, xlab = "", ylab = "",
             yaxt = "n", xaxt = "n", col = color, lwd = 2, lty = 1)
#       if(!is.null(perm)){
#         lines(dist, conf_hi, lty = 4, col = color)
#         lines(dist, conf_low, lty = 4, col = color)
#       }

    } else{

      plot(dist, estimate, type = "b", pch = symbols, las = 1,
           ylab = "Kinship", xlab = "Distance (m)", lwd = 2,
           ylim = c(min(estimate, conf_low, na.rm = TRUE) * 1.1,
                    max(estimate, conf_hi, na.rm = TRUE) * 1.1),
           xlim = c(0, max_distance), col = color)
      abline(h = 0, lty = 1, col = "grey50")
      if(!is.null(conf_hi)){
         lines(dist, conf_hi, lty = 4, col = color)
         lines(dist, conf_low, lty = 4, col = color)
       }
    }

  } ## End plotting function




