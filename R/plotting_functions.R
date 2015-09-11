


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

    ## Include once permutation is included
#     if(!is.null(perm)){
#       perm <- perm[, c(2:(ncol(perm) - 4))] # chop off last 3 columns and 1st col
#
#       estimate <- as.numeric(perm["Obs val", ])
#       conf_hi <- as.numeric(perm["95%CI-sup", ])
#       conf_low <- as.numeric(perm["95%CI-inf", ])
#
#       # Closed symbol if permutation was significant
#       sig <- apply(perm[c(8,9,10), ], 2, FUN = function(x) {any(x < 0.05)})
#       symbols <- rep(1, ncol(perm))
#       symbols[sig] <- 19


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
#       if(!is.null(perm)){
#         lines(dist, conf_hi, lty = 4, col = color)
#         lines(dist, conf_low, lty = 4, col = color)
#       }
    }

  } ## End plotting function




