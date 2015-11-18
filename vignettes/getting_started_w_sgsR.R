## ---- fig.height = 3, fig.width = 5--------------------------------------

library(sgsR)

## Simulate genetic data
    Nind = 100
    Nloci = 10
    Nallele = 10
    n =  Nind * 2 # Number of gene copies

    ## Initialize data frame
    dat <- data.frame(id = 0:(Nind-1))
    dat$x = runif(Nind, 0, 100)
    dat$y = runif(Nind, 0, 100)

    ## Simulate Random genetic data
    for(loci in 1:Nloci){
      loci_name_a = paste("Loc", loci, "_A", sep = "")
      loci_name_b = paste("Loc", loci, "_B", sep = "")
      dat[loci_name_a] <- sample.int(Nallele, Nind, replace = TRUE)
      dat[loci_name_b] <- sample.int(Nallele, Nind, replace = TRUE)
    }

## Convert to sgsObj
sgsObj = createSgsObj(sample_ids = dat$id, 
                      genotype_data = dat[, 4:(Nloci*2 + 3)],
                      ploidy = 2,
                      x_coords = dat$x, 
                      y_coords = dat$y)

# Display genetic data
head(sgsObj$gen_data)


## Run analysis
distance_intervals = seq(10, 110, 10) # Set distance intervals

out1 = sgs(sgsObj = sgsObj, distance_intervals = distance_intervals, nperm = 99)




## Plotting results

## Solid line is Fij estimate for each distance class
## Dashed lines are the 2.5 % and 97.5 % quantiles of the permuted values
plot(out1)

# Summary of information on distance classes
out1$di

# Summary of information on estimated Kinship coefficient for each distance class (columns)
round(out1$fij_obs, 3)

 

