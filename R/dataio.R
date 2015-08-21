




### Building basic data structure


# Issues - how to deal with labeling of loci - genotype_data should have labeled columns?
  # Should it be a matrix or data frame?

## Set it up to be able to deal with varying number of alleles per locus

createSgsObj <- function(sample_ids,
                        groups = NULL,
                        genotype_data,
                        ploidy,
                        x_coords,
                        y_coords){

  ### Run checks to make sure that inputs are the proper format

  df <- structure(list(), class = "sgsObj")

  # Set sample ids and groups
  df$ids = sample_ids
  df$groups <- groups

  # Set genetic data
  df$gen_data <- genotype_data

  # Set spatial coordinates
  df$x <- x_coords
  df$y <- y_coords

  # Set various attributes
  df$ploidy = ploidy #  ploidy
  df$Nind = length(sample_ids) # Set number of individuals
  df$Nloci = ncol(genotype_data) / ploidy # Set number of loci
  df$Ngenecopies = Nind * ploidy # Will need to change when missing data is a thing

  # Find max number of alleles across all loci
  df$Nallele = max(sapply(genotype_data, FUN = function(x) length(table(x)))) # Might not work well with missing data

  return(df)

}

