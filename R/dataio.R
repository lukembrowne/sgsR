#################
##### Building basic data structure - sgsObj
###########

### sgsObj is the data structure that will hold all the information necessary for the sgs analysis
# It includes the genotype data, the location data, information about ploidy, groups, names of loci
# Can be created from scratch through the function below, or by reading in text input files created by other programs

createSgsObj <- function(sample_ids,
                        groups = NULL,
                        genotype_data,
                        ploidy,
                        loci_names = NULL,
                        x_coords,
                        y_coords){

  ### Run checks to make sure that inputs are the proper format

  df <- structure(list(), class = "sgsObj")

  # Set sample ids and groups
  df$ids = sample_ids
  df$groups <- groups ## Groups functionality not currently implemented..

  # Set spatial coordinates
  df$x <- x_coords
  df$y <- y_coords

  # Set various attributes
  df$ploidy = ploidy #  ploidy
  df$Nind = length(sample_ids) # Set number of individuals
  df$Nloci = ncol(genotype_data) / ploidy # Set number of loci

  # Set genetic data and loci names
  df$gen_data <- genotype_data
  if(is.null(loci_names)) df$loci_names <- colnames(genotype_data)[seq(1, Nloci*2, 2)]
  if(!is.null(loci_names)) df$loci_names <- loci_names

    # Format genetic data so that uses integers instead of raw numbers,
    # used later for indexing in C++

    # Exclude missing data!
  df$gen_data_int = df$gen_data

  for(col in seq(1, df$Nloci * df$ploidy, df$ploidy)){

    ## Save indices of missing data, should be same for both columns of the locus
    miss_ind <- which(df$gen_data_int[ , col ] == -999)

      # Save levels of alleles
    lev = levels(as.factor(c(df$gen_data_int[, col], df$gen_data_int[, col + 1])))

     if(any(lev == "-999")){ # If any missing data present
      lev <- lev[-which(lev == "-999")] # Remove missing data from levels list
     }

      # Reformat as integer
    df$gen_data_int[, col] = match(df$gen_data_int[, col], lev) - 1 # Minus 1 for 0 indexing
    df$gen_data_int[, col + 1] = match(df$gen_data_int[, col + 1], lev) - 1

      # Add back in missing data
    df$gen_data_int[miss_ind, col] <- -999
    df$gen_data_int[miss_ind, col + 1] <- -999
  }


  # Find max number of alleles across all loci
  i = 1
  for(col in seq(1, df$Nloci * df$ploidy, df$ploidy)){

    inc_missing <- table(c(genotype_data[, col], genotype_data[, col + 1]))

    if(any(names(inc_missing) == "-999")){
      exc_missing = inc_missing[-which(names(inc_missing) == "-999")]
    } else {
      exc_missing = inc_missing
    }

    df$Nallele[i] = length(exc_missing)

    ## Save number of gene copies separately for each locus, excluding missing data
    df$Ngenecopies[i] <- sum(exc_missing)
    i = i + 1
  }

  names(df$Nallele) = df$loci_names

  return(df)

}





#################
##### Read spagedi input file and convert to sgsObj
###########

## This function reads in a spagedi text input file and converts it to an sgsObj

## Make it clear that data should follow format in example Spagedi data..
## Alleles with no non-numeric characters in between
# Only tested for diploid and 2d coordinate system

readSpagedi <- function(path_to_spagedi_file, missing_val = "-999"){

  lines <- readLines(path_to_spagedi_file) # Read lines in output file
  lines <- lines[-grep("^/", lines)] # Remove comment lines - start with / character
  lines <- lines[1:(min(grep("END", lines))-1)] # All info before first END statement


  ## 1st non comment line : set of 6 numbers separated by a tabulation and representing :		#individuals	#categories	#coordinates	#loci	#digits/allele	max ploidy

    first_line <- as.numeric(strsplit(lines[1], "\t")[[1]]) ## Returns vector with info from first non-comment line

    Nind =          first_line[1]
    cat("Number of individuals:", Nind, "\n" )
    Ncats =         first_line[2]
    cat("Number of categories:", Ncats, "\n" )
    Ncoords =       first_line[3]
    cat("Number of spatial dimensions:", Ncoords, "\n" )
    Nloci =         first_line[4]
    cat("Number of loci:", Nloci, "\n" )
    Ndigits =       first_line[5]
    cat("Number of digits per allele:", Ndigits, "\n" )
    Ploidy =        first_line[6]
    cat("Ploidy level:", Ploidy, "\n" )

   # 2nd non comment line : # of distance intervals followed by the upper distance of each interval. Here a negative # of intervals is given so that upper distances are choosen to obtain intervals with the same # of pairwise comparisons.

    second_line <- na.omit(as.numeric(strsplit(lines[2], "\t")[[1]]))

    Ndis = second_line[1] # Assign number of distance intervals
    cat("Number of distance intervals:", Ndis, "\n")
    # If specific distance intervals are assigned, save those
    if(length(second_line) > 1) {
      Dis = second_line[2:length(second_line)]
      cat("Distance intervals:", Dis, "\n")
    }

  ## 3rd non comment line : column labels (<=15 characters).
  ## 4th and next lines : data for each individual in the following order: Ind name, category, coordinates (here 0), genotypes

    third_line <- na.omit(strsplit(lines[3], "\t")[[1]])
    data = strsplit(lines[4:length(lines)], split = "\t")

    genotype_data = data.frame(rep(NA, Nind)) ## Initialize genotype data frame


  ## If no categories or spatial information
    if(Ncats == 0 & Ncoords == 0){
      cat("No category or spatial information detected...\n")
      loci_names = third_line[2:length(third_line)]
      ids = sapply(data, "[[", 1) ## First column is id information

      ## Fill in genotype data
      col2 = 1
      for(col in 2:(length(data[[1]]))){
        genotype_data[, col2] = sapply(data, "[[", col)
        col2 = col2 + 1
      }
    }

  ## If categories but not spatial information
    if(Ncats > 0 & Ncoords == 0){
      cat("Category but no spatial information detected...\n")
      loci_names = third_line[3:length(third_line)]
      ids = sapply(data, "[[", 1) ## First column is id information
      cats = sapply(data, "[[", 2) ## Category labels

      ## Fill in Genotype data
      col2 = 1
      for(col in 3:(length(data[[1]]))){
        genotype_data[, col2] = sapply(data, "[[", col)
        col2 = col2 + 1
      }
    }

  ## With categories and spatial information
    if(Ncats > 0 & Ncoords > 0){
      cat("Categories and spatial information detected...\n")
      loci_names = third_line[(3 + Ncoords):length(third_line)]
      ids = sapply(data, "[[", 1) ## First column is id information
      cats = sapply(data, "[[", 2) ## Category labels
      x = sapply(data, "[[", 3) ## X coordinates
      y = sapply(data, "[[", 4) ## X coordinates

      ## Fill in Genotype data
      col2 = 1
      for(col in 5:(length(data[[1]]))){
        genotype_data[, col2] = sapply(data, "[[", col)
        col2 = col2 + 1
      }
    }


    ## Split genotype data into separate columns for each allele
    split_gen_data <- data.frame(rep(NA, Nind))
    for(col in 1:ncol(genotype_data)){

      split_gen_data <- cbind(split_gen_data,
                              split_alleles(genotype_data[, col], Ndigits, Ploidy))
    }

    # Remove first column of split gen_data (All Nas..)
    split_gen_data <- split_gen_data[, -1]

    ## Convert to numeric
    for(col in 1:ncol(split_gen_data)){
      split_gen_data[, col] <- as.numeric(split_gen_data[, col])
    }

    ## Make sure number of columns matche up with Nloci and Ploidy
    if(ncol(split_gen_data) != (Nloci * Ploidy)){
      stop("Error in splitting genotype data.. number of columns do not match number of loci and ploidy level \n")
    }

    ## Replace with -999 for missing values
    missing_val = as.character(missing_val) ## In case it's inputted as a numeric
    split_gen_data[split_gen_data == missing_val] <- -999


    ### Assemble sgsObj
    out <- createSgsObj(sample_ids = ids,
                        groups = cats,
                        genotype_data = split_gen_data,
                        ploidy = Ploidy,
                        loci_names = loci_names,
                        x_coords = as.numeric(x),
                        y_coords = as.numeric(y))

    ## Final error checking and success message
    # Check loci names matches up with N loci
    if(length(out$loci_names) != Nloci){
      stop("Error: Number of loci does not match with length of loci names... \n")
    }

    if(length(out$ids) != Nind){
      stop("Error: Length of sample ids does not match number of individuals... \n")
    }

    cat("Successfully read SPAGeDi input file!\n")
   return(out)
}

########
##This function reads in a Genepop text input file and converts it to an sgsObj
##File must be in correct Genepop format
##Currently only works with ploidy of 2, as well as having individual id in first column

readGenepop <- function(path_to_genefile, missing_vals= "-999"){
  lines <- readLines(path_to_genefile) #Reads in file, stores in lines
  lines <- lines[-1] #Deletes title row
  Ncoords= 0 #No spatial data included in Genepop file format
  Ncats = 1 #Default 1 category
  Ploidy= 2 #Assuming all individuals have ploidy of 2
  Loci= 0

  #Calculates loci number
  Break= "Pop|pop|POP"
  for(x in 1:length(lines)){
    if (grepl(lines[x], Break) ){
      Loci= Loci +1
    }
  }

  lines <- lines[-grep('^P', lines)] #removes pop separators
  dataz = strsplit(lines[(Loci+1):length(lines)], split = "\\s+") #Need variable for 4, number of loci+1
  Nind= length(dataz) #The total number of individuals in data set

  loci_names <- '' #Initializes Loci names
  catz<- '' #Initializes Categories



  genotype_data = data.frame(rep(NA, Nind)) ## Initialize genotype data frame

  #Fill in genotype data
    if(Ncats > 0  & Ncoords == 0){
      cat("No spatial information detected...\n")
      loci_names = lines[1:Loci]
      ids = sapply(dataz, "[[", 1) ## First column is id information
      catz = sapply(dataz, "[[", 2) ## Category labels
      ## Fill in genotype data
      col2 = 1
      for(x in 1:(length(dataz[[1]]))){
        genotype_data[, col2] = sapply(dataz, "[[", x)
        col2 = col2 + 1
      }
    }

  #Removes Individual numbers and Category
  genotype_data<- genotype_data[,-2]
  genotype_data<- genotype_data[,-1]

  #Number of digits per allel
  Ndigits= (nchar(genotype_data[1,1]))/2

  ## Split genotype data into separate columns for each allele
  split_gen_dataz <- data.frame(rep(NA, Nind))
  for(col in 1:ncol(genotype_data)){
    split_gen_dataz <- cbind(split_gen_dataz,
                            split_alleles(genotype_data[, col], Ndigits, Ploidy))
  }

  # Remove first column of split gen_data (All Nas..)
  split_gen_dataz <- split_gen_dataz[, -1]



  ## Convert to numeric
  for(col in 1:ncol(split_gen_dataz)){
    split_gen_dataz[, col] <- as.numeric(split_gen_dataz[, col])
  }



  ## Replace with -999 for missing values
  missing_vals = as.character(missing_vals) ## In case it's inputted as a numeric
  split_gen_dataz[split_gen_dataz == missing_vals] <- -999


  ### Assemble sgsObj
  out <- createSgsObj(sample_ids = ids,
                      groups = catz,
                      genotype_data = split_gen_dataz,
                      ploidy = Ploidy,
                      loci_names = loci_names,
                      x_coords = NULL,
                      y_coords = NULL)


  return(out)

  }












########
### Function to split apart concatenated alleles - 0606 into 6 and 6 in separate columns
##########


split_alleles <- function(column_to_split, Ndigits, Ploidy){

  ## Have to maintain missing data
  ## Leading 0s are left off in Spagedi input data sometimes


  ## Add leading zeros for those that are missing
  # And find missing data

  char_lengths = sapply(column_to_split, nchar)

  missing_indices = which(char_lengths == 1)

  add_ld_zero = which(char_lengths > 1 & char_lengths < Ndigits*Ploidy)
  cat("Adding leading '0' to ", length(add_ld_zero), " genotypes...\n")
  column_to_split[add_ld_zero] <- paste("0", column_to_split[add_ld_zero], sep = "")

  # Initialize both columns with blank data
  out <- data.frame(col1 = rep("0", length(column_to_split)))

  # Split into individual characters
  split <- sapply(column_to_split, function(x) strsplit(x, split = NULL))

  ## Reconcatonate them
  col1 = 1 # For split column
  col2 = 1 # For output
  for(col1 in seq(1, Ploidy*Ndigits, Ndigits)){

    ## Separate into individual digits
    dig1 = sapply(split, '[', col1)

    dig2 = NULL
    dig3 = NULL

    if(Ndigits > 1) dig2 = sapply(split, '[', col1 + 1)
    if(Ndigits > 2) dig3 = sapply(split, '[', col1 + 2)

   out[, col2] = paste(dig1, dig2, dig3, sep="")
   col2 = col2 + 1
  }

  ## Replace missing data
  out[missing_indices, ] <- "0"


  return(out)
}





