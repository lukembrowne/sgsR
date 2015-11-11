#################
##### Building basic data structure - sgsObj
###########

#' Create basic data structure for sgs analysis
#'
#'This function creates the basic data structure that will hold all the information necessary for sgs analyses - an object of the class sgsObj. It includes genotype data, spatial data, information about ploidy, groups, names of loci, etc.
#' sgsObjs can be created from scratch through the function below, or by reading in text input files created by other programs using the \code{\link{readSpagedi}} and \code{\link{readGenepop }} functions
#'
#' @param sample_ids A vector of numbers or names giving individual IDs of samples.
#' @param genotype_data A matrix or dataframe with genetic data. See details for more information.
#' @param x_coords A numeric vector with X coordinates of samples.
#' @param y_coords A numeric vector with Y coordinates of samples.
#' @param missing_val Numeric or character value indicating value used for missing data. Default is -999.
#' @param groups An optional vector of  of grouping classfication. \emph{Note:} Analyses on groups have not been implemented yet.
#' @param ploidy Ploidy level. \emph{Note:} Analyses have only been tested with ploidy = 2.
#' @param loci_names Optional vector of loci names.
#'
#'@section Details:
#'
#'\emph{\strong{Genotype data structure}}
#'
#'Genotype data should be formatted as a dataframe or matrix with individual samples as rows and loci as columns. Each locus should have two columns (in the case of diploid organisms), and these pairs of columns should be adjacent to each other.
#'
#'\tabular{rrrrrr}{
#'  Loc1_A \tab Loc1_B \tab Loc2_A \tab Loc2_B \tab Loc3_A \tab Loc3_B \cr
#' 1 \tab 0 \tab 7 \tab 0 \tab 1 \tab 0 \cr
#' 2 \tab 3 \tab 4 \tab 5 \tab 4 \tab 1 \cr
#' 3 \tab 2 \tab 2 \tab 3 \tab 2 \tab 4 \cr
#' 2 \tab 1 \tab 1 \tab 2 \tab 2 \tab 2 \cr
#' }
#' @return
#'
#' Returns an object of the class sgsObj.
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
#'                      y_coords = dat$y)
#'attributes(sgsObj)
createSgsObj <- function(sample_ids,
                         genotype_data,
                         x_coords,
                         y_coords,
                         missing_val = -999,
                        groups = NULL,
                        ploidy = 2,
                        loci_names = NULL
                        ){

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
  if(is.null(loci_names)) df$loci_names <- colnames(genotype_data)[seq(1, df$Nloci*2, 2)]
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

  ## Add column names to genotype data

  colnames(df$gen_data)<- colnames(genotype_data)
  colnames(df$gen_data_int) <- colnames(genotype_data)

  ## Look for loci that only have one allele and remove from dataframe

  rem_locus <- which(df$Nallele == 1)

  if(length(rem_locus) > 0 ){
        cat("Removing locus -- ", df$loci_names[rem_locus],
            "-- because it only has one allele. \n")

    # Find column index
        col_locus <- match(df$loci_names[rem_locus], colnames(df$gen_data))
        print(col_locus)

          # Remove from datasets
        df$gen_data <- df$gen_data[, - c(col_locus, col_locus + 1)]
        df$gen_data_int <- df$gen_data_int[, - c(col_locus, col_locus + 1)]
        df$Nallele <- df$Nallele[-rem_locus]
        df$loci_names <- df$loci_names[-rem_locus]
        df$Ngenecopies <- df$Ngenecopies[-rem_locus]
        df$Nloci = df$Nloci - length(rem_locus)
  }

  #Calculates missing % of data for all loci
  missingsum= sum(df$gen_data_int== -999)
  missingpercent= (missingsum / ((ncol(df$gen_data_int))*nrow(df$gen_data_int)))*100
  cat("This data is missing ", missingpercent ,"% of loci data \n")
  
  
  #Calculates missing % of data per locus
  xy=1
  for(col in seq(1, df$Nloci * df$ploidy, df$ploidy)){
    missingsumper= sum(df$gen_data_int[,col]== -999)+sum(df$gen_data_int[,col+1]== -999)
    missingpercentper= (missingsumper/(nrow(df$gen_data_int))*100)
    cat("Locus: ", loci_names[xy], " is missing ", missingpercentper, "% of locus data \n")
    xy= xy+1
  }
  

  names(df$Nallele) = df$loci_names

  return(df)

}

#Summary method for sgsObj class
summary.sgsObj <- function(x){
  cat("Number of individuals: ", x$Nind, "\n")
  cat("Number of categories: ", length(unique(x$groups)), "\n")
  cat("Number of loci: ",x$Nloci, "\n" )
  cat("Number of alleles per locus: ", x$Nallele, "\n")
  cat("Number of gene copies per loci: ",x$Ngenecopies ,"\n")
  xz=1
  for(col in seq(1, x$Nloci * x$ploidy, x$ploidy)){
    missingsumper= sum(x$gen_data_int[,col]== -999)+sum(x$gen_data_int[,col+1]== -999)
    missingpercentper= (missingsumper/(nrow(x$gen_data_int))*100)
    cat("Locus: ", loci_names[xz], " is missing ", missingpercentper, "% of locus data \n")
    xz= xz+1
  }
  
}




#################
##### Read spagedi input file and convert to sgsObj
###########

#' Read SPAGeDi input from file and convert to sgsObj object
#'
#' This function reads a SPAGeDi formatted text file (tab delimited text) and converts it to an sgsObj.
#'
#' @param path_to_spagedi_file A file path pointing to SPAGeDi input file. Must point to a tab delimited text file in following format in SPAGeDi manual (see details below).
#' @param missing_val Value indidicating missing data. Default is -999
#'
#' @return An object of the class sgsObj, which is used for futher sgs analyses.
#' @export
#'
#'@section Details:
#'
#'SPAGeDi (Spatial Pattern Analysis of Genetic Diversity) is a popular program for calculating spatial genetic structure, among other things. More information about SPAGeDi can be found at their hompage \href{http://ebe.ulb.ac.be/ebe/SPAGeDi.html}{here}.
#'
#'Data format should follow the style described in section 3.1 of the SPAGeDi manual.
#'
#'File should be a tab delimited text file.
#'
#'Briefly, the first line is a series of 6 numbers separated by tabs: number of individuals, number of categories, number of spatial coordinates, number of loci, number of digits used to code one allele, and ploidy.
#'
#'The second line indicates the distance intervals (which can be changed after inputting file).
#'
#'The third line lists the column names for individuals, categories, spatial coordinates, and loci.
#'
#'The fourth line begins the sample ID, category, spatial coordinate, and genotype data, with each line representing an individual.
#'
#'END is entered in the line after the last individual.
#'
#'\emph{Note:} This function as been tested only for diploid genetic data in 2d coordinate systems.
#'
#'
#' @examples
#'## Read in example data provided with package
#' path_to_example_data <- system.file("extdata", "obataua_spagedi.txt", package = "sgsR")
#'
#'## Use readSpagedi to read and convert data to sgsObj
#' dat <- readSpagedi(path_to_example_data, missing_val = "0")
#'
#' dat
#'
#'
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

    ## Fancy up column names for genotype data
    n_cols <- ncol(split_gen_data)
    new_col_names <- rep(NA, n_cols)

    new_col_names[seq(1, n_cols, 2)] <- paste(loci_names, "_A", sep = "")
    new_col_names[seq(2, n_cols, 2)] <- paste(loci_names, "_B", sep = "")

    colnames(split_gen_data) <- new_col_names


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





