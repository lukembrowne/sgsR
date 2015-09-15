### sgsR

*sgsR* is a package for calculating spatial genetic structure in R. The aim is to implement analyses, similar to those found in [SPAGeDi](http://ebe.ulb.ac.be/ebe/SPAGeDi.html) and [GenAlEx](http://biology-assets.anu.edu.au/GenAlEx/Welcome.html), that estimate the degree of spatial autocorrelation in genetic data.

Some key features of sgsR are:

1.  Calculating relatedness among individuals based on set distance intervals
2.  Conducting permutation tests
3.  Creating spatial autocorrelation plots
4.  Reading and converting from SPAGeDi data format

===============

sgsR is very much still in development, and certainly contains many bugs. If you're interested in contributing, have any comments or suggestions, please get in touch via github or email - <lukembrowne@gmail.com>

===============

To install sgsR, you must first make sure the package 'devtools' is installed. This will allow you to install sgsR directly from github.

``` r

install.packages("devtools")

devtools::install_github("lukembrowne/sgsR")
```

=================

Here's an example of a typical workflow...

-   Input data into the sgs data structure using createSgsObj()
-   Set desired distance intervals and number of permutations
-   Run SGS analysis with sgs(), currently configured to run the kinship coefficient of Loiselle et al. 1995
-   use plotSgs() to make an autorcorrelation plot

``` r

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
sgsObj = createSgsObj(sample_ids = dat$id, genotype_data = dat[, 4:(Nloci*2 + 3)],
                    ploidy = 2,
                    x_coords = dat$x, y_coords = dat$y)
attributes(sgsObj)
#> $names
#>  [1] "ids"         "x"           "y"           "ploidy"      "Nind"       
#>  [6] "Nloci"       "Ngenecopies" "gen_data"    "loci_names"  "gen_data_f" 
#> [11] "Nallele"    
#> 
#> $class
#> [1] "sgsObj"

# Display genetic data
head(sgsObj$gen_data)
#>   Loc1_A Loc1_B Loc2_A Loc2_B Loc3_A Loc3_B Loc4_A Loc4_B Loc5_A Loc5_B
#> 1      7      5      5      4     10      1      5      4      8      8
#> 2      8      5      4      2      9      1      1     10     10      6
#> 3      3      8      4      3      2     10      9      9      4      7
#> 4      9      7      2      4      2      8      2      2     10     10
#> 5      4      2      7     10      4      2      2      1     10      4
#> 6      3      9      2      5      2      3      9      8      2      7
#>   Loc6_A Loc6_B Loc7_A Loc7_B Loc8_A Loc8_B Loc9_A Loc9_B Loc10_A Loc10_B
#> 1      1      7     10      9      7      9      8      6       5       6
#> 2     10      2      2      8      3      8      4      9       8       4
#> 3      4      7      7      6      1      1      5      4       1       5
#> 4      6      6      4      9      2      4      6      1       6       6
#> 5      8      8      1      1      2      5      4     10       3       9
#> 6      9      7      6      3      5     10      3      3      10      10


## Run analysis
distance_intervals = seq(10, 110, 10) # Set distance intervals

out1 = sgs(sgsObj = sgsObj, distance_intervals = distance_intervals, nperm = 99)
#> Adding an aditional distance interval -- 128.0867 -- to encompass all pairwise distances.. 
#> Working on permutation: 0...

## Solid line is Fij estimate for each distance class
## Dashed lines are the 2.5 % and 97.5 % quantiles for permuted values
plotSgs(out1)
```

![](README-unnamed-chunk-3-1.png)

``` r

# Summary of information on distance classes
out1$DIsummary
#>                        [,1]      [,2]      [,3]      [,4]     [,5]
#> Distance class     0.000000   1.00000   2.00000   3.00000   4.0000
#> Max distance      10.000000  20.00000  30.00000  40.00000  50.0000
#> Average distance   6.241301  15.40722  25.35529  35.23155  45.0197
#> Number of pairs  155.000000 376.00000 568.00000 663.00000 751.0000
#>                       [,6]      [,7]     [,8]      [,9]    [,10]    [,11]
#> Distance class     5.00000   6.00000   7.0000   8.00000   9.0000  10.0000
#> Max distance      60.00000  70.00000  80.0000  90.00000 100.0000 110.0000
#> Average distance  54.93701  65.11901  74.8578  84.53017  94.3508 104.1537
#> Number of pairs  662.00000 635.00000 496.0000 342.00000 192.0000  78.0000
#>                     [,12]
#> Distance class    11.0000
#> Max distance     128.0867
#> Average distance 116.7789
#> Number of pairs   32.0000

# Summary of information on estimated Kinship coefficient for each distance class (columns)
out1$Fijsummary
#>                   [,1]         [,2]          [,3]         [,4]
#> Loc1_A    0.0237558428 -0.002462208 -0.0094959792  0.001924795
#> Loc2_A    0.0198626816 -0.002304020  0.0123678287  0.004540861
#> Loc3_A   -0.0003407009 -0.010649785 -0.0002491854  0.006998091
#> Loc4_A    0.0119032776  0.007029250  0.0009330290  0.002315798
#> Loc5_A   -0.0058055278 -0.006823549 -0.0037423435  0.010052198
#> Loc6_A    0.0018560255  0.014216469  0.0001606579 -0.005835506
#> Loc7_A   -0.0028061776 -0.007722525  0.0130401161 -0.009435172
#> Loc8_A    0.0050973310  0.002283696  0.0143047823 -0.010299228
#> Loc9_A    0.0038945405 -0.007067842  0.0017300559 -0.001095822
#> Loc10_A   0.0066441614  0.024944881  0.0010885546 -0.013090606
#> ALL LOCI  0.0064061447  0.001144432  0.0030137538 -0.001392459
#>                   [,5]          [,6]          [,7]          [,8]
#> Loc1_A   -7.709755e-05 -0.0028043154  0.0062882029 -0.0008117849
#> Loc2_A   -4.783158e-03 -0.0023113398 -0.0122727090 -0.0068304008
#> Loc3_A   -2.303472e-03  0.0010664653 -0.0027740463 -0.0039491039
#> Loc4_A   -4.733878e-03 -0.0007023494 -0.0078601046 -0.0024615531
#> Loc5_A   -2.774844e-03  0.0032171591 -0.0049831741  0.0003328249
#> Loc6_A   -1.347182e-03  0.0007636534  0.0095534744  0.0027207553
#> Loc7_A   -3.885350e-03 -0.0030436206  0.0061925938 -0.0036150196
#> Loc8_A    3.872400e-03 -0.0075903297 -0.0091410568  0.0058972407
#> Loc9_A   -7.741791e-03  0.0030272054  0.0062378077  0.0052760653
#> Loc10_A  -9.726908e-04 -0.0115862982  0.0011726123  0.0108754151
#> ALL LOCI -2.474702e-03 -0.0019963761 -0.0007586374  0.0007434437
#>                   [,9]         [,10]        [,11]         [,12]
#> Loc1_A   -0.0020166556  0.0108379526 -0.027779466 -2.790866e-03
#> Loc2_A    0.0105467476  0.0065459241  0.010221006  1.212296e-02
#> Loc3_A    0.0111853369  0.0068393745 -0.008687842  1.363918e-02
#> Loc4_A    0.0024122228  0.0155052962  0.007187929 -1.996369e-02
#> Loc5_A    0.0053088753  0.0012169111  0.001349859 -1.817272e-02
#> Loc6_A   -0.0193969738 -0.0063013234 -0.034363575 -7.194886e-05
#> Loc7_A    0.0153912622  0.0030012184 -0.014595323 -1.598545e-03
#> Loc8_A    0.0056812624  0.0025798029 -0.012774367  2.078128e-02
#> Loc9_A   -0.0004138907  0.0001198648  0.010273156 -1.187416e-02
#> Loc10_A  -0.0023475701 -0.0073651900  0.020770593 -7.281457e-04
#> ALL LOCI  0.0026350562  0.0032979837 -0.004839804 -8.656670e-04
```

=================
