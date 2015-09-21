// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calcFijPairwiseCpp
std::vector<float> calcFijPairwiseCpp(List ref_gen, NumericMatrix alfreq1, NumericMatrix alfreq2, int Nloci, NumericVector Nallele, NumericVector Ngenecopies);
RcppExport SEXP sgsR_calcFijPairwiseCpp(SEXP ref_genSEXP, SEXP alfreq1SEXP, SEXP alfreq2SEXP, SEXP NlociSEXP, SEXP NalleleSEXP, SEXP NgenecopiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type ref_gen(ref_genSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alfreq1(alfreq1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alfreq2(alfreq2SEXP);
    Rcpp::traits::input_parameter< int >::type Nloci(NlociSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nallele(NalleleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ngenecopies(NgenecopiesSEXP);
    __result = Rcpp::wrap(calcFijPairwiseCpp(ref_gen, alfreq1, alfreq2, Nloci, Nallele, Ngenecopies));
    return __result;
END_RCPP
}
// calcFijPopCpp
NumericMatrix calcFijPopCpp(NumericMatrix Mcij, NumericMatrix Mdij, NumericVector distance_intervals, NumericMatrix genotype_data, List ref_gen, int Nloci, NumericVector Nallele, int Nind, NumericVector Ngenecopies, int Nperm, NumericVector x_coord, NumericVector y_coord);
RcppExport SEXP sgsR_calcFijPopCpp(SEXP McijSEXP, SEXP MdijSEXP, SEXP distance_intervalsSEXP, SEXP genotype_dataSEXP, SEXP ref_genSEXP, SEXP NlociSEXP, SEXP NalleleSEXP, SEXP NindSEXP, SEXP NgenecopiesSEXP, SEXP NpermSEXP, SEXP x_coordSEXP, SEXP y_coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Mcij(McijSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Mdij(MdijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_intervals(distance_intervalsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type genotype_data(genotype_dataSEXP);
    Rcpp::traits::input_parameter< List >::type ref_gen(ref_genSEXP);
    Rcpp::traits::input_parameter< int >::type Nloci(NlociSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nallele(NalleleSEXP);
    Rcpp::traits::input_parameter< int >::type Nind(NindSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ngenecopies(NgenecopiesSEXP);
    Rcpp::traits::input_parameter< int >::type Nperm(NpermSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_coord(x_coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_coord(y_coordSEXP);
    __result = Rcpp::wrap(calcFijPopCpp(Mcij, Mdij, distance_intervals, genotype_data, ref_gen, Nloci, Nallele, Nind, Ngenecopies, Nperm, x_coord, y_coord));
    return __result;
END_RCPP
}
// calcAlleleFreqPop
NumericVector calcAlleleFreqPop(NumericVector alleles_1, NumericVector alleles_2, int Nallele, int Ngenecopies);
RcppExport SEXP sgsR_calcAlleleFreqPop(SEXP alleles_1SEXP, SEXP alleles_2SEXP, SEXP NalleleSEXP, SEXP NgenecopiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type alleles_1(alleles_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alleles_2(alleles_2SEXP);
    Rcpp::traits::input_parameter< int >::type Nallele(NalleleSEXP);
    Rcpp::traits::input_parameter< int >::type Ngenecopies(NgenecopiesSEXP);
    __result = Rcpp::wrap(calcAlleleFreqPop(alleles_1, alleles_2, Nallele, Ngenecopies));
    return __result;
END_RCPP
}
// calcAlleleFreqCppInd
std::vector<float> calcAlleleFreqCppInd(int alleles_1, int alleles_2, int Nallele);
RcppExport SEXP sgsR_calcAlleleFreqCppInd(SEXP alleles_1SEXP, SEXP alleles_2SEXP, SEXP NalleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type alleles_1(alleles_1SEXP);
    Rcpp::traits::input_parameter< int >::type alleles_2(alleles_2SEXP);
    Rcpp::traits::input_parameter< int >::type Nallele(NalleleSEXP);
    __result = Rcpp::wrap(calcAlleleFreqCppInd(alleles_1, alleles_2, Nallele));
    return __result;
END_RCPP
}
// calcPairwiseDist
NumericMatrix calcPairwiseDist(NumericVector x, NumericVector y, int Nind);
RcppExport SEXP sgsR_calcPairwiseDist(SEXP xSEXP, SEXP ySEXP, SEXP NindSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type Nind(NindSEXP);
    __result = Rcpp::wrap(calcPairwiseDist(x, y, Nind));
    return __result;
END_RCPP
}
// findDIs
NumericMatrix findDIs(NumericMatrix Mdij, NumericVector distance_intervals, int Nind);
RcppExport SEXP sgsR_findDIs(SEXP MdijSEXP, SEXP distance_intervalsSEXP, SEXP NindSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Mdij(MdijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_intervals(distance_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type Nind(NindSEXP);
    __result = Rcpp::wrap(findDIs(Mdij, distance_intervals, Nind));
    return __result;
END_RCPP
}
// findEqualDIs
NumericVector findEqualDIs(NumericMatrix Mdij, NumericVector distance_intervals, int Nind);
RcppExport SEXP sgsR_findEqualDIs(SEXP MdijSEXP, SEXP distance_intervalsSEXP, SEXP NindSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Mdij(MdijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_intervals(distance_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type Nind(NindSEXP);
    __result = Rcpp::wrap(findEqualDIs(Mdij, distance_intervals, Nind));
    return __result;
END_RCPP
}
// summarizeDIs
NumericMatrix summarizeDIs(NumericMatrix Mdij, NumericMatrix Mcij, NumericVector distance_intervals, int Nind);
RcppExport SEXP sgsR_summarizeDIs(SEXP MdijSEXP, SEXP McijSEXP, SEXP distance_intervalsSEXP, SEXP NindSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Mdij(MdijSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Mcij(McijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_intervals(distance_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type Nind(NindSEXP);
    __result = Rcpp::wrap(summarizeDIs(Mdij, Mcij, distance_intervals, Nind));
    return __result;
END_RCPP
}
