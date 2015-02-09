#include <Rcpp.h>

using namespace Rcpp;

// date_original
CharacterVector date_original(SEXP file);
RcppExport SEXP phenotyping_date_original(SEXP fileSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type file(fileSEXP );
        CharacterVector __result = date_original(file);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// count_outside
NumericVector count_outside(SEXP s_im, SEXP s_rad, SEXP s_dx, SEXP s_dy, SEXP s_thresh);
RcppExport SEXP phenotyping_count_outside(SEXP s_imSEXP, SEXP s_radSEXP, SEXP s_dxSEXP, SEXP s_dySEXP, SEXP s_threshSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type s_im(s_imSEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_rad(s_radSEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_dx(s_dxSEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_dy(s_dySEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_thresh(s_threshSEXP );
        NumericVector __result = count_outside(s_im, s_rad, s_dx, s_dy, s_thresh);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// round_binary
NumericMatrix round_binary(SEXP s_im);
RcppExport SEXP phenotyping_round_binary(SEXP s_imSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type s_im(s_imSEXP );
        NumericMatrix __result = round_binary(s_im);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// closest_point
NumericMatrix closest_point(SEXP s_im, SEXP s_points, SEXP s_nfeatures);
RcppExport SEXP phenotyping_closest_point(SEXP s_imSEXP, SEXP s_pointsSEXP, SEXP s_nfeaturesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type s_im(s_imSEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_points(s_pointsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type s_nfeatures(s_nfeaturesSEXP );
        NumericMatrix __result = closest_point(s_im, s_points, s_nfeatures);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// normalpostp
List normalpostp(SEXP x, SEXP mu, SEXP sigma, SEXP lambda);
RcppExport SEXP phenotyping_normalpostp(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type x(xSEXP );
        Rcpp::traits::input_parameter< SEXP >::type mu(muSEXP );
        Rcpp::traits::input_parameter< SEXP >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< SEXP >::type lambda(lambdaSEXP );
        List __result = normalpostp(x, mu, sigma, lambda);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
