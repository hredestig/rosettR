#include <stdio.h>
#include <R.h>
#include <math.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

static double msqrtpi = 1 / sqrt(2 * 3.141592653589793);

// [[Rcpp::export]]
NumericVector count_outside(SEXP s_im, SEXP s_rad, SEXP s_dx, SEXP s_dy,
							  SEXP s_thresh) {
  // count pixels that are outside a circle with radius rad, shifted
  // dx and dy from center of image and intensity below thresh
  try{
	Rcpp::NumericMatrix im(s_im);
	Rcpp::NumericVector rad(s_rad);
	Rcpp::NumericVector thresh(s_thresh);
	Rcpp::NumericVector dx(s_dx);
	Rcpp::NumericVector dy(s_dy);
	Rcpp::NumericVector soutside(1);
	double r = 0;
	double c = 0;
	int nr = im.nrow();
	int nc = im.ncol();
	int outside = 0;
	// shift (0,0) to center of image
	double hnr = nr / 2;
	double hnc = nc / 2;
	double r2 = rad(0) * rad(0);

	for(int i = 0; i < nr; i++) {
	  for(int j = 0; j < nc; j++) {
		r = (((double)i + 1) - dx(0)) - hnr;
		c = (((double)j + 1)- dy(0)) - hnc;
		if((r * r + c * c) > r2) {
		  if(im(i,j) < thresh(0)) {
			outside++;
		  }
		}
	  }
	}
	soutside(0) = outside;
	return soutside;
  }catch(...) {
	::Rf_error("unknown error in count_outside");
  }
  return(0);
}

// [[Rcpp::export]]
NumericMatrix round_binary(SEXP s_im) {
  try{
	Rcpp::NumericMatrix im(s_im);
	double midpoint = 0.1;
	for(int i = 0; i < im.nrow(); i++) {
	  for(int j = 0; j < im.ncol(); j++) {
		if(im(i,j) < midpoint) {
		  im(i,j) = 0;
		} else {
		  im(i,j) = 1;
		}
	  }
	}
	return im;
  }catch(...) {
	::Rf_error("unknown error in round_binary");
  }
  return(0);
}

// [[Rcpp::export]]
NumericMatrix closest_point(SEXP s_im, SEXP s_points, SEXP s_nfeatures) {
  try{
	Rcpp::NumericMatrix im(s_im);				  // labeled image
	Rcpp::NumericMatrix points(s_points);
	Rcpp::NumericVector nfeatures(s_nfeatures);
	Rcpp::NumericMatrix distances((int)nfeatures(0), points.nrow());
	int nr = im.nrow();
	int nc = im.ncol();
	int current;
	int previous;
	int next;
	int row;
	int col;
	double this_row;
	double this_col;
	double xd;
	double yd;
	double new_distance;

	for(int i = 0; i < distances.nrow(); i++) {
	  for(int j = 0; j < distances.ncol(); j++) {
		distances(i,j) = nr;
	  }
	}

	for(int i = 0; i < points.nrow(); i++) {
	  row = (int)points(i,0);
	  col = (int)points(i,1);
	  current = (int)im(row, col);
	  if(current != 0) {
		distances(current - 1, i) = (double)0;
	  }
	}

	for(int i = 1; i < (nr - 1); i++) {
	  for(int j = 1; j < (nc - 1); j++) {
		current = (int)im(i,j);
		if(current != 0) {
		  previous = (int)im(i, j - 1);
		  next = (int)im(i, j + 1);
		  if(previous == 0 || next == 0) {
			for(int k = 0; k < points.nrow(); k++) {
			  this_row = (double)i;
			  this_col = (double)j;
			  xd = this_row - points(k,0);
			  yd = this_col - points(k,1);
			  new_distance = sqrt((xd * xd) + (yd * yd));
			  if(distances(current - 1, k) > new_distance) {
				distances(current - 1, k) = new_distance;
			  }
			}
		  }
		}
	  }
	}
	return distances;
  }catch(...) {
	::Rf_error("unknown error in closest_point");
  }
  return(0);
}

// [[Rcpp::export]]
List normalpostp(SEXP x, SEXP mu, SEXP sigma, SEXP lambda) {
  try{
	Rcpp::NumericVector cx(x);
	Rcpp::NumericVector cmu(mu);
	Rcpp::NumericVector csigma(sigma);
	Rcpp::NumericVector clambda(lambda);
	int N = cx.size();
	int K = cmu.size();
	Rcpp::NumericMatrix resp(N, K);
	Rcpp::NumericMatrix scaled_resp(N, K);
	Rcpp::NumericVector sresp(N);
	Rcpp::NumericVector ll(1);
	double v = 0;
	double u = 0;

	for(int k = 0; k < K; k++) {
	  for(int n = 0; n < N; n++) {
		v = cx(n) - cmu(k);
		resp(n,k) = clambda(k) * msqrtpi * (1 / csigma(k)) *
		  exp(-(v*v) / (2 * csigma(k) * csigma(k)));
	  }
	}

	ll(0) = 0;
	for(int n = 0; n < N; n++) {
	  u = 0;
	  for(int k = 0; k < K; k++) {
		u += resp(n,k);
	  }
	  ll(0) += log(u);
	}

	for(int n = 0; n < N; n++) {
	  sresp(n) = 0;
	  for(int k = 0; k < K; k++) {
		sresp(n) += resp(n,k);
	  }
	}

	for(int k = 0; k < K; k++) {
	  for(int n = 0; n < N; n++) {
		scaled_resp(n,k) = resp(n,k) / sresp(n);
		if(isnan(scaled_resp(n,k))) {
		  if(K == 1) {
			scaled_resp(n,k) = 1;
		  } else {
			scaled_resp(n,k) = 0;
		  }
		}
	  }
	}

	return Rcpp::List::create(Rcpp::Named("resp", resp),
							  Rcpp::Named("scaled_resp", scaled_resp),
							  Rcpp::Named("loglik", ll));
  }catch(...) {
	::Rf_error("unknown error in normalpostp");
  }
  return(0);
}

