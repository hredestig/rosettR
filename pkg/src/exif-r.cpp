#include <stdio.h>
#include "exif.h"
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector date_original(SEXP file) {
  try{
    Rcpp::CharacterVector cfile(file);
    Rcpp::CharacterVector res(cfile.size());
    
    for(int i = 0; i < cfile.size(); i++) {
	  FILE *fp = fopen(cfile(i), "rb");
      if(!fp) { 
		throw std::invalid_argument("can't open image file.\n"); 
      }
      fseek(fp, 0, SEEK_END);
      unsigned long fsize = ftell(fp);
      rewind(fp);
      unsigned char *buf = new unsigned char[fsize];
      if(fread(buf, 1, fsize, fp) != fsize) {
		throw std::invalid_argument("Can't read file.\n"); 
      }
      fclose(fp);
      // Parse exif
      try{
		if(fsize > 0.1) {
		  EXIFInfo result;
		  int retval;
		  retval = ParseEXIF(buf, fsize, result);
		  if(retval != 0) {
			::Rf_error("failed to read exif tag");
		  }
		  res(i) = result.dateTimeOriginal;
		} else {
		  res(i) = "NA";
		}
      } catch (...){
		return R_NilValue;
      }
    }
    return res;
  } catch(std::exception &ex) {
		forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("unknown error in date_original");
  }
  return(0);
}

