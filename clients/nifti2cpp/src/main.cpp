#include <Rcpp.h>

#define RNIFTI_NIFTILIB_VERSION 2

#include "RNifti.h"
#include "RNiftiAPI.h"

SEXP test (SEXP _path)
{
    RNifti::NiftiImage image(Rcpp::as<std::string>(_path));
    nifti_1_header header;
    nifti_convert_nim2n1hdr(image, &header);
    int status = disp_nifti_1_header(NULL, &header);
    return Rcpp::wrap(status);
}
