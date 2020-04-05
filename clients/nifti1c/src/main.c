#include <R.h>

#include "RNifti.h"
#include "RNiftiAPI.h"

SEXP test (SEXP _path)
{
    nifti_image *image = nifti_image_read(CHAR(STRING_ELT(_path, 0)), 0);
    nifti_1_header header = nifti_convert_nim2nhdr(image);
    int status = disp_nifti_1_header(NULL, &header);
    nifti_image_free(image);
    return ScalarInteger(status);
}
