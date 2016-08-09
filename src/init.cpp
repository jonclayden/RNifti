#include <R_ext/Rdynload.h>
#include <Rcpp.h>

#include "niftilib/nifti1_io.h"

extern "C" {

void R_init_RNifti (DllInfo *info)
{
    R_RegisterCCallable("RNifti",   "nii_make_new_header",  (DL_FUNC) &nifti_make_new_header);
    R_RegisterCCallable("RNifti",   "nii_make_new_nim",     (DL_FUNC) &nifti_make_new_nim);
    R_RegisterCCallable("RNifti",   "nii_convert_nhdr2nim", (DL_FUNC) &nifti_convert_nhdr2nim);
    R_RegisterCCallable("RNifti",   "nii_convert_nim2nhdr", (DL_FUNC) &nifti_convert_nim2nhdr);
    R_RegisterCCallable("RNifti",   "nii_copy_nim_info",    (DL_FUNC) &nifti_copy_nim_info);
    R_RegisterCCallable("RNifti",   "nii_copy_extensions",  (DL_FUNC) &nifti_copy_extensions);
    R_RegisterCCallable("RNifti",   "nii_image_free",       (DL_FUNC) &nifti_image_free);
    
    R_RegisterCCallable("RNifti",   "nii_datatype_sizes",   (DL_FUNC) &nifti_datatype_sizes);
    R_RegisterCCallable("RNifti",   "nii_datatype_string",  (DL_FUNC) &nifti_datatype_string);
    R_RegisterCCallable("RNifti",   "nii_units_string",     (DL_FUNC) &nifti_units_string);
    R_RegisterCCallable("RNifti",   "nii_get_volsize",      (DL_FUNC) &nifti_get_volsize);
    R_RegisterCCallable("RNifti",   "nii_update_dims_from_array", (DL_FUNC) &nifti_update_dims_from_array);
    
    R_RegisterCCallable("RNifti",   "nii_set_filenames",    (DL_FUNC) &nifti_set_filenames);
    R_RegisterCCallable("RNifti",   "nii_image_read",       (DL_FUNC) &nifti_image_read);
    R_RegisterCCallable("RNifti",   "nii_image_write",      (DL_FUNC) &nifti_image_write);
    
    R_RegisterCCallable("RNifti",   "nii_mat33_determ",     (DL_FUNC) &nifti_mat33_determ);
    R_RegisterCCallable("RNifti",   "nii_mat33_inverse",    (DL_FUNC) &nifti_mat33_inverse);
    R_RegisterCCallable("RNifti",   "nii_mat33_mul",        (DL_FUNC) &nifti_mat33_mul);
    R_RegisterCCallable("RNifti",   "nii_mat33_polar",      (DL_FUNC) &nifti_mat33_polar);
    R_RegisterCCallable("RNifti",   "nii_mat44_inverse",    (DL_FUNC) &nifti_mat44_inverse);
    R_RegisterCCallable("RNifti",   "nii_mat44_to_quatern", (DL_FUNC) &nifti_mat44_to_quatern);
    R_RegisterCCallable("RNifti",   "nii_quatern_to_mat44", (DL_FUNC) &nifti_quatern_to_mat44);
}

}
