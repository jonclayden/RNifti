#include <R_ext/Rdynload.h>
#include <Rcpp.h>

#include "niftilib/nifti1_io.h"
#include "lib/NiftiImage.h"

using namespace Rcpp;
using namespace RNifti;

typedef std::vector<float> float_vector;

RcppExport SEXP readNifti (SEXP _object, SEXP _internal)
{
BEGIN_RCPP
    NiftiImage image(_object);
    return image.toArrayOrPointer(as<bool>(_internal), "NIfTI image");
END_RCPP
}

RcppExport SEXP writeNifti (SEXP _image, SEXP _file, SEXP _datatype)
{
BEGIN_RCPP
    NiftiImage image(_image);
    image.toFile(as<std::string>(_file), as<std::string>(_datatype));
    return R_NilValue;
END_RCPP
}

RcppExport SEXP updateNifti (SEXP _image, SEXP _reference)
{
BEGIN_RCPP
    const NiftiImage reference(_reference);
    RObject image(_image);
    
    if (!reference.isNull())
    {
        NiftiImage updatedImage = reference;
        updatedImage.update(image);
        return updatedImage.toArray();
    }
    else
        return image;
END_RCPP
}

RcppExport SEXP dumpNifti (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, false);
    return image.headerToList();
END_RCPP
}

RcppExport SEXP getXform (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    const NiftiImage image(_image, false);
    const bool preferQuaternion = as<bool>(_preferQuaternion);
    
    mat44 xform = image.xform(preferQuaternion);
    NumericMatrix matrix(4,4);
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            matrix(i,j) = static_cast<double>(xform.m[i][j]);
    }
    
    return matrix;
END_RCPP
}

RcppExport SEXP setXform (SEXP _image, SEXP _matrix, SEXP _isQform)
{
BEGIN_RCPP
    NiftiImage image(_image);
    NumericMatrix matrix(_matrix);
    
    if (matrix.cols() != 4 || matrix.rows() != 4)
        throw std::runtime_error("Specified affine matrix does not have dimensions of 4x4");
    mat44 xform;
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            xform.m[i][j] = static_cast<float>(matrix(i,j));
    }
    
    int code = -1;
    if (!Rf_isNull(matrix.attr("code")))
        code = as<int>(matrix.attr("code"));
    
    if (!image.isNull())
    {
        if (as<bool>(_isQform))
        {
            image->qto_xyz = xform;
            image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
            nifti_mat44_to_quatern(image->qto_xyz, &image->quatern_b, &image->quatern_c, &image->quatern_d, &image->qoffset_x, &image->qoffset_y, &image->qoffset_z, NULL, NULL, NULL, &image->qfac);
            
            if (code >= 0)
                image->qform_code = code;
        }
        else
        {
            image->sto_xyz = xform;
            image->sto_ijk = nifti_mat44_inverse(image->sto_xyz);
            
            if (code >= 0)
                image->sform_code = code;
        }
    }
    
    if (image.isPersistent())
        return _image;
    else
        return image.toArray();
END_RCPP
}

RcppExport SEXP rescaleImage (SEXP _image, SEXP _scales)
{
BEGIN_RCPP
    const float_vector scales = as<float_vector>(_scales);
    const NiftiImage image(_image);
    
    NiftiImage newImage(nifti_copy_nim_info(image));
    newImage.rescale(scales);
    return newImage.toPointer("NIfTI image");
END_RCPP
}

RcppExport SEXP pointerToArray (SEXP _image)
{
BEGIN_RCPP
    NiftiImage image(_image);
    return image.toArray();
END_RCPP
}

extern "C" {

static R_CallMethodDef callMethods[] = {
    { "readNifti",      (DL_FUNC) &readNifti,       2 },
    { "writeNifti",     (DL_FUNC) &writeNifti,      3 },
    { "updateNifti",    (DL_FUNC) &updateNifti,     2 },
    { "dumpNifti",      (DL_FUNC) &dumpNifti,       1 },
    { "getXform",       (DL_FUNC) &getXform,        2 },
    { "setXform",       (DL_FUNC) &setXform,        3 },
    { "rescaleImage",   (DL_FUNC) &rescaleImage,    2 },
    { "pointerToArray", (DL_FUNC) &pointerToArray,  1 },
    { NULL, NULL, 0 }
};

void R_init_RNifti (DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    
    R_RegisterCCallable("RNifti",   "nii_make_new_header",  (DL_FUNC) &nifti_make_new_header);
    R_RegisterCCallable("RNifti",   "nii_make_new_nim",     (DL_FUNC) &nifti_make_new_nim);
    R_RegisterCCallable("RNifti",   "nii_convert_nhdr2nim", (DL_FUNC) &nifti_convert_nhdr2nim);
    R_RegisterCCallable("RNifti",   "nii_convert_nim2nhdr", (DL_FUNC) &nifti_convert_nim2nhdr);
    R_RegisterCCallable("RNifti",   "nii_copy_nim_info",    (DL_FUNC) &nifti_copy_nim_info);
    R_RegisterCCallable("RNifti",   "nii_copy_extensions",  (DL_FUNC) &nifti_copy_extensions);
    R_RegisterCCallable("RNifti",   "nii_image_unload",     (DL_FUNC) &nifti_image_unload);
    R_RegisterCCallable("RNifti",   "nii_image_free",       (DL_FUNC) &nifti_image_free);
    
    R_RegisterCCallable("RNifti",   "nii_datatype_sizes",   (DL_FUNC) &nifti_datatype_sizes);
    R_RegisterCCallable("RNifti",   "nii_datatype_string",  (DL_FUNC) &nifti_datatype_string);
    R_RegisterCCallable("RNifti",   "nii_units_string",     (DL_FUNC) &nifti_units_string);
    R_RegisterCCallable("RNifti",   "nii_get_volsize",      (DL_FUNC) &nifti_get_volsize);
    R_RegisterCCallable("RNifti",   "nii_update_dims_from_array", (DL_FUNC) &nifti_update_dims_from_array);
    
    R_RegisterCCallable("RNifti",   "nii_set_filenames",    (DL_FUNC) &nifti_set_filenames);
    R_RegisterCCallable("RNifti",   "nii_image_read",       (DL_FUNC) &nifti_image_read);
    R_RegisterCCallable("RNifti",   "nii_image_write",      (DL_FUNC) &nifti_image_write);
    
    R_RegisterCCallable("RNifti",   "nii_mat33_rownorm",    (DL_FUNC) &nifti_mat33_rownorm);
    R_RegisterCCallable("RNifti",   "nii_mat33_colnorm",    (DL_FUNC) &nifti_mat33_colnorm);
    R_RegisterCCallable("RNifti",   "nii_mat33_determ",     (DL_FUNC) &nifti_mat33_determ);
    R_RegisterCCallable("RNifti",   "nii_mat33_inverse",    (DL_FUNC) &nifti_mat33_inverse);
    R_RegisterCCallable("RNifti",   "nii_mat33_mul",        (DL_FUNC) &nifti_mat33_mul);
    R_RegisterCCallable("RNifti",   "nii_mat33_polar",      (DL_FUNC) &nifti_mat33_polar);
    R_RegisterCCallable("RNifti",   "nii_mat44_inverse",    (DL_FUNC) &nifti_mat44_inverse);
    R_RegisterCCallable("RNifti",   "nii_mat44_to_quatern", (DL_FUNC) &nifti_mat44_to_quatern);
    R_RegisterCCallable("RNifti",   "nii_quatern_to_mat44", (DL_FUNC) &nifti_quatern_to_mat44);
    
    R_RegisterCCallable("RNifti",   "nii_swap_2bytes",      (DL_FUNC) &nifti_swap_2bytes);
    R_RegisterCCallable("RNifti",   "nii_swap_4bytes",      (DL_FUNC) &nifti_swap_4bytes);
    R_RegisterCCallable("RNifti",   "nii_swap_8bytes",      (DL_FUNC) &nifti_swap_8bytes);
    R_RegisterCCallable("RNifti",   "nii_swap_16bytes",     (DL_FUNC) &nifti_swap_16bytes);
}

}
