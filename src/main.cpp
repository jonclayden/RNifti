#include <Rcpp.h>

#include "lib/NiftiImage.h"

using namespace Rcpp;

typedef std::vector<float> float_vector;

std::map<std::string,short> NiftiImage::DatatypeCodes = NiftiImage::buildDatatypeCodes();

RcppExport SEXP retrieveImage (SEXP _image)
{
BEGIN_RCPP
    NiftiImage image(_image);
    return image.toPointer("NIfTI image");
END_RCPP
}

RcppExport SEXP readNifti (SEXP _file, SEXP _internal)
{
BEGIN_RCPP
    NiftiImage image(_file);
    return image.toArrayOrPointer(as<bool>(_internal), "NIfTI image");
END_RCPP
}

RcppExport SEXP writeNifti (SEXP _image, SEXP _file, SEXP _datatype)
{
BEGIN_RCPP
    NiftiImage image(_image);
    std::string datatypeString = as<std::string>(_datatype);
    if (NiftiImage::DatatypeCodes.count(datatypeString) == 0)
    {
        std::ostringstream message;
        message << "Datatype \"" << datatypeString << "\" is not valid";
        Rf_warning(message.str().c_str());
        
        datatypeString = "auto";
    }
    image.toFile(as<std::string>(_file), NiftiImage::DatatypeCodes[datatypeString]);
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
