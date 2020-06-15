#include <R_ext/Rdynload.h>
#include <Rcpp.h>

#include "niftilib/nifti1_io.h"
#include "niftilib/nifti2_io.h"
#include "RNifti/NiftiImage.h"

using namespace Rcpp;
using namespace RNifti;

typedef std::vector<int> int_vector;
typedef std::vector<NiftiImage::dim_t> dim_vector;
typedef std::vector<NiftiImage::pixdim_t> pixdim_vector;

inline bool isXformMatrix (const SEXP object)
{
    if (!Rf_isMatrix(object))
        return false;
    NumericMatrix matrix(object);
    return (matrix.cols() == 4 && matrix.rows() == 4);
}

inline unsigned char clip (const double &value)
{
    unsigned char result;
    if (value < 0.0)
        result = 0;
    else if (value > 1.0)
        result = 255;
    else
        result = RNifti::internal::roundEven(value * 255.0);
    return result;
}

RcppExport SEXP packRgb (SEXP _object, SEXP _channels, SEXP _maxValue)
{
BEGIN_RCPP
    const size_t length = size_t(Rf_length(_object));
    const int channels = as<int>(_channels);
    const size_t pixels = length / size_t(channels);
    const double maxValue = as<double>(_maxValue);
    
    if (pixels * channels != length)
    {
        std::ostringstream message;
        message << "Data length (" << length << ") is not a multiple of the number of channels (" << channels << ")";
        Rf_error(message.str().c_str());
    }
    
    NumericVector source(_object);
    IntegerVector result(pixels);
    rgba32_t rgba;
    for (size_t i=0; i<pixels; i++)
    {
        if (channels > 2)
        {
            for (int j=0; j<channels; j++)
                rgba.value.bytes[j] = clip(source[i + pixels*j] / maxValue);
            for (int j=channels; j<4; j++)
                rgba.value.bytes[j] = 0;
        }
        else
        {
            for (int j=0; j<3; j++)
                rgba.value.bytes[j] = clip(source[i] / maxValue);
            rgba.value.bytes[3] = clip(source[i + pixels] / maxValue);
        }
        result[i] = rgba.value.packed;
    }
    
    return result;
END_RCPP
}

RcppExport SEXP rgbToStrings (SEXP _object)
{
BEGIN_RCPP
    const NiftiImage image(_object, true, true);
    const NiftiImageData data = image.data();
    CharacterVector result(image.nVoxels());
    for (size_t i=0; i<image.nVoxels(); i++)
    {
        rgba32_t rgba;
        rgba.value.packed = int(data[i]);
        std::ostringstream value;
        value << "#" << std::hex << std::uppercase;
        for (int j=0; j<image.nChannels(); j++)
            value << std::setw(2) << std::setfill('0') << int(rgba.value.bytes[j]);
        result[i] = value.str();
    }
    return result;
END_RCPP
}

RcppExport SEXP unpackRgb (SEXP _object, SEXP _channels)
{
BEGIN_RCPP
    const NiftiImage image(_object, true, true);
    const NiftiImageData data = image.data();
    const int_vector channels = as<int_vector>(_channels);
    dim_vector dim = image.dim();
    dim.push_back(channels.size());
    
    const size_t len = image.nVoxels();
    RawVector result(len * channels.size());
    rgba32_t rgba;
    for (size_t i=0; i<len; i++)
    {
        rgba.value.packed = int(data[i]);
        for (int j=0; j<channels.size(); j++)
        {
            const int channel = channels[j] - 1;
            result[i + len * j] = int(rgba.value.bytes[channel]);
        }
    }
    result.attr("dim") = dim;
    return result;
END_RCPP
}

RcppExport SEXP asNifti (SEXP _image, SEXP _reference, SEXP _datatype, SEXP _internal)
{
BEGIN_RCPP
    const std::string datatype = as<std::string>(_datatype);
    const bool willChangeDatatype = (datatype != "auto");
    const int internal = as<int>(_internal);
    const bool usePointer = (internal == 1 || (internal == NA_LOGICAL && Rf_inherits(_image,"internalImage")) || willChangeDatatype);
    
    NiftiImage image;
    if (Rf_isVectorList(_reference) && Rf_length(_reference) < 36)
    {
        image = NiftiImage(_image);
        image.update(_reference);
    }
    else if (Rf_isNull(_reference))
        image = NiftiImage(_image);
    else
    {
        image = NiftiImage(_reference);
        image.update(_image);
    }
    
    if (willChangeDatatype)
        image.changeDatatype(datatype);
    return image.toArrayOrPointer(usePointer, "NIfTI image");
END_RCPP
}

RcppExport SEXP niftiVersion (SEXP _path)
{
BEGIN_RCPP
    int version = NiftiImage::fileVersion(as<std::string>(_path));
    return wrap(version);
END_RCPP
}

RcppExport SEXP readNifti (SEXP _object, SEXP _internal, SEXP _volumes)
{
BEGIN_RCPP
    if (Rf_isNull(_volumes))
    {
        const NiftiImage image(as<std::string>(_object));
        return image.toArrayOrPointer(as<bool>(_internal), "NIfTI image");
    }
    else
    {
        std::vector<int> volumes;
        IntegerVector volumesR(_volumes);
        for (int i=0; i<volumesR.length(); i++)
            volumes.push_back(volumesR[i] - 1);
        const NiftiImage image(as<std::string>(_object), volumes);
        return image.toArrayOrPointer(as<bool>(_internal), "NIfTI image");
    }
END_RCPP
}

RcppExport SEXP writeNifti (SEXP _image, SEXP _file, SEXP _datatype, SEXP _filetype)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    const std::string filetypeString = as<std::string>(_filetype);
    int filetype = NIFTI_FTYPE_NIFTI1_1;
    if (filetypeString == "analyze")
        filetype = NIFTI_FTYPE_ANALYZE;
    const std::pair<std::string,std::string> paths = image.toFile(as<std::string>(_file), as<std::string>(_datatype), filetype);
    return CharacterVector::create(Named("header")=paths.first, Named("image")=paths.second);
END_RCPP
}

RcppExport SEXP niftiHeader (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, false, true);
    return image.headerToList();
END_RCPP
}

RcppExport SEXP analyzeHeader (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, false, true);
    nifti_1_header header = nifti_convert_nim2nhdr(image);
    nifti_analyze75 *analyze = (nifti_analyze75 *) &header;
    List result;
    
    result["sizeof_hdr"] = analyze->sizeof_hdr;
    result["data_type"] = std::string(analyze->data_type, 10);
    result["db_name"] = std::string(analyze->db_name, 18);
    result["extents"] = analyze->extents;
    result["session_error"] = analyze->session_error;
    result["regular"] = std::string(&analyze->regular, 1);
    
    result["dim"] = std::vector<short>(analyze->dim, analyze->dim+8);
    result["datatype"] = analyze->datatype;
    result["bitpix"] = analyze->bitpix;
    result["pixdim"] = std::vector<float>(analyze->pixdim, analyze->pixdim+8);
    
    result["vox_offset"] = analyze->vox_offset;
    result["cal_max"] = analyze->cal_max;
    result["cal_min"] = analyze->cal_min;
    result["compressed"] = analyze->compressed;
    result["verified"] = analyze->verified;
    result["glmax"] = analyze->glmax;
    result["glmin"] = analyze->glmin;
    
    result["descrip"] = std::string(analyze->descrip, 80);
    result["aux_file"] = std::string(analyze->aux_file, 24);
    result["orient"] = int(analyze->orient);
    result["originator"] = std::string(analyze->originator, 10);
    
    // SPM and FSL use the originator field to store a coordinate origin
    short *origin = (short *) analyze->originator;
    result["origin"] = std::vector<short>(origin, origin+5);
    
    result["generated"] = std::string(analyze->generated, 10);
    result["scannum"] = std::string(analyze->scannum, 10);
    result["patient_id"] = std::string(analyze->patient_id, 10);
    result["exp_date"] = std::string(analyze->exp_date, 10);
    result["exp_time"] = std::string(analyze->exp_time, 10);
    
    result["views"] = analyze->views;
    result["vols_added"] = analyze->vols_added;
    result["start_field"] = analyze->start_field;
    result["field_skip"] = analyze->field_skip;
    result["omax"] = analyze->omax;
    result["omin"] = analyze->omin;
    result["smax"] = analyze->smax;
    result["smin"] = analyze->smin;
    
    result.attr("class") = CharacterVector::create("analyzeHeader");
    
    return result;
END_RCPP
}

RcppExport SEXP getXform (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    if (isXformMatrix(_image))
        return _image;
    else
    {
        const NiftiImage image(_image, false, true);
        const bool preferQuaternion = as<bool>(_preferQuaternion);
        NumericMatrix matrix = wrap(image.xform(preferQuaternion).matrix());
        if (image.isNull())
            matrix.attr("code") = 0;
        else
            matrix.attr("code") = ((preferQuaternion && image->qform_code > 0) || image->sform_code <= 0) ? image->qform_code : image->sform_code;
        return matrix;
    }
END_RCPP
}

RcppExport SEXP setXform (SEXP _image, SEXP _matrix, SEXP _isQform)
{
BEGIN_RCPP
    NumericMatrix matrix(_matrix);
    NiftiImage::Xform xform(_matrix);
    
    int code = -1;
    if (!Rf_isNull(matrix.attr("code")))
        code = as<int>(matrix.attr("code"));
    
    if (Rf_isVectorList(_image) && Rf_inherits(_image,"niftiHeader"))
    {
        // Header only
        List image(_image);
        if (MAYBE_SHARED(_image))
            image = Rf_duplicate(image);
        
        float qbcd[3], qxyz[3], dxyz[3], qfac;
        nifti_mat44_to_quatern(xform, &qbcd[0], &qbcd[1], &qbcd[2], &qxyz[0], &qxyz[1], &qxyz[2], &dxyz[0], &dxyz[1], &dxyz[2], &qfac);
        
        if (as<bool>(_isQform))
        {
            *REAL(image["quatern_b"]) = static_cast<double>(qbcd[0]);
            *REAL(image["quatern_c"]) = static_cast<double>(qbcd[1]);
            *REAL(image["quatern_d"]) = static_cast<double>(qbcd[2]);
            *REAL(image["qoffset_x"]) = static_cast<double>(qxyz[0]);
            *REAL(image["qoffset_y"]) = static_cast<double>(qxyz[1]);
            *REAL(image["qoffset_z"]) = static_cast<double>(qxyz[2]);
            REAL(image["pixdim"])[0] = static_cast<double>(qfac);
            if (code >= 0)
                *INTEGER(image["qform_code"]) = code;
        }
        else
        {
            for (int i=0; i<4; i++)
            {
                REAL(image["srow_x"])[i] = matrix(0,i);
                REAL(image["srow_y"])[i] = matrix(1,i);
                REAL(image["srow_z"])[i] = matrix(2,i);
            }
            if (code >= 0)
                *INTEGER(image["sform_code"]) = code;
        }
        
        const int dimensionality = INTEGER(image["dim"])[0];
        for (int i=0; i<std::min(3,dimensionality); i++)
            REAL(image["pixdim"])[i+1] = static_cast<double>(dxyz[i]);
        
        return image;
    }
    else
    {
        // From here, we assume we have a proper image
        NiftiImage image(_image);
        
        if (!image.isNull())
        {
            if (as<bool>(_isQform))
            {
                image.qform() = xform;
                if (code >= 0)
                    image->qform_code = code;
            }
            else
            {
                image.sform() = xform;
                if (code >= 0)
                    image->sform_code = code;
            }
        }
    
        // If the image was copied above it will have been marked nonpersistent
        return image.toArrayOrPointer(Rf_inherits(_image,"internalImage"), "NIfTI image");
    }
END_RCPP
}

RcppExport SEXP getOrientation (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    std::string orientation;
    
    if (isXformMatrix(_image))
        orientation = NiftiImage::Xform(_image).orientation();
    else
    {
        const NiftiImage image(_image, false, true);
        orientation = image.xform(as<bool>(_preferQuaternion)).orientation();
    }
    
    return wrap(orientation);
END_RCPP
}

RcppExport SEXP setOrientation (SEXP _image, SEXP _axes)
{
BEGIN_RCPP
    if (isXformMatrix(_image))
    {
        // Create an empty image for temporary purposes
        nifti_image *ptr = nifti_make_new_nim(NULL, DT_UNSIGNED_CHAR, 0);
        NiftiImage image(ptr);
        
        // Set the qform matrix
        image.qform() = NiftiImage::Xform(_image);
        image->qform_code = 2;
        
        image.reorient(as<std::string>(_axes));
        return image.qform().matrix();
    }
    else
    {
        NiftiImage image(_image);
        image.reorient(as<std::string>(_axes));
        return image.toArrayOrPointer(Rf_inherits(_image,"internalImage"), "NIfTI image");
    }
END_RCPP
}

RcppExport SEXP getRotation (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    if (isXformMatrix(_image))
        return NiftiImage::Xform(_image).rotation();
    else
    {
        const NiftiImage image(_image, false, true);
        return image.xform(as<bool>(_preferQuaternion)).rotation();
    }
END_RCPP
}

RcppExport SEXP getAddresses (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    if (image.isNull())
        return R_NilValue;
    else
    {
        std::ostringstream imageString, dataString;
        imageString << (const nifti_image *) image;
        dataString << image->data;
        return CharacterVector::create(Named("image")=imageString.str(), Named("data")=dataString.str());
    }
END_RCPP
}

RcppExport SEXP hasData (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    return wrap(image->data != NULL);
END_RCPP
}

RcppExport SEXP indexVector (SEXP _image, SEXP _indices)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    if (image.isNull())
        Rf_error("Cannot index into a NULL image");
    else if (image->data == NULL)
        return LogicalVector(Rf_length(_indices), NA_LOGICAL);
    else
    {
        const IntegerVector indices(_indices);
        const NiftiImageData data = image.data();
        if (data.isComplex())
        {
            ComplexVector result(indices.length());
            const Rcomplex naValue = { NA_REAL, NA_REAL };
            for (int i=0; i<indices.length(); i++)
                result[i] = (size_t(indices[i]) > data.size() ? naValue : data[indices[i] - 1]);
            return result;
        }
        else if (data.isFloatingPoint() || data.isScaled())
        {
            NumericVector result(indices.length());
            for (int i=0; i<indices.length(); i++)
                result[i] = (size_t(indices[i]) > data.size() ? NA_REAL : data[indices[i] - 1]);
            return result;
        }
        else
        {
            IntegerVector result(indices.length());
            for (int i=0; i<indices.length(); i++)
                result[i] = (size_t(indices[i]) > data.size() ? NA_INTEGER : data[indices[i] - 1]);
            return result;
        }
    }
END_RCPP
}

RcppExport SEXP indexList (SEXP _image, SEXP _indices)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    if (image.isNull())
        Rf_error("Cannot index into a NULL image");
    else if (image->data == NULL)
        return LogicalVector(1, NA_LOGICAL);
    else
    {
        const List indices(_indices);
        const std::vector<int> dim = image.dim();
        const int nDims = indices.length();
        std::vector<size_t> strides(nDims);
        std::vector<dim_vector> locs(nDims);
        int_vector sizes(nDims);
        std::vector<size_t> cumulativeSizes(nDims);
        size_t count = 1;
        for (int i=0; i<nDims; i++)
        {
            strides[i] = (i == 0 ? 1 : strides[i-1] * dim[i-1]);
            locs[i] = as<dim_vector>(indices[i]);
            sizes[i] = locs[i].size();
            cumulativeSizes[i] = (i == 0 ? 1 : cumulativeSizes[i-1] * sizes[i-1]);
            count *= sizes[i];
        }
        
        const NiftiImageData data = image.data();
        if (data.isComplex())
        {
            ComplexVector result(count);
            const Rcomplex naValue = { NA_REAL, NA_REAL };
            for (size_t j=0; j<count; j++)
            {
                size_t loc = 0;
                for (int i=0; i<nDims; i++)
                    loc += (locs[i][(j / cumulativeSizes[i]) % sizes[i]] - 1) * strides[i];
                result[j] = (loc >= data.size() ? naValue : data[loc]);
            }
            return result;
        }
        else if (data.isFloatingPoint() || data.isScaled())
        {
            NumericVector result(count);
            for (size_t j=0; j<count; j++)
            {
                size_t loc = 0;
                for (int i=0; i<nDims; i++)
                    loc += (locs[i][(j / cumulativeSizes[i]) % sizes[i]] - 1) * strides[i];
                result[j] = (loc >= data.size() ? NA_REAL : data[loc]);
            }
            return result;
        }
        else
        {
            IntegerVector result(count);
            for (size_t j=0; j<count; j++)
            {
                size_t loc = 0;
                for (int i=0; i<nDims; i++)
                    loc += (locs[i][(j / cumulativeSizes[i]) % sizes[i]] - 1) * strides[i];
                result[j] = (loc >= data.size() ? NA_INTEGER : data[loc]);
            }
            return result;
        }
    }
END_RCPP
}

RcppExport SEXP rescaleImage (SEXP _image, SEXP _scales)
{
BEGIN_RCPP
    NiftiImage image(_image);
    image.rescale(as<pixdim_vector>(_scales));
    return image.toPointer("NIfTI image");
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
    { "packRgb",        (DL_FUNC) &packRgb,         3 },
    { "rgbToStrings",   (DL_FUNC) &rgbToStrings,    1 },
    { "unpackRgb",      (DL_FUNC) &unpackRgb,       2 },
    { "asNifti",        (DL_FUNC) &asNifti,         4 },
    { "niftiVersion",   (DL_FUNC) &niftiVersion,    1 },
    { "readNifti",      (DL_FUNC) &readNifti,       3 },
    { "writeNifti",     (DL_FUNC) &writeNifti,      4 },
    { "niftiHeader",    (DL_FUNC) &niftiHeader,     1 },
    { "analyzeHeader",  (DL_FUNC) &analyzeHeader,   1 },
    { "getXform",       (DL_FUNC) &getXform,        2 },
    { "setXform",       (DL_FUNC) &setXform,        3 },
    { "getOrientation", (DL_FUNC) &getOrientation,  2 },
    { "setOrientation", (DL_FUNC) &setOrientation,  2 },
    { "getRotation",    (DL_FUNC) &getRotation,     2 },
    { "getAddresses",   (DL_FUNC) &getAddresses,    1 },
    { "hasData",        (DL_FUNC) &hasData,         1 },
    { "indexVector",    (DL_FUNC) &indexVector,     2 },
    { "indexList",      (DL_FUNC) &indexList,       2 },
    { "rescaleImage",   (DL_FUNC) &rescaleImage,    2 },
    { "pointerToArray", (DL_FUNC) &pointerToArray,  1 },
    { NULL, NULL, 0 }
};

void R_init_RNifti (DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    
    R_RegisterCCallable("RNifti", "nii_datatype_string", (DL_FUNC) &nifti_datatype_string);
    R_RegisterCCallable("RNifti", "nii_units_string", (DL_FUNC) &nifti_units_string);
    R_RegisterCCallable("RNifti", "nii_intent_string", (DL_FUNC) &nifti_intent_string);
    R_RegisterCCallable("RNifti", "nii_xform_string", (DL_FUNC) &nifti_xform_string);
    R_RegisterCCallable("RNifti", "nii_slice_string", (DL_FUNC) &nifti_slice_string);
    R_RegisterCCallable("RNifti", "nii_orientation_string", (DL_FUNC) &nifti_orientation_string);
    R_RegisterCCallable("RNifti", "nii_is_inttype", (DL_FUNC) &nifti_is_inttype);
    R_RegisterCCallable("RNifti", "nii_mat44_inverse", (DL_FUNC) &nifti_mat44_inverse);
    R_RegisterCCallable("RNifti", "nii_mat33_inverse", (DL_FUNC) &nifti_mat33_inverse);
    R_RegisterCCallable("RNifti", "nii_mat33_polar", (DL_FUNC) &nifti_mat33_polar);
    R_RegisterCCallable("RNifti", "nii_mat33_rownorm", (DL_FUNC) &nifti_mat33_rownorm);
    R_RegisterCCallable("RNifti", "nii_mat33_colnorm", (DL_FUNC) &nifti_mat33_colnorm);
    R_RegisterCCallable("RNifti", "nii_mat33_determ", (DL_FUNC) &nifti_mat33_determ);
    R_RegisterCCallable("RNifti", "nii_mat33_mul", (DL_FUNC) &nifti_mat33_mul);
    R_RegisterCCallable("RNifti", "nii_swap_2bytes", (DL_FUNC) &nifti_swap_2bytes);
    R_RegisterCCallable("RNifti", "nii_swap_4bytes", (DL_FUNC) &nifti_swap_4bytes);
    R_RegisterCCallable("RNifti", "nii_swap_8bytes", (DL_FUNC) &nifti_swap_8bytes);
    R_RegisterCCallable("RNifti", "nii_swap_16bytes", (DL_FUNC) &nifti_swap_16bytes);
    R_RegisterCCallable("RNifti", "nii_swap_Nbytes", (DL_FUNC) &nifti_swap_Nbytes);
    R_RegisterCCallable("RNifti", "nii_datatype_is_valid", (DL_FUNC) &nifti_datatype_is_valid);
    R_RegisterCCallable("RNifti", "nii_datatype_from_string", (DL_FUNC) &nifti_datatype_from_string);
    R_RegisterCCallable("RNifti", "nii_datatype_to_string", (DL_FUNC) &nifti_datatype_to_string);
    R_RegisterCCallable("RNifti", "nii_get_filesize", (DL_FUNC) &nifti_get_filesize);
    R_RegisterCCallable("RNifti", "swap_nii_header", (DL_FUNC) &swap_nifti_header);
    R_RegisterCCallable("RNifti", "old_swap_nii_header", (DL_FUNC) &old_swap_nifti_header);
    R_RegisterCCallable("RNifti", "nii_swap_as_analyze", (DL_FUNC) &nifti_swap_as_analyze);
    R_RegisterCCallable("RNifti", "nii_image_read_bricks", (DL_FUNC) &nifti_image_read_bricks);
    R_RegisterCCallable("RNifti", "nii_image_load_bricks", (DL_FUNC) &nifti_image_load_bricks);
    R_RegisterCCallable("RNifti", "nii_free_NBL", (DL_FUNC) &nifti_free_NBL);
    R_RegisterCCallable("RNifti", "nii_image_read", (DL_FUNC) &nifti_image_read);
    R_RegisterCCallable("RNifti", "nii_image_load", (DL_FUNC) &nifti_image_load);
    R_RegisterCCallable("RNifti", "nii_image_unload", (DL_FUNC) &nifti_image_unload);
    R_RegisterCCallable("RNifti", "nii_image_free", (DL_FUNC) &nifti_image_free);
    R_RegisterCCallable("RNifti", "nii_read_collapsed_image", (DL_FUNC) &nifti_read_collapsed_image);
    R_RegisterCCallable("RNifti", "nii_read_subregion_image", (DL_FUNC) &nifti_read_subregion_image);
    R_RegisterCCallable("RNifti", "nii_image_write", (DL_FUNC) &nifti_image_write);
    R_RegisterCCallable("RNifti", "nii_image_write_bricks", (DL_FUNC) &nifti_image_write_bricks);
    R_RegisterCCallable("RNifti", "nii_image_infodump", (DL_FUNC) &nifti_image_infodump);
    R_RegisterCCallable("RNifti", "nii_disp_lib_hist", (DL_FUNC) &nifti_disp_lib_hist);
    R_RegisterCCallable("RNifti", "nii_disp_lib_version", (DL_FUNC) &nifti_disp_lib_version);
    R_RegisterCCallable("RNifti", "nii_disp_matrix_orient", (DL_FUNC) &nifti_disp_matrix_orient);
    R_RegisterCCallable("RNifti", "nii_disp_type_list", (DL_FUNC) &nifti_disp_type_list);
    R_RegisterCCallable("RNifti", "nii_image_to_ascii", (DL_FUNC) &nifti_image_to_ascii);
    R_RegisterCCallable("RNifti", "nii_image_from_ascii", (DL_FUNC) &nifti_image_from_ascii);
    R_RegisterCCallable("RNifti", "nii_get_volsize", (DL_FUNC) &nifti_get_volsize);
    R_RegisterCCallable("RNifti", "nii_set_filenames", (DL_FUNC) &nifti_set_filenames);
    R_RegisterCCallable("RNifti", "nii_makehdrname", (DL_FUNC) &nifti_makehdrname);
    R_RegisterCCallable("RNifti", "nii_makeimgname", (DL_FUNC) &nifti_makeimgname);
    R_RegisterCCallable("RNifti", "is_nii_file", (DL_FUNC) &is_nifti_file);
    R_RegisterCCallable("RNifti", "nii_find_file_extension", (DL_FUNC) &nifti_find_file_extension);
    R_RegisterCCallable("RNifti", "nii_is_complete_filename", (DL_FUNC) &nifti_is_complete_filename);
    R_RegisterCCallable("RNifti", "nii_validfilename", (DL_FUNC) &nifti_validfilename);
    R_RegisterCCallable("RNifti", "disp_nii_1_header", (DL_FUNC) &disp_nifti_1_header);
    R_RegisterCCallable("RNifti", "nii_set_debug_level", (DL_FUNC) &nifti_set_debug_level);
    R_RegisterCCallable("RNifti", "nii_set_skip_blank_ext", (DL_FUNC) &nifti_set_skip_blank_ext);
    R_RegisterCCallable("RNifti", "nii_set_allow_upper_fext", (DL_FUNC) &nifti_set_allow_upper_fext);
    R_RegisterCCallable("RNifti", "valid_nii_brick_list", (DL_FUNC) &valid_nifti_brick_list);
    R_RegisterCCallable("RNifti", "nii_image_open", (DL_FUNC) &nifti_image_open);
    R_RegisterCCallable("RNifti", "nii_image_write_hdr_img", (DL_FUNC) &nifti_image_write_hdr_img);
    R_RegisterCCallable("RNifti", "nii_image_write_hdr_img2", (DL_FUNC) &nifti_image_write_hdr_img2);
    R_RegisterCCallable("RNifti", "nii_read_buffer", (DL_FUNC) &nifti_read_buffer);
    R_RegisterCCallable("RNifti", "nii_write_all_data", (DL_FUNC) &nifti_write_all_data);
    R_RegisterCCallable("RNifti", "nii_write_buffer", (DL_FUNC) &nifti_write_buffer);
    R_RegisterCCallable("RNifti", "nii_read_ascii_image", (DL_FUNC) &nifti_read_ascii_image);
    R_RegisterCCallable("RNifti", "nii_write_ascii_image", (DL_FUNC) &nifti_write_ascii_image);
    R_RegisterCCallable("RNifti", "nii_datatype_sizes", (DL_FUNC) &nifti_datatype_sizes);
    R_RegisterCCallable("RNifti", "nii_mat44_to_quatern", (DL_FUNC) &nifti_mat44_to_quatern);
    R_RegisterCCallable("RNifti", "nii_quatern_to_mat44", (DL_FUNC) &nifti_quatern_to_mat44);
    R_RegisterCCallable("RNifti", "nii_make_orthog_mat44", (DL_FUNC) &nifti_make_orthog_mat44);
    R_RegisterCCallable("RNifti", "nii_short_order", (DL_FUNC) &nifti_short_order);
    R_RegisterCCallable("RNifti", "nii_mat44_to_orientation", (DL_FUNC) &nifti_mat44_to_orientation);
    R_RegisterCCallable("RNifti", "nii_findhdrname", (DL_FUNC) &nifti_findhdrname);
    R_RegisterCCallable("RNifti", "nii_findimgname", (DL_FUNC) &nifti_findimgname);
    R_RegisterCCallable("RNifti", "nii_is_gzfile", (DL_FUNC) &nifti_is_gzfile);
    R_RegisterCCallable("RNifti", "nii_makebasename", (DL_FUNC) &nifti_makebasename);
    R_RegisterCCallable("RNifti", "nii_convert_nim2nhdr", (DL_FUNC) &nifti_convert_nim2nhdr);
    R_RegisterCCallable("RNifti", "nii_make_new_header", (DL_FUNC) &nifti_make_new_header);
    R_RegisterCCallable("RNifti", "nii_read_header", (DL_FUNC) &nifti_read_header);
    R_RegisterCCallable("RNifti", "nii_copy_nim_info", (DL_FUNC) &nifti_copy_nim_info);
    R_RegisterCCallable("RNifti", "nii_make_new_nim", (DL_FUNC) &nifti_make_new_nim);
    R_RegisterCCallable("RNifti", "nii_simple_init_nim", (DL_FUNC) &nifti_simple_init_nim);
    R_RegisterCCallable("RNifti", "nii_convert_nhdr2nim", (DL_FUNC) &nifti_convert_nhdr2nim);
    R_RegisterCCallable("RNifti", "nii_hdr_looks_good", (DL_FUNC) &nifti_hdr_looks_good);
    R_RegisterCCallable("RNifti", "nii_is_valid_datatype", (DL_FUNC) &nifti_is_valid_datatype);
    R_RegisterCCallable("RNifti", "nii_is_valid_ecode", (DL_FUNC) &nifti_is_valid_ecode);
    R_RegisterCCallable("RNifti", "nii_nim_is_valid", (DL_FUNC) &nifti_nim_is_valid);
    R_RegisterCCallable("RNifti", "nii_nim_has_valid_dims", (DL_FUNC) &nifti_nim_has_valid_dims);
    R_RegisterCCallable("RNifti", "is_valid_nii_type", (DL_FUNC) &is_valid_nifti_type);
    R_RegisterCCallable("RNifti", "nii_test_datatype_sizes", (DL_FUNC) &nifti_test_datatype_sizes);
    R_RegisterCCallable("RNifti", "nii_type_and_names_match", (DL_FUNC) &nifti_type_and_names_match);
    R_RegisterCCallable("RNifti", "nii_update_dims_from_array", (DL_FUNC) &nifti_update_dims_from_array);
    R_RegisterCCallable("RNifti", "nii_set_iname_offset", (DL_FUNC) &nifti_set_iname_offset);
    R_RegisterCCallable("RNifti", "nii_set_type_from_names", (DL_FUNC) &nifti_set_type_from_names);
    R_RegisterCCallable("RNifti", "nii_add_extension", (DL_FUNC) &nifti_add_extension);
    R_RegisterCCallable("RNifti", "nii_compiled_with_zlib", (DL_FUNC) &nifti_compiled_with_zlib);
    R_RegisterCCallable("RNifti", "nii_copy_extensions", (DL_FUNC) &nifti_copy_extensions);
    R_RegisterCCallable("RNifti", "nii_free_extensions", (DL_FUNC) &nifti_free_extensions);
    R_RegisterCCallable("RNifti", "nii_get_intlist", (DL_FUNC) &nifti_get_intlist);
    R_RegisterCCallable("RNifti", "nii_strdup", (DL_FUNC) &nifti_strdup);
    R_RegisterCCallable("RNifti", "valid_nii_extensions", (DL_FUNC) &valid_nifti_extensions);
    R_RegisterCCallable("RNifti", "nii_mat44_mul", (DL_FUNC) &nifti_mat44_mul);
    R_RegisterCCallable("RNifti", "nii_dmat44_inverse", (DL_FUNC) &nifti_dmat44_inverse);
    R_RegisterCCallable("RNifti", "nii_mat44_to_dmat44", (DL_FUNC) &nifti_mat44_to_dmat44);
    R_RegisterCCallable("RNifti", "nii_dmat44_to_mat44", (DL_FUNC) &nifti_dmat44_to_mat44);
    R_RegisterCCallable("RNifti", "nii_dmat44_mul", (DL_FUNC) &nifti_dmat44_mul);
    R_RegisterCCallable("RNifti", "nii_dmat33_inverse", (DL_FUNC) &nifti_dmat33_inverse);
    R_RegisterCCallable("RNifti", "nii_dmat33_polar", (DL_FUNC) &nifti_dmat33_polar);
    R_RegisterCCallable("RNifti", "nii_dmat33_rownorm", (DL_FUNC) &nifti_dmat33_rownorm);
    R_RegisterCCallable("RNifti", "nii_dmat33_colnorm", (DL_FUNC) &nifti_dmat33_colnorm);
    R_RegisterCCallable("RNifti", "nii_dmat33_determ", (DL_FUNC) &nifti_dmat33_determ);
    R_RegisterCCallable("RNifti", "nii_dmat33_mul", (DL_FUNC) &nifti_dmat33_mul);
    R_RegisterCCallable("RNifti", "nii_header_version", (DL_FUNC) &nifti_header_version);
    R_RegisterCCallable("RNifti", "nii_swap_as_nifti1", (DL_FUNC) &nifti_swap_as_nifti1);
    R_RegisterCCallable("RNifti", "nii_swap_as_nifti2", (DL_FUNC) &nifti_swap_as_nifti2);
    R_RegisterCCallable("RNifti", "nii2_image_read_bricks", (DL_FUNC) &nifti2_image_read_bricks);
    R_RegisterCCallable("RNifti", "nii2_image_load_bricks", (DL_FUNC) &nifti2_image_load_bricks);
    R_RegisterCCallable("RNifti", "nii2_free_NBL", (DL_FUNC) &nifti2_free_NBL);
    R_RegisterCCallable("RNifti", "nii2_image_read", (DL_FUNC) &nifti2_image_read);
    R_RegisterCCallable("RNifti", "nii2_image_load", (DL_FUNC) &nifti2_image_load);
    R_RegisterCCallable("RNifti", "nii2_image_unload", (DL_FUNC) &nifti2_image_unload);
    R_RegisterCCallable("RNifti", "nii2_image_free", (DL_FUNC) &nifti2_image_free);
    R_RegisterCCallable("RNifti", "nii2_read_collapsed_image", (DL_FUNC) &nifti2_read_collapsed_image);
    R_RegisterCCallable("RNifti", "nii2_read_subregion_image", (DL_FUNC) &nifti2_read_subregion_image);
    R_RegisterCCallable("RNifti", "nii2_image_write", (DL_FUNC) &nifti2_image_write);
    R_RegisterCCallable("RNifti", "nii2_image_write_bricks", (DL_FUNC) &nifti2_image_write_bricks);
    R_RegisterCCallable("RNifti", "nii2_image_infodump", (DL_FUNC) &nifti2_image_infodump);
    R_RegisterCCallable("RNifti", "nii2_disp_matrix_orient", (DL_FUNC) &nifti2_disp_matrix_orient);
    R_RegisterCCallable("RNifti", "nii2_image_to_ascii", (DL_FUNC) &nifti2_image_to_ascii);
    R_RegisterCCallable("RNifti", "nii2_image_from_ascii", (DL_FUNC) &nifti2_image_from_ascii);
    R_RegisterCCallable("RNifti", "nii2_get_volsize", (DL_FUNC) &nifti2_get_volsize);
    R_RegisterCCallable("RNifti", "nii2_set_filenames", (DL_FUNC) &nifti2_set_filenames);
    R_RegisterCCallable("RNifti", "disp_nii_2_header", (DL_FUNC) &disp_nifti_2_header);
    R_RegisterCCallable("RNifti", "nii_get_alter_cifti", (DL_FUNC) &nifti_get_alter_cifti);
    R_RegisterCCallable("RNifti", "nii_set_alter_cifti", (DL_FUNC) &nifti_set_alter_cifti);
    R_RegisterCCallable("RNifti", "nii_alter_cifti_dims", (DL_FUNC) &nifti_alter_cifti_dims);
    R_RegisterCCallable("RNifti", "valid_nii2_brick_list", (DL_FUNC) &valid_nifti2_brick_list);
    R_RegisterCCallable("RNifti", "nii2_image_open", (DL_FUNC) &nifti2_image_open);
    R_RegisterCCallable("RNifti", "nii2_image_write_hdr_img", (DL_FUNC) &nifti2_image_write_hdr_img);
    R_RegisterCCallable("RNifti", "nii2_image_write_hdr_img2", (DL_FUNC) &nifti2_image_write_hdr_img2);
    R_RegisterCCallable("RNifti", "nii2_read_buffer", (DL_FUNC) &nifti2_read_buffer);
    R_RegisterCCallable("RNifti", "nii2_write_all_data", (DL_FUNC) &nifti2_write_all_data);
    R_RegisterCCallable("RNifti", "nii2_write_buffer", (DL_FUNC) &nifti2_write_buffer);
    R_RegisterCCallable("RNifti", "nii2_read_ascii_image", (DL_FUNC) &nifti2_read_ascii_image);
    R_RegisterCCallable("RNifti", "nii2_write_ascii_image", (DL_FUNC) &nifti2_write_ascii_image);
    R_RegisterCCallable("RNifti", "nii_dmat44_to_quatern", (DL_FUNC) &nifti_dmat44_to_quatern);
    R_RegisterCCallable("RNifti", "nii_quatern_to_dmat44", (DL_FUNC) &nifti_quatern_to_dmat44);
    R_RegisterCCallable("RNifti", "nii_make_orthog_dmat44", (DL_FUNC) &nifti_make_orthog_dmat44);
    R_RegisterCCallable("RNifti", "nii_dmat44_to_orientation", (DL_FUNC) &nifti_dmat44_to_orientation);
    R_RegisterCCallable("RNifti", "nii_convert_nim2n1hdr", (DL_FUNC) &nifti_convert_nim2n1hdr);
    R_RegisterCCallable("RNifti", "nii_convert_nim2n2hdr", (DL_FUNC) &nifti_convert_nim2n2hdr);
    R_RegisterCCallable("RNifti", "nii_make_new_n1_header", (DL_FUNC) &nifti_make_new_n1_header);
    R_RegisterCCallable("RNifti", "nii_make_new_n2_header", (DL_FUNC) &nifti_make_new_n2_header);
    R_RegisterCCallable("RNifti", "nii2_read_header", (DL_FUNC) &nifti2_read_header);
    R_RegisterCCallable("RNifti", "nii_read_n1_hdr", (DL_FUNC) &nifti_read_n1_hdr);
    R_RegisterCCallable("RNifti", "nii_read_n2_hdr", (DL_FUNC) &nifti_read_n2_hdr);
    R_RegisterCCallable("RNifti", "nii2_copy_nim_info", (DL_FUNC) &nifti2_copy_nim_info);
    R_RegisterCCallable("RNifti", "nii2_make_new_nim", (DL_FUNC) &nifti2_make_new_nim);
    R_RegisterCCallable("RNifti", "nii2_simple_init_nim", (DL_FUNC) &nifti2_simple_init_nim);
    R_RegisterCCallable("RNifti", "nii_convert_n1hdr2nim", (DL_FUNC) &nifti_convert_n1hdr2nim);
    R_RegisterCCallable("RNifti", "nii_convert_n2hdr2nim", (DL_FUNC) &nifti_convert_n2hdr2nim);
    R_RegisterCCallable("RNifti", "nii_looks_like_cifti", (DL_FUNC) &nifti_looks_like_cifti);
    R_RegisterCCallable("RNifti", "nii_hdr1_looks_good", (DL_FUNC) &nifti_hdr1_looks_good);
    R_RegisterCCallable("RNifti", "nii_hdr2_looks_good", (DL_FUNC) &nifti_hdr2_looks_good);
    R_RegisterCCallable("RNifti", "nii2_nim_is_valid", (DL_FUNC) &nifti2_nim_is_valid);
    R_RegisterCCallable("RNifti", "nii2_nim_has_valid_dims", (DL_FUNC) &nifti2_nim_has_valid_dims);
    R_RegisterCCallable("RNifti", "is_valid_nii2_type", (DL_FUNC) &is_valid_nifti2_type);
    R_RegisterCCallable("RNifti", "nii2_type_and_names_match", (DL_FUNC) &nifti2_type_and_names_match);
    R_RegisterCCallable("RNifti", "nii2_update_dims_from_array", (DL_FUNC) &nifti2_update_dims_from_array);
    R_RegisterCCallable("RNifti", "nii2_set_iname_offset", (DL_FUNC) &nifti2_set_iname_offset);
    R_RegisterCCallable("RNifti", "nii2_set_type_from_names", (DL_FUNC) &nifti2_set_type_from_names);
    R_RegisterCCallable("RNifti", "nii2_add_extension", (DL_FUNC) &nifti2_add_extension);
    R_RegisterCCallable("RNifti", "nii2_copy_extensions", (DL_FUNC) &nifti2_copy_extensions);
    R_RegisterCCallable("RNifti", "nii2_free_extensions", (DL_FUNC) &nifti2_free_extensions);
    R_RegisterCCallable("RNifti", "nii_get_int64list", (DL_FUNC) &nifti_get_int64list);
    R_RegisterCCallable("RNifti", "valid_nii2_extensions", (DL_FUNC) &valid_nifti2_extensions);
    R_RegisterCCallable("RNifti", "nii_valid_header_size", (DL_FUNC) &nifti_valid_header_size);
}

}
