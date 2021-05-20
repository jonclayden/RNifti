#include <Rcpp.h>

#define RNIFTI_NIFTILIB_VERSION 2
#include "RNifti.h"

using namespace Rcpp;
using namespace RNifti;

typedef std::vector<int> int_vector;
typedef std::vector<NiftiImage::dim_t> dim_vector;
typedef std::vector<NiftiImage::pixdim_t> pixdim_vector;
typedef std::list<NiftiImage::Extension> extension_list;

inline static bool isXformMatrix (const SEXP object)
{
    if (!Rf_isMatrix(object))
        return false;
    NumericMatrix matrix(object);
    return (matrix.cols() == 4 && matrix.rows() == 4);
}

inline static unsigned char clip (const double &value)
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
        dim_vector volumes;
        IntegerVector volumesR(_volumes);
        for (int i=0; i<volumesR.length(); i++)
            volumes.push_back(volumesR[i] - 1);
        const NiftiImage image(as<std::string>(_object), volumes);
        return image.toArrayOrPointer(as<bool>(_internal), "NIfTI image");
    }
END_RCPP
}

RcppExport SEXP readNiftiBlob (SEXP _file, SEXP _length, SEXP _datatype, SEXP _offset, SEXP _swap)
{
BEGIN_RCPP
    int datatype;
    if (Rf_isInteger(_datatype))
        datatype = as<int>(_datatype);
    else
        datatype = RNifti::internal::stringToDatatype(as<std::string>(_datatype));
    
    const std::string filename = as<std::string>(_file);
    const size_t length = as<size_t>(_length);
    const size_t offset = Rf_isNull(_offset) ? 0 : as<size_t>(_offset);
    
    int nbyper;
    nifti_datatype_sizes(datatype, &nbyper, NULL);
    
    const bool gzExtension = filename.length() > 3 && filename.substr(filename.length()-3,3) == ".gz";
    znzFile file = znzopen(filename.c_str(), "rb", gzExtension);
    if (znz_isnull(file))
        Rf_error("Failed to open file %s", filename.c_str());
    if (offset > 0)
        znzseek(file, offset, SEEK_SET);
    char *buffer = (char *) calloc(length, nbyper);
    znzread(buffer, nbyper, length, file);
    znzclose(file);
    
    if (as<bool>(_swap))
        nifti_swap_Nbytes(length, nbyper, buffer);
    
    NiftiImageData data(buffer, length, datatype);
    RObject result;
    if (data.isComplex())
        result = ComplexVector(data.begin(), data.end());
    else if (data.isFloatingPoint())
        result = NumericVector(data.begin(), data.end());
    else
        result = IntegerVector(data.begin(), data.end());
    free(buffer);
    return result;
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
    else if (filetypeString == "nifti2")
        filetype = NIFTI_FTYPE_NIFTI2_1;
    const std::pair<std::string,std::string> paths = image.toFile(as<std::string>(_file), as<std::string>(_datatype), filetype);
    return CharacterVector::create(Named("header")=paths.first, Named("image")=paths.second);
END_RCPP
}

template <typename Header>
static List niftiHeaderToList (const Header &header)
{
    List result;
    
    result["sizeof_hdr"] = header.sizeof_hdr;
    
    result["dim_info"] = int(header.dim_info);
    result["dim"] = std::vector<int>(header.dim, header.dim+8);
    
    result["intent_p1"] = header.intent_p1;
    result["intent_p2"] = header.intent_p2;
    result["intent_p3"] = header.intent_p3;
    result["intent_code"] = header.intent_code;
    
    result["datatype"] = header.datatype;
    result["bitpix"] = header.bitpix;
    
    result["slice_start"] = header.slice_start;
    result["pixdim"] = std::vector<double>(header.pixdim, header.pixdim+8);
    result["vox_offset"] = header.vox_offset;
    result["scl_slope"] = header.scl_slope;
    result["scl_inter"] = header.scl_inter;
    result["slice_end"] = header.slice_end;
    result["slice_code"] = int(header.slice_code);
    result["xyzt_units"] = int(header.xyzt_units);
    result["cal_max"] = header.cal_max;
    result["cal_min"] = header.cal_min;
    result["slice_duration"] = header.slice_duration;
    result["toffset"] = header.toffset;
    result["descrip"] = std::string(header.descrip, 80);
    result["aux_file"] = std::string(header.aux_file, 24);
    
    result["qform_code"] = header.qform_code;
    result["sform_code"] = header.sform_code;
    result["quatern_b"] = header.quatern_b;
    result["quatern_c"] = header.quatern_c;
    result["quatern_d"] = header.quatern_d;
    result["qoffset_x"] = header.qoffset_x;
    result["qoffset_y"] = header.qoffset_y;
    result["qoffset_z"] = header.qoffset_z;
    result["srow_x"] = std::vector<float>(header.srow_x, header.srow_x+4);
    result["srow_y"] = std::vector<float>(header.srow_y, header.srow_y+4);
    result["srow_z"] = std::vector<float>(header.srow_z, header.srow_z+4);
    
    result["intent_name"] = std::string(header.intent_name, 16);
    result["magic"] = std::string(header.magic, 4);
    
    List strings;
    strings["datatype"] = nifti_datatype_string(header.datatype);
    strings["intent_code"] = nifti_intent_string(header.intent_code);
    strings["qform_code"] = nifti_xform_string(header.qform_code);
    strings["sform_code"] = nifti_xform_string(header.sform_code);
    strings["slice_code"] = nifti_slice_string(header.slice_code);
    
    result.attr("class") = CharacterVector::create("niftiHeader");
    result.attr("strings") = strings;
    
    return result;
}

RcppExport SEXP niftiHeader (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image, false, true);
    if (image.isNull())
        return R_NilValue;
    
    const int version = (image->nifti_type == NIFTI_FTYPE_NIFTI2_1 || image->nifti_type == NIFTI_FTYPE_NIFTI2_2) ? 2 : 1;
    List result;
    
    if (version == 1)
    {
        nifti_1_header header;
        nifti_convert_nim2n1hdr(image, &header);
        result = niftiHeaderToList(header);
    }
    else if (version == 2)
    {
        nifti_2_header header;
        nifti_convert_nim2n2hdr(image, &header);
        result = niftiHeaderToList(header);
    }
    
    RNifti::internal::addAttributes(result, image, false, false);
    result.attr("version") = version;
    
    return result;
END_RCPP
}

RcppExport SEXP analyzeHeader (SEXP _image)
{
BEGIN_RCPP
    RObject object(_image);
    nifti_1_header header;
    
    // This special-case treatment of strings is important, because converting
    // ANALYZE files into NIfTI objects and back is somewhat destructive. It
    // mustn't pick up internal images by accident, though
    if (Rf_isString(object) && !object.hasAttribute(".nifti_image_ptr"))
    {
        const std::string path = as<std::string>(object);
        int version;
        void *ptr = nifti2_read_header(RNifti::internal::stringToPath(path), &version, true);
        if (ptr == NULL)
            return R_NilValue;
        else if (version < 0 || version > 1)
            Rf_error("File is not in ANALYZE-7.5 or NIfTI-1 format");
        header = *((nifti_1_header *) ptr);
        free(ptr);
    }
    else
    {
        const NiftiImage image(_image, false, true);
        if (image.isNull())
            return R_NilValue;
        nifti_convert_nim2n1hdr(image, &header);
    }
    
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
        
        double qbcd[3], qxyz[3], dxyz[3], qfac;
        nifti_dmat44_to_quatern(xform, &qbcd[0], &qbcd[1], &qbcd[2], &qxyz[0], &qxyz[1], &qxyz[2], &dxyz[0], &dxyz[1], &dxyz[2], &qfac);
        
        if (as<bool>(_isQform))
        {
            *REAL(image["quatern_b"]) = qbcd[0];
            *REAL(image["quatern_c"]) = qbcd[1];
            *REAL(image["quatern_d"]) = qbcd[2];
            *REAL(image["qoffset_x"]) = qxyz[0];
            *REAL(image["qoffset_y"]) = qxyz[1];
            *REAL(image["qoffset_z"]) = qxyz[2];
            REAL(image["pixdim"])[0] = qfac;
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
            REAL(image["pixdim"])[i+1] = dxyz[i];
        
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
        nifti2_image *ptr = nifti2_make_new_nim(NULL, DT_UNSIGNED_CHAR, 0);
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
        imageString << (const nifti2_image *) image;
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
        const dim_vector dim = image.dim();
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
    // Not flagging as read-only, because a copy must happen or the objects will be tied together
    NiftiImage image(_image);
    return image.toArray();
END_RCPP
}

// Extract an external pointer to a nifti_image, for use in plain C client code
// This involves a copy, and currently only produces version 2 structs
RcppExport SEXP unwrapPointer (SEXP _image, SEXP _disown)
{
BEGIN_RCPP
    const NiftiImage image(_image, true, true);
    nifti2_image *result = nifti2_copy_nim_info(image);
    if (image->data != NULL)
    {
        const size_t dataSize = nifti2_get_volsize(image);
        result->data = calloc(1, dataSize);
        memcpy(result->data, image->data, dataSize);
    }
    
    // Create the bare pointer, and set nifti2_image_free() as a finaliser function if we're not disowning
    XPtr<nifti2_image,PreserveStorage,nifti2_image_free,true> pointer(result, !as<bool>(_disown));
    return pointer;
END_RCPP
}

RcppExport SEXP wrapPointer (SEXP _image)
{
BEGIN_RCPP
    // Finalisers are not set when the pointer is retrieved, so the other template parameters aren't needed
    XPtr<nifti2_image> pointer(_image);
    const NiftiImage image(pointer.get(), true);
    return image.toPointer("NIfTI image");
END_RCPP
}

RcppExport SEXP getExtensions (SEXP _image, SEXP _code)
{
BEGIN_RCPP
    const NiftiImage image(_image, false, true);
    const extension_list extensions = image.extensions(as<int>(_code));
    List result(extensions.begin(), extensions.end());
    return result;
END_RCPP
}

RcppExport SEXP setExtensions (SEXP _image, SEXP _extensions, SEXP _code)
{
BEGIN_RCPP
    NiftiImage image(_image);
    const int code = as<int>(_code);
    extension_list extensions = image.extensions();
    
    if (Rf_isNull(_extensions))
    {
        if (code < 0)
            extensions.clear();
        else
        {
            for (extension_list::iterator it=extensions.begin(); it!=extensions.end(); ++it)
            {
                if (it->code() == code)
                    extensions.erase(it);
            }
        }
    }
    else if (Rf_isVectorList(_extensions))
        extensions = as<extension_list>(_extensions);
    else
        extensions.push_back(NiftiImage::Extension(_extensions, code));
    
    return image.replaceExtensions(extensions).toArrayOrPointer(Rf_inherits(_image,"internalImage"), "NIfTI image");
END_RCPP
}

RcppExport SEXP setDebugLevel (SEXP _level)
{
BEGIN_RCPP
    nifti_set_debug_level(as<int>(_level));
    return R_NilValue;
END_RCPP
}

extern "C" {

R_CallMethodDef callMethods[] = {
    { "packRgb",        (DL_FUNC) &packRgb,         3 },
    { "rgbToStrings",   (DL_FUNC) &rgbToStrings,    1 },
    { "unpackRgb",      (DL_FUNC) &unpackRgb,       2 },
    { "asNifti",        (DL_FUNC) &asNifti,         4 },
    { "niftiVersion",   (DL_FUNC) &niftiVersion,    1 },
    { "readNifti",      (DL_FUNC) &readNifti,       3 },
    { "readNiftiBlob",  (DL_FUNC) &readNiftiBlob,   5 },
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
    { "unwrapPointer",  (DL_FUNC) &unwrapPointer,   2 },
    { "wrapPointer",    (DL_FUNC) &wrapPointer,     1 },
    { "getExtensions",  (DL_FUNC) &getExtensions,   2 },
    { "setExtensions",  (DL_FUNC) &setExtensions,   3 },
    { "setDebugLevel",  (DL_FUNC) &setDebugLevel,   1 },
    { NULL, NULL, 0 }
};

}
