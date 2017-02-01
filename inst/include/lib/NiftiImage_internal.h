#ifndef _NIFTI_IMAGE_INTERNAL_H_
#define _NIFTI_IMAGE_INTERNAL_H_

namespace internal {

template <typename TargetType>
inline void copyIfPresent (const Rcpp::List &list, const std::set<std::string> names, const std::string &name, TargetType &target)
{
    if (names.count(name) == 1)
        target = Rcpp::as<TargetType>(list[name]);
}

// Special case for char, because Rcpp tries to be too clever and convert it to a string
template <>
inline void copyIfPresent (const Rcpp::List &list, const std::set<std::string> names, const std::string &name, char &target)
{
    if (names.count(name) == 1)
        target = static_cast<char>(Rcpp::as<int>(list[name]));
}

inline mat33 topLeftCorner (const mat44 &matrix)
{
    mat33 newMatrix;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
            newMatrix.m[i][j] = matrix.m[i][j];
    }
    
    return newMatrix;
}

template <typename SourceType, typename TargetType>
inline TargetType convertValue (SourceType value)
{
    return static_cast<TargetType>(value);
}

template <typename SourceType, typename TargetType>
inline void convertArray (const SourceType *source, const size_t length, TargetType *target)
{
    std::transform(source, source + length, target, convertValue<SourceType,TargetType>);
}

template <typename SourceType, typename TargetType>
inline void convertVector (const SourceType *source, const size_t length, std::vector<TargetType> &target)
{
    std::transform(source, source + length, target.begin(), convertValue<SourceType,TargetType>);
}

template <typename TargetType>
inline void changeDatatype (nifti_image *image, const short datatype)
{
    TargetType *data;
    const size_t dataSize = image->nvox * sizeof(TargetType);
    data = static_cast<TargetType *>(calloc(1, dataSize));

    switch (image->datatype)
    {
        case DT_UINT8:
        convertArray(static_cast<uint8_t *>(image->data), image->nvox, data);
        break;
        
        case DT_INT16:
        convertArray(static_cast<int16_t *>(image->data), image->nvox, data);
        break;
        
        case DT_INT32:
        convertArray(static_cast<int32_t *>(image->data), image->nvox, data);
        break;
        
        case DT_FLOAT32:
        convertArray(static_cast<float *>(image->data), image->nvox, data);
        break;
        
        case DT_FLOAT64:
        convertArray(static_cast<double *>(image->data), image->nvox, data);
        break;
        
        case DT_INT8:
        convertArray(static_cast<int8_t *>(image->data), image->nvox, data);
        break;
        
        case DT_UINT16:
        convertArray(static_cast<uint16_t *>(image->data), image->nvox, data);
        break;
        
        case DT_UINT32:
        convertArray(static_cast<uint32_t *>(image->data), image->nvox, data);
        break;
        
        case DT_INT64:
        convertArray(static_cast<int64_t *>(image->data), image->nvox, data);
        break;
        
        case DT_UINT64:
        convertArray(static_cast<uint64_t *>(image->data), image->nvox, data);
        break;
        
        default:
        throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(image->datatype)) + ")");
    }
    
    free(image->data);
    image->data = data;
    image->datatype = datatype;
    nifti_datatype_sizes(datatype, &image->nbyper, &image->swapsize);
}

template <typename SourceType, int SexpType>
inline Rcpp::RObject imageDataToArray (const nifti_image *source)
{
    if (source == NULL)
        return Rcpp::RObject();
    else if (source->data == NULL)
    {
        Rf_warning("Internal image contains no data - filling array with NAs");
        Rcpp::Vector<SexpType> array(static_cast<int>(source->nvox));
        // Rcpp's proxy infrastructure should handle converting NA_REAL to the appropriate NA
        std::fill(array.begin(), array.end(), NA_REAL);
        return array;
    }
    else
    {
        SourceType *original = static_cast<SourceType *>(source->data);
        Rcpp::Vector<SexpType> array(static_cast<int>(source->nvox));
        
        if (SexpType == INTSXP || SexpType == LGLSXP)
            std::transform(original, original + source->nvox, array.begin(), convertValue<SourceType,int>);
        else if (SexpType == REALSXP)
            std::transform(original, original + source->nvox, array.begin(), convertValue<SourceType,double>);
        else
            throw std::runtime_error("Only numeric arrays can be created");
        
        return array;
    }
}

inline void finaliseNiftiImage (SEXP xptr)
{
    NiftiImage *object = (NiftiImage *) R_ExternalPtrAddr(xptr);
    object->setPersistence(false);
    delete object;
    R_ClearExternalPtr(xptr);
}

inline void addAttributes (Rcpp::RObject &object, nifti_image *source, const bool realDim = true)
{
    const int nDims = source->dim[0];
    Rcpp::IntegerVector dim(source->dim+1, source->dim+1+nDims);

    if (realDim)
        object.attr("dim") = dim;
    else
        object.attr("imagedim") = dim;
    
    Rcpp::DoubleVector pixdim(source->pixdim+1, source->pixdim+1+nDims);
    object.attr("pixdim") = pixdim;
    
    if (source->xyz_units == NIFTI_UNITS_UNKNOWN && source->time_units == NIFTI_UNITS_UNKNOWN)
        object.attr("pixunits") = "Unknown";
    else
    {
        Rcpp::CharacterVector pixunits(2);
        pixunits[0] = nifti_units_string(source->xyz_units);
        pixunits[1] = nifti_units_string(source->time_units);
        object.attr("pixunits") = pixunits;
    }
    
    NiftiImage *wrappedSource = new NiftiImage(source, true);
    wrappedSource->setPersistence(true);
    Rcpp::XPtr<NiftiImage> xptr(wrappedSource);
    R_RegisterCFinalizerEx(SEXP(xptr), &finaliseNiftiImage, FALSE);
    object.attr(".nifti_image_ptr") = xptr;
}

} // namespace

#endif
