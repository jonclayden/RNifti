#ifndef _NIFTI_IMAGE_INTERNAL_H_
#define _NIFTI_IMAGE_INTERNAL_H_

namespace internal {

template <typename ElementType>
struct DataRescaler
{
    float slope, intercept;
    
    DataRescaler (const float slope, const float intercept)
        : slope(slope), intercept(intercept) {}
    
    inline ElementType operator() (const ElementType value)
    {
        return static_cast<ElementType>(value * slope + intercept);
    }
};

template <typename SourceType, typename TargetType>
inline TargetType convertValue (SourceType value)
{
    if (std::numeric_limits<SourceType>::has_quiet_NaN && !std::numeric_limits<TargetType>::has_quiet_NaN && ISNAN(value))
        return TargetType(0);
    else
        return static_cast<TargetType>(value);
}

template <typename TargetType, class OutputIterator>
inline void convertData (void *source, const short datatype, const size_t length, OutputIterator target, const ptrdiff_t offset = 0)
{
    if (source == NULL)
        return;
    
    switch (datatype)
    {
        case DT_UINT8:
        {
            uint8_t *castSource = static_cast<uint8_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<uint8_t,TargetType>);
            break;
        }
        
        case DT_INT16:
        {
            int16_t *castSource = static_cast<int16_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<int16_t,TargetType>);
            break;
        }
        
        case DT_INT32:
        {
            int32_t *castSource = static_cast<int32_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<int32_t,TargetType>);
            break;
        }
        
        case DT_FLOAT32:
        {
            float *castSource = static_cast<float *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<float,TargetType>);
            break;
        }
        
        case DT_FLOAT64:
        {
            double *castSource = static_cast<double *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<double,TargetType>);
            break;
        }
        
        case DT_INT8:
        {
            int8_t *castSource = static_cast<int8_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<int8_t,TargetType>);
            break;
        }
        
        case DT_UINT16:
        {
            uint16_t *castSource = static_cast<uint16_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<uint16_t,TargetType>);
            break;
        }
        
        case DT_UINT32:
        {
            uint32_t *castSource = static_cast<uint32_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<uint32_t,TargetType>);
            break;
        }
        
        case DT_INT64:
        {
            int64_t *castSource = static_cast<int64_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<int64_t,TargetType>);
            break;
        }
        
        case DT_UINT64:
        {
            uint64_t *castSource = static_cast<uint64_t *>(source) + offset;
            std::transform(castSource, castSource + length, target, convertValue<uint64_t,TargetType>);
            break;
        }
        
        default:
        throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(datatype)) + ")");
    }
}

template <typename SourceType, class InputIterator>
inline void replaceData (InputIterator begin, InputIterator end, void *target, const short datatype)
{
    if (target == NULL)
        return;
    
    switch (datatype)
    {
        case DT_UINT8:
        std::transform(begin, end, static_cast<uint8_t *>(target), convertValue<SourceType,uint8_t>);
        break;
        
        case DT_INT16:
        std::transform(begin, end, static_cast<int16_t *>(target), convertValue<SourceType,int16_t>);
        break;
        
        case DT_INT32:
        std::transform(begin, end, static_cast<int32_t *>(target), convertValue<SourceType,int32_t>);
        break;
        
        case DT_FLOAT32:
        std::transform(begin, end, static_cast<float *>(target), convertValue<SourceType,float>);
        break;
        
        case DT_FLOAT64:
        std::transform(begin, end, static_cast<double *>(target), convertValue<SourceType,double>);
        break;
        
        case DT_INT8:
        std::transform(begin, end, static_cast<int8_t *>(target), convertValue<SourceType,int8_t>);
        break;
        
        case DT_UINT16:
        std::transform(begin, end, static_cast<uint16_t *>(target), convertValue<SourceType,uint16_t>);
        break;
        
        case DT_UINT32:
        std::transform(begin, end, static_cast<uint32_t *>(target), convertValue<SourceType,uint32_t>);
        break;
        
        case DT_INT64:
        std::transform(begin, end, static_cast<int64_t *>(target), convertValue<SourceType,int64_t>);
        break;
        
        case DT_UINT64:
        std::transform(begin, end, static_cast<uint64_t *>(target), convertValue<SourceType,uint64_t>);
        break;
        
        default:
        throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(datatype)) + ")");
    }
}

inline short stringToDatatype (const std::string &datatype)
{
    static std::map<std::string,short> datatypeCodes;
    if (datatypeCodes.empty())
    {
        datatypeCodes["auto"] = DT_NONE;
        datatypeCodes["none"] = DT_NONE;
        datatypeCodes["unknown"] = DT_NONE;
        datatypeCodes["uint8"] = DT_UINT8;
        datatypeCodes["char"] = DT_UINT8;
        datatypeCodes["int16"] = DT_INT16;
        datatypeCodes["short"] = DT_INT16;
        datatypeCodes["int32"] = DT_INT32;
        datatypeCodes["int"] = DT_INT32;
        datatypeCodes["float32"] = DT_FLOAT32;
        datatypeCodes["float"] = DT_FLOAT32;
        datatypeCodes["float64"] = DT_FLOAT64;
        datatypeCodes["double"] = DT_FLOAT64;
        datatypeCodes["int8"] = DT_INT8;
        datatypeCodes["uint16"] = DT_UINT16;
        datatypeCodes["uint32"] = DT_UINT32;
        datatypeCodes["int64"] = DT_INT64;
        datatypeCodes["uint64"] = DT_UINT64;
    }
    
    std::locale locale;
    std::string lowerCaseDatatype = datatype;
    for (std::string::size_type i=0; i<lowerCaseDatatype.length(); i++)
        lowerCaseDatatype[i] = std::tolower(lowerCaseDatatype[i], locale);
    
    if (datatypeCodes.count(lowerCaseDatatype) == 0)
    {
        std::ostringstream message;
        message << "Datatype \"" << datatype << "\" is not valid";
        Rf_warning(message.str().c_str());
        return DT_NONE;
    }
    else
        return datatypeCodes[lowerCaseDatatype];
}

#ifndef _NO_R__

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

#endif // _NO_R__

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

#ifndef _NO_R__

template <int SexpType>
inline Rcpp::RObject imageDataToArray (const nifti_image *source)
{
    if (source == NULL)
        return Rcpp::RObject();
    else if (source->data == NULL)
    {
        Rf_warning("Internal image contains no data - filling array with NAs");
        Rcpp::Vector<SexpType> array(static_cast<int>(source->nvox));
        if (SexpType == INTSXP || SexpType == LGLSXP)
            std::fill(array.begin(), array.end(), NA_INTEGER);
        else if (SexpType == REALSXP)
            std::fill(array.begin(), array.end(), NA_REAL);
        else
            throw std::runtime_error("Only numeric arrays can be created");
        return array;
    }
    else
    {
        Rcpp::Vector<SexpType> array(static_cast<int>(source->nvox));
        if (SexpType == INTSXP || SexpType == LGLSXP)
            convertData<int>(source->data, source->datatype, source->nvox, array.begin());
        else if (SexpType == REALSXP)
            convertData<double>(source->data, source->datatype, source->nvox, array.begin());
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

#endif // _NO_R__

} // namespace

#endif
