#ifndef _NIFTI_IMAGE_INTERNAL_H_
#define _NIFTI_IMAGE_INTERNAL_H_

namespace internal {

template <typename SourceType, typename TargetType>
inline TargetType convertValue (SourceType value)
{
    return static_cast<TargetType>(value);
}

// Base class type responsible for handling image data buffers
struct DataHandler
{
    virtual short code () { return DT_NONE; }
    
    template <typename TargetType>
    void convertToArray (void *source, const size_t length, TargetType *target) {}
    
    template <typename TargetType>
    void convertToVector (void *source, const size_t length, std::vector<TargetType> &target, const size_t offset = 0) {}
    
    template <typename SourceType>
    void convertFromVector (const std::vector<SourceType> &source, void *target) {}
    
    template <int SexpType>
    void convertToRcppVector (void *source, const size_t length, Rcpp::Vector<SexpType> &target) {}
};

template <typename Type>
struct TypedDataHandler : public DataHandler
{
    template <typename TargetType>
    void convertToArray (void *source, const size_t length, TargetType *target)
    {
        Type *castSource = static_cast<Type *>(source);
        std::transform(castSource, castSource + length, target, convertValue<Type,TargetType>);
    }
    
    template <typename TargetType>
    void convertToVector (void *source, const size_t length, std::vector<TargetType> &target, const size_t offset = 0)
    {
        Type *castSource = static_cast<Type *>(source) + offset;
        std::transform(castSource, castSource + length, target.begin(), convertValue<Type,TargetType>);
    }
    
    template <typename SourceType>
    void convertFromVector (const std::vector<SourceType> &source, void *target)
    {
        Type *castTarget = static_cast<Type *>(target);
        std::transform(source.begin(), source.end(), castTarget, convertValue<SourceType,Type>);
    }
    
    template <int SexpType>
    void convertToRcppVector (void *source, const size_t length, Rcpp::Vector<SexpType> &target)
    {
        Type *castSource = static_cast<Type *>(source);
        if (SexpType == INTSXP || SexpType == LGLSXP)
            std::transform(castSource, castSource + length, target.begin(), convertValue<Type,int>);
        else if (SexpType == REALSXP)
            std::transform(castSource, castSource + length, target.begin(), convertValue<Type,double>);
        else
            throw std::runtime_error("Only numeric arrays can be created");
    }
};

template <> struct TypedDataHandler<uint8_t> : public DataHandler  { short code () { return DT_UINT8; } };
template <> struct TypedDataHandler<int16_t> : public DataHandler  { short code () { return DT_INT16; } };
template <> struct TypedDataHandler<int32_t> : public DataHandler  { short code () { return DT_INT32; } };
template <> struct TypedDataHandler<float> : public DataHandler    { short code () { return DT_FLOAT32; } };
template <> struct TypedDataHandler<double> : public DataHandler   { short code () { return DT_FLOAT64; } };
template <> struct TypedDataHandler<int8_t> : public DataHandler   { short code () { return DT_INT8; } };
template <> struct TypedDataHandler<uint16_t> : public DataHandler { short code () { return DT_UINT16; } };
template <> struct TypedDataHandler<uint32_t> : public DataHandler { short code () { return DT_UINT32; } };
template <> struct TypedDataHandler<int64_t> : public DataHandler  { short code () { return DT_INT64; } };
template <> struct TypedDataHandler<uint64_t> : public DataHandler { short code () { return DT_UINT64; } };

inline DataHandler * getDataHandler (const short typeCode)
{
    typedef std::auto_ptr<DataHandler> pointer_type;
    static std::map<short,pointer_type> typeMap;
    if (typeMap.empty())
    {
        typeMap[DT_UINT8] = pointer_type(new TypedDataHandler<uint8_t>);
        typeMap[DT_INT16] = pointer_type(new TypedDataHandler<int16_t>);
        typeMap[DT_INT32] = pointer_type(new TypedDataHandler<uint32_t>);
        typeMap[DT_FLOAT32] = pointer_type(new TypedDataHandler<float>);
        typeMap[DT_FLOAT64] = pointer_type(new TypedDataHandler<double>);
        typeMap[DT_INT8] = pointer_type(new TypedDataHandler<int8_t>);
        typeMap[DT_UINT16] = pointer_type(new TypedDataHandler<uint16_t>);
        typeMap[DT_UINT32] = pointer_type(new TypedDataHandler<uint32_t>);
        typeMap[DT_INT64] = pointer_type(new TypedDataHandler<int64_t>);
        typeMap[DT_UINT64] = pointer_type(new TypedDataHandler<uint64_t>);
    }
    
    if (typeMap.count(typeCode) == 0)
        throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(typeCode)) + ")");
    else
        return typeMap[typeCode].get();
};

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

template <typename TargetType>
inline void changeDatatype (nifti_image *image, const short datatype)
{
    TargetType *data;
    const size_t dataSize = image->nvox * sizeof(TargetType);
    data = static_cast<TargetType *>(calloc(1, dataSize));
    
    DataHandler *handler = getDataHandler(image->datatype);
    handler->convertToArray(image->data, image->nvox, data);

    free(image->data);
    image->data = data;
    image->datatype = datatype;
    nifti_datatype_sizes(datatype, &image->nbyper, &image->swapsize);
}

template <int SexpType>
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
        Rcpp::Vector<SexpType> array(static_cast<int>(source->nvox));
        DataHandler *handler = getDataHandler(source->datatype);
        handler->convertToRcppVector(source->data, source->nvox, array);
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
