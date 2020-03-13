#ifndef _NIFTI_IMAGE_IMPL_H_
#define _NIFTI_IMAGE_IMPL_H_

namespace internal {

struct vec3
{
    float v[3];

    vec3 operator-() const
    {
        vec3 r;
        r.v[0] = -v[0];
        r.v[1] = -v[1];
        r.v[2] = -v[2];
        return r;
    }
};

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

inline vec3 matrixVectorProduct (const mat33 &matrix, const vec3 &vector)
{
    vec3 newVector;
    for (int i=0; i<3; i++)
    {
        newVector.v[i] = 0.0;
        for (int j=0; j<3; j++)
            newVector.v[i] += matrix.m[i][j] * vector.v[j];
    }
    return newVector;
}

// A poor man's NaN check, but should work whenever proper IEEE arithmetic is being used
template <typename Type>
inline bool isNaN (const Type x) { return (x != x); }

#ifdef USING_R
// R offers the portable ISNAN macro for doubles, which is more robust
template <>
inline bool isNaN<double> (const double x) { return bool(ISNAN(x)); }

// For R specifically, we have to catch NA_INTEGER (a.k.a. INT_MIN)
template <>
inline bool isNaN<int> (const int x) { return (x == NA_INTEGER); }

template <>
inline bool isNaN<rgba32_t> (const rgba32_t x) { return (x.value.packed == NA_INTEGER); }
#endif

template <typename Type>
inline bool lessThan (Type a, Type b) { return (!isNaN(a) && !isNaN(b) && a < b); }

inline double roundEven (const double value)
{
    if (isNaN(value))
        return value;

    double whole;
    double frac = std::fabs(std::modf(value, &whole));
    double sign = (value < 0.0 ? -1.0 : 1.0);

    if (frac < 0.5)
        return whole;
    else if (frac > 0.5)
        return whole + sign;
    else if (std::fmod(whole, 2.0) < 0.0001)
        return whole;
    else
        return whole + sign;
}

inline int stringToDatatype (const std::string &datatype)
{
    static std::map<std::string,int> datatypeCodes;
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
        datatypeCodes["complex64"] = DT_COMPLEX64;
        datatypeCodes["complex128"] = DT_COMPLEX128;
        datatypeCodes["complex"] = DT_COMPLEX128;
        datatypeCodes["rgb24"] = DT_RGB24;
        datatypeCodes["rgb"] = DT_RGB24;
        datatypeCodes["rgba32"] = DT_RGBA32;
        datatypeCodes["rgba"] = DT_RGBA32;
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

#ifdef USING_R
inline const char * stringToPath (const std::string &str) { return R_ExpandFileName(str.c_str()); }
#else
inline const char * stringToPath (const std::string &str) { return str.c_str(); }
#endif

#ifdef USING_R

template <typename TargetType>
inline void copyIfPresent (const Rcpp::List &list, const std::set<std::string> names, const std::string &name, TargetType &target)
{
    if (names.count(name) == 1)
    {
        const Rcpp::RObject object = list[name];
        const int length = Rf_length(object);
        if (length == 0)
        {
            std::ostringstream message;
            message << "Field \"" << name << "\" is empty and will be ignored";
            Rf_warning(message.str().c_str());
        }
        else if (length > 1)
        {
            std::ostringstream message;
            message << "Field \"" << name << "\" has " << length << "elements, but only the first will be used";
            Rf_warning(message.str().c_str());
            target = Rcpp::as< std::vector<TargetType> >(object)[0];
        }
        else
            target = Rcpp::as<TargetType>(object);
    }
}

// Special case for char, because Rcpp tries to be too clever and convert it to a string
template <>
inline void copyIfPresent (const Rcpp::List &list, const std::set<std::string> names, const std::string &name, char &target)
{
    if (names.count(name) == 1)
    {
        int intValue = 0;
        copyIfPresent<int>(list, names, name, intValue);
        target = static_cast<char>(intValue);
    }
}

inline void updateHeader (nifti_1_header *header, const Rcpp::List &list, const bool ignoreDatatype = false)
{
    if (header == NULL || Rf_isNull(list.names()))
        return;
    
    const Rcpp::CharacterVector _names = list.names();
    std::set<std::string> names;
    for (Rcpp::CharacterVector::const_iterator it=_names.begin(); it!=_names.end(); it++)
        names.insert(Rcpp::as<std::string>(*it));
    
    copyIfPresent(list, names, "sizeof_hdr", header->sizeof_hdr);
    
    copyIfPresent(list, names, "dim_info", header->dim_info);
    if (names.count("dim") == 1)
    {
        std::vector<short> dim = list["dim"];
        if (dim.size() != 8)
            throw std::runtime_error("Field \"dim\" must contain 8 elements");
        for (size_t i=0; i<8; i++)
            header->dim[i] = dim[i];
    }
    
    copyIfPresent(list, names, "intent_p1", header->intent_p1);
    copyIfPresent(list, names, "intent_p2", header->intent_p2);
    copyIfPresent(list, names, "intent_p3", header->intent_p3);
    copyIfPresent(list, names, "intent_code", header->intent_code);
    
    if (!ignoreDatatype)
    {
        copyIfPresent(list, names, "datatype", header->datatype);
        copyIfPresent(list, names, "bitpix", header->bitpix);
    }
    
    copyIfPresent(list, names, "slice_start", header->slice_start);
    if (names.count("pixdim") == 1)
    {
        std::vector<float> pixdim = list["pixdim"];
        if (pixdim.size() != 8)
            throw std::runtime_error("Field \"pixdim\" must contain 8 elements");
        for (size_t i=0; i<8; i++)
            header->pixdim[i] = pixdim[i];
    }
    copyIfPresent(list, names, "vox_offset", header->vox_offset);
    copyIfPresent(list, names, "scl_slope", header->scl_slope);
    copyIfPresent(list, names, "scl_inter", header->scl_inter);
    copyIfPresent(list, names, "slice_end", header->slice_end);
    copyIfPresent(list, names, "slice_code", header->slice_code);
    copyIfPresent(list, names, "xyzt_units", header->xyzt_units);
    copyIfPresent(list, names, "cal_max", header->cal_max);
    copyIfPresent(list, names, "cal_min", header->cal_min);
    copyIfPresent(list, names, "slice_duration", header->slice_duration);
    copyIfPresent(list, names, "toffset", header->toffset);
    
    if (names.count("descrip") == 1)
        strcpy(header->descrip, Rcpp::as<std::string>(list["descrip"]).substr(0,79).c_str());
    if (names.count("aux_file") == 1)
        strcpy(header->aux_file, Rcpp::as<std::string>(list["aux_file"]).substr(0,23).c_str());
    
    copyIfPresent(list, names, "qform_code", header->qform_code);
    copyIfPresent(list, names, "sform_code", header->sform_code);
    copyIfPresent(list, names, "quatern_b", header->quatern_b);
    copyIfPresent(list, names, "quatern_c", header->quatern_c);
    copyIfPresent(list, names, "quatern_d", header->quatern_d);
    copyIfPresent(list, names, "qoffset_x", header->qoffset_x);
    copyIfPresent(list, names, "qoffset_y", header->qoffset_y);
    copyIfPresent(list, names, "qoffset_z", header->qoffset_z);
    
    if (names.count("srow_x") == 1)
    {
        std::vector<float> srow_x = list["srow_x"];
        if (srow_x.size() != 4)
            throw std::runtime_error("Field \"srow_x\" must contain 4 elements");
        for (size_t i=0; i<4; i++)
            header->srow_x[i] = srow_x[i];
    }
    if (names.count("srow_y") == 1)
    {
        std::vector<float> srow_y = list["srow_y"];
        if (srow_y.size() != 4)
            throw std::runtime_error("Field \"srow_y\" must contain 4 elements");
        for (size_t i=0; i<4; i++)
            header->srow_y[i] = srow_y[i];
    }
    if (names.count("srow_z") == 1)
    {
        std::vector<float> srow_z = list["srow_z"];
        if (srow_z.size() != 4)
            throw std::runtime_error("Field \"srow_z\" must contain 4 elements");
        for (size_t i=0; i<4; i++)
            header->srow_z[i] = srow_z[i];
    }
    
    if (names.count("intent_name") == 1)
        strcpy(header->intent_name, Rcpp::as<std::string>(list["intent_name"]).substr(0,15).c_str());
    if (names.count("magic") == 1)
        strcpy(header->magic, Rcpp::as<std::string>(list["magic"]).substr(0,3).c_str());
}

inline void addAttributes (const SEXP pointer, const NiftiImage &source, const bool realDim = true, const bool includeXptr = true, const bool keepData = true)
{
    const int nDims = source->dim[0];
    Rcpp::RObject object(pointer);
    Rcpp::IntegerVector dim(source->dim+1, source->dim+1+nDims);

    if (realDim)
        object.attr("dim") = dim;
    else
        object.attr("imagedim") = dim;
    
    Rcpp::DoubleVector pixdim(nDims);
    for (int i=0; i<nDims; i++)
        pixdim[i] = std::abs(static_cast<double>(source->pixdim[i+1]));
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
    
    if (includeXptr)
    {
        NiftiImage *imagePtr = new NiftiImage(source, false);
        if (!keepData)
            imagePtr->dropData();
        Rcpp::XPtr<NiftiImage> xptr(imagePtr);
        object.attr(".nifti_image_ptr") = xptr;
    }
}

#endif  // USING_R

}       // internal namespace

template <typename Type, bool alpha>
inline void NiftiImageData::ConcreteTypeHandler<Type,alpha>::minmax (void *ptr, const size_t length, double *min, double *max) const
{
    if (ptr == NULL || length < 1)
    {
        *min = static_cast<double>(std::numeric_limits<Type>::min());
        *max = static_cast<double>(std::numeric_limits<Type>::max());
    }
    else
    {
        Type *loc = static_cast<Type*>(ptr);
        Type currentMin = *loc, currentMax = *loc;
        for (size_t i=1; i<length; i++)
        {
            loc++;
            if (internal::lessThan(*loc, currentMin))
                currentMin = *loc;
            if (internal::lessThan(currentMax, *loc))
                currentMax = *loc;
        }
        *min = static_cast<double>(currentMin);
        *max = static_cast<double>(currentMax);
    }
}

template <typename ElementType>
inline void NiftiImageData::ConcreteTypeHandler<std::complex<ElementType>,false>::minmax (void *ptr, const size_t length, double *min, double *max) const
{
    if (ptr == NULL || length < 1)
    {
        *min = static_cast<double>(std::numeric_limits<ElementType>::min());
        *max = static_cast<double>(std::numeric_limits<ElementType>::max());
    }
    else
    {
        ElementType *loc = static_cast<ElementType*>(ptr);
        ElementType currentMin = *loc, currentMax = *loc;
        for (size_t i=1; i<(2*length); i++)
        {
            loc++;
            if (internal::lessThan(*loc, currentMin))
                currentMin = *loc;
            if (internal::lessThan(currentMax, *loc))
                currentMax = *loc;
        }
        *min = static_cast<double>(currentMin);
        *max = static_cast<double>(currentMax);
    }
}

template <typename SourceType>
inline NiftiImageData::Element & NiftiImageData::Element::operator= (const SourceType &value)
{
    if (internal::isNaN(value))
    {
        if (!parent.handler->hasNaN())
        {
            const double zeroValue = parent.isScaled() ? (-parent.intercept / parent.slope) : 0.0;
            if (parent.isFloatingPoint())
                parent.handler->setDouble(ptr, zeroValue);
            else
                parent.handler->setInt(ptr, static_cast<int>(internal::roundEven(zeroValue)));
        }
#ifdef USING_R
        // Only happens for integer types that admit an NaN/NA value.
        // In practice this means int specifically for R, so we don't
        // need to worry about the effect of casting INT_MIN to a wider
        // or narrower type
        else if (parent.isInteger())
            parent.handler->setInt(ptr, NA_INTEGER);
#endif
        else
            parent.handler->setDouble(ptr, std::numeric_limits<double>::quiet_NaN());
    }
    else if (parent.isScaled())
    {
        double reverseScaledValue = (static_cast<double>(value) - parent.intercept) / parent.slope;
        if (parent.isFloatingPoint())
            parent.handler->setDouble(ptr, reverseScaledValue);
        else
            parent.handler->setInt(ptr, static_cast<int>(internal::roundEven(reverseScaledValue)));
    }
    else if (std::numeric_limits<SourceType>::is_integer)
        parent.handler->setInt(ptr, static_cast<int>(value));
    else
        parent.handler->setDouble(ptr, static_cast<double>(value));
    return *this;
}

inline NiftiImageData::Element & NiftiImageData::Element::operator= (const NiftiImageData::Element &other)
{
    if (other.parent.isScaled() || other.parent.isFloatingPoint())
    {
        const double value = other;
        *this = value;
    }
    else
    {
        const int value = other;
        *this = value;
    }
    return *this;
}

inline mat33 NiftiImage::xformToRotation (const mat44 matrix)
{
    float qb, qc, qd, qfac;
    nifti_mat44_to_quatern(matrix, &qb, &qc, &qd, NULL, NULL, NULL, NULL, NULL, NULL, &qfac);
    mat44 rotationMatrix = nifti_quatern_to_mat44(qb, qc, qd, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, qfac);
    return internal::topLeftCorner(rotationMatrix);
}

inline std::string NiftiImage::xformToString (const mat44 matrix)
{
    int icode, jcode, kcode;
    nifti_mat44_to_orientation(matrix, &icode, &jcode, &kcode);
    
    int codes[3] = { icode, jcode, kcode };
    std::string result("---");
    for (int i=0; i<3; i++)
    {
        switch (codes[i])
        {
            case NIFTI_L2R: result[i] = 'R'; break;
            case NIFTI_R2L: result[i] = 'L'; break;
            case NIFTI_P2A: result[i] = 'A'; break;
            case NIFTI_A2P: result[i] = 'P'; break;
            case NIFTI_I2S: result[i] = 'S'; break;
            case NIFTI_S2I: result[i] = 'I'; break;
        }
    }
    return result;
}

inline int NiftiImage::fileVersion (const std::string &path)
{
#if RNIFTI_NIFTILIB_VERSION == 1
    nifti_1_header *header = nifti_read_header(internal::stringToPath(path), NULL, false);
    if (header == NULL)
        return -1;
    else
    {
        int version = NIFTI_VERSION(*header);
        if (version == 0)
        {
            // NIfTI-2 has a 540-byte header - check for this or its byte-swapped equivalent
            if (header->sizeof_hdr == 540 || header->sizeof_hdr == 469893120)
            {
                // The magic number has moved in NIfTI-2, so find it by byte offset
                const char *magic = (char *) header + 4;
                if (strncmp(magic,"ni2",3) == 0 || strncmp(magic,"n+2",3) == 0)
                    version = 2;
            }
            else if (!nifti_hdr_looks_good(header))
            {
                // Not plausible as ANALYZE, so return -1
                version = -1;
            }
        }
        free(header);
        return version;
    }
#elif RNIFTI_NIFTILIB_VERSION == 2
    int version;
    void *header = nifti2_read_header(internal::stringToPath(path), &version, true);
    free(header);
    return version;
#endif
}

inline void NiftiImage::acquire (nifti_image * const image)
{
    // If we're taking ownership of a new image, release the old one
    if (this->image != NULL && this->image != image)
        release();
    
    // Set the internal pointer and create or update the reference counter
    this->image = image;
    if (image != NULL)
    {
        if (this->refCount == NULL)
            this->refCount = new int(1);
        else
            (*this->refCount)++;
        
#ifndef NDEBUG
        Rc_printf("Acquiring pointer %p (reference count is %d)\n", this->image, *this->refCount);
#endif
    }
}

inline void NiftiImage::release ()
{
    if (this->image != NULL)
    {
        if (this->refCount != NULL)
        {
            (*this->refCount)--;
#ifndef NDEBUG
            Rc_printf("Releasing pointer %p (reference count is %d)\n", this->image, *this->refCount);
#endif
            if (*this->refCount < 1)
            {
                nifti_image_free(this->image);
                this->image = NULL;
                delete this->refCount;
                this->refCount = NULL;
            }
        }
        else
            Rc_printf("Releasing untracked object %p", this->image);
    }
}

inline void NiftiImage::copy (const nifti_image *source)
{
    if (source == NULL)
        acquire(NULL);
    else
    {
        acquire(nifti_copy_nim_info(source));
        if (source->data != NULL)
        {
            size_t dataSize = nifti_get_volsize(source);
            image->data = calloc(1, dataSize);
            memcpy(image->data, source->data, dataSize);
        }
    }
}

inline void NiftiImage::copy (const NiftiImage &source)
{
    const nifti_image *sourceStruct = source;
    copy(sourceStruct);
}

inline void NiftiImage::copy (const Block &source)
{
    const nifti_image *sourceStruct = source.image;
    if (sourceStruct == NULL)
        acquire(NULL);
    else
    {
        acquire(nifti_copy_nim_info(sourceStruct));
        image->dim[0] = source.image->dim[0] - 1;
        image->dim[source.dimension] = 1;
        image->pixdim[source.dimension] = 1.0;
        nifti_update_dims_from_array(image);
        
        if (sourceStruct->data != NULL)
        {
            size_t blockSize = nifti_get_volsize(image);
            image->data = calloc(1, blockSize);
            memcpy(image->data, static_cast<char*>(source.image->data) + blockSize*source.index, blockSize);
        }
    }
}

#ifdef USING_R

// Convert an S4 "nifti" object, as defined in the oro.nifti package, to a "nifti_image" struct
inline void NiftiImage::initFromNiftiS4 (const Rcpp::RObject &object, const bool copyData)
{
    nifti_1_header header;
    header.sizeof_hdr = 348;
    
    const std::vector<int> dims = object.slot("dim_");
    for (int i=0; i<8; i++)
        header.dim[i] = dims[i];
    
    header.intent_p1 = object.slot("intent_p1");
    header.intent_p2 = object.slot("intent_p2");
    header.intent_p3 = object.slot("intent_p3");
    header.intent_code = object.slot("intent_code");
    
    header.datatype = object.slot("datatype");
    header.bitpix = object.slot("bitpix");
    
    header.slice_start = object.slot("slice_start");
    header.slice_end = object.slot("slice_end");
    header.slice_code = Rcpp::as<int>(object.slot("slice_code"));
    header.slice_duration = object.slot("slice_duration");
    
    const std::vector<float> pixdims = object.slot("pixdim");
    for (int i=0; i<8; i++)
        header.pixdim[i] = pixdims[i];
    header.xyzt_units = Rcpp::as<int>(object.slot("xyzt_units"));
    
    header.vox_offset = object.slot("vox_offset");
    
    // oro.nifti does its own data rescaling, so we ignore the slope and intercept fields
    header.scl_slope = 0.0;
    header.scl_inter = 0.0;
    header.toffset = object.slot("toffset");
    
    header.cal_max = object.slot("cal_max");
    header.cal_min = object.slot("cal_min");
    header.glmax = header.glmin = 0;
    
    strncpy(header.descrip, Rcpp::as<std::string>(object.slot("descrip")).c_str(), 79);
    header.descrip[79] = '\0';
    strncpy(header.aux_file, Rcpp::as<std::string>(object.slot("aux_file")).c_str(), 23);
    header.aux_file[23] = '\0';
    strncpy(header.intent_name, Rcpp::as<std::string>(object.slot("intent_name")).c_str(), 15);
    header.intent_name[15] = '\0';
    strncpy(header.magic, Rcpp::as<std::string>(object.slot("magic")).c_str(), 3);
    header.magic[3] = '\0';
    
    header.qform_code = object.slot("qform_code");
    header.sform_code = object.slot("sform_code");
    
    header.quatern_b = object.slot("quatern_b");
    header.quatern_c = object.slot("quatern_c");
    header.quatern_d = object.slot("quatern_d");
    header.qoffset_x = object.slot("qoffset_x");
    header.qoffset_y = object.slot("qoffset_y");
    header.qoffset_z = object.slot("qoffset_z");
    
    const std::vector<float> srow_x = object.slot("srow_x");
    const std::vector<float> srow_y = object.slot("srow_y");
    const std::vector<float> srow_z = object.slot("srow_z");
    for (int i=0; i<4; i++)
    {
        header.srow_x[i] = srow_x[i];
        header.srow_y[i] = srow_y[i];
        header.srow_z[i] = srow_z[i];
    }
    
    // Ignoring complex and RGB types here because oro.nifti doesn't yet support them
    if (header.datatype == DT_UINT8 || header.datatype == DT_INT16 || header.datatype == DT_INT32 || header.datatype == DT_INT8 || header.datatype == DT_UINT16 || header.datatype == DT_UINT32)
        header.datatype = DT_INT32;
    else if (header.datatype == DT_FLOAT32 || header.datatype == DT_FLOAT64)
        header.datatype = DT_FLOAT64;
    else
        throw std::runtime_error("Data type is not supported");
    
    acquire(nifti_convert_nhdr2nim(header, NULL));
    
    const SEXP data = PROTECT(object.slot(".Data"));
    if (!copyData || Rf_length(data) <= 1)
        this->image->data = NULL;
    else if (header.datatype == DT_INT32)
    {
        Rcpp::IntegerVector intData(data);
        replaceData(NiftiImageData(intData.begin(), intData.end(), DT_INT32));
    }
    else
    {
        Rcpp::DoubleVector doubleData(data);
        replaceData(NiftiImageData(doubleData.begin(), doubleData.end(), DT_FLOAT64));
    }
    UNPROTECT(1);
}

inline void NiftiImage::initFromMriImage (const Rcpp::RObject &object, const bool copyData)
{
    Rcpp::Reference mriImage(object);
    Rcpp::Function getXform = mriImage.field("getXform");
    Rcpp::NumericMatrix xform = getXform();
    
    acquire(NULL);
    
    if (Rf_length(mriImage.field("tags")) > 0)
        initFromList(mriImage.field("tags"));
    
    Rcpp::RObject data = mriImage.field("data");
    if (data.inherits("SparseArray"))
    {
        Rcpp::Language call("as.array", data);
        data = call.eval();
    }
    
    const int datatype = (Rf_isNull(data) ? DT_INT32 : sexpTypeToNiftiType(data.sexp_type()));
    
    int dims[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    const std::vector<int> dimVector = mriImage.field("imageDims");
    const int nDims = std::min(7, int(dimVector.size()));
    dims[0] = nDims;
    size_t nVoxels = 1;
    for (int i=0; i<nDims; i++)
    {
        dims[i+1] = dimVector[i];
        nVoxels *= dimVector[i];
    }
    
    if (this->image == NULL)
        acquire(nifti_make_new_nim(dims, datatype, FALSE));
    else
    {
        std::copy(dims, dims+8, this->image->dim);
        this->image->datatype = datatype;
        nifti_datatype_sizes(image->datatype, &image->nbyper, NULL);
    }
    
    if (copyData && !Rf_isNull(data))
    {
        // NB: nifti_get_volsize() will not be right here if there were tags
        const size_t dataSize = nVoxels * image->nbyper;
        this->image->data = calloc(1, dataSize);
        if (datatype == DT_INT32)
            memcpy(this->image->data, INTEGER(data), dataSize);
        else
            memcpy(this->image->data, REAL(data), dataSize);
    }
    else
        this->image->data = NULL;
    
    const std::vector<float> pixdimVector = mriImage.field("voxelDims");
    const int pixdimLength = pixdimVector.size();
    for (int i=0; i<std::min(pixdimLength,nDims); i++)
        this->image->pixdim[i+1] = std::abs(pixdimVector[i]);
    
    const std::vector<std::string> pixunitsVector = mriImage.field("voxelDimUnits");
    setPixunits(pixunitsVector);
    
    if (xform.rows() != 4 || xform.cols() != 4)
        this->image->qform_code = this->image->sform_code = 0;
    else
    {
        mat44 matrix;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
                matrix.m[i][j] = static_cast<float>(xform(i,j));
        }
        
        this->image->qto_xyz = matrix;
        this->image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
        nifti_mat44_to_quatern(image->qto_xyz, &image->quatern_b, &image->quatern_c, &image->quatern_d, &image->qoffset_x, &image->qoffset_y, &image->qoffset_z, NULL, NULL, NULL, &image->qfac);
        
        this->image->sto_xyz = matrix;
        this->image->sto_ijk = nifti_mat44_inverse(image->sto_xyz);
        
        this->image->qform_code = this->image->sform_code = 2;
    }
}

inline void NiftiImage::initFromList (const Rcpp::RObject &object)
{
    Rcpp::List list(object);
    nifti_1_header *header = nifti_make_new_header(NULL, DT_FLOAT64);
    
    internal::updateHeader(header, list);
    
    acquire(nifti_convert_nhdr2nim(*header, NULL));
    this->image->data = NULL;
    free(header);
}

inline void NiftiImage::initFromArray (const Rcpp::RObject &object, const bool copyData)
{
    int dims[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    const std::vector<int> dimVector = object.attr("dim");
    
    const int nDims = std::min(7, int(dimVector.size()));
    dims[0] = nDims;
    for (int i=0; i<nDims; i++)
        dims[i+1] = dimVector[i];
    
    int datatype = sexpTypeToNiftiType(object.sexp_type());
    if (object.inherits("rgbArray"))
    {
        const int channels = (object.hasAttribute("channels") ? object.attr("channels") : 3);
        datatype = (channels == 4 ? DT_RGBA32 : DT_RGB24);
    }
    acquire(nifti_make_new_nim(dims, datatype, int(copyData)));
    
    if (copyData)
    {
        const size_t dataSize = nifti_get_volsize(image);
        if (datatype == DT_INT32 || datatype == DT_RGBA32)
            memcpy(this->image->data, INTEGER(object), dataSize);
        else if (datatype == DT_RGB24)
        {
            NiftiImageData data(image);
            std::copy(INTEGER(object), INTEGER(object)+image->nvox, data.begin());
        }
        else if (datatype == DT_COMPLEX128)
            memcpy(this->image->data, COMPLEX(object), dataSize);
        else
            memcpy(this->image->data, REAL(object), dataSize);
    }
    else
        this->image->data = NULL;
    
    if (object.hasAttribute("pixdim"))
    {
        const std::vector<float> pixdimVector = object.attr("pixdim");
        const int pixdimLength = pixdimVector.size();
        for (int i=0; i<std::min(pixdimLength,nDims); i++)
            this->image->pixdim[i+1] = pixdimVector[i];
    }
    
    if (object.hasAttribute("pixunits"))
    {
        const std::vector<std::string> pixunitsVector = object.attr("pixunits");
        setPixunits(pixunitsVector);
    }
}

inline void NiftiImage::initFromDims (const std::vector<int> &dim, const int datatype)
{
    const int nDims = std::min(7, int(dim.size()));
    int dims[8] = { nDims, 0, 0, 0, 0, 0, 0, 0 };
    std::copy(dim.begin(), dim.begin() + nDims, &dims[1]);
    acquire(nifti_make_new_nim(dims, datatype, 1));
    
    if (image == NULL)
        throw std::runtime_error("Failed to create image from scratch");
}

inline NiftiImage::NiftiImage (const SEXP object, const bool readData, const bool readOnly)
    : image(NULL), refCount(NULL)
{
    Rcpp::RObject imageObject(object);
    bool resolved = false;
    
    if (imageObject.hasAttribute(".nifti_image_ptr"))
    {
        Rcpp::XPtr<NiftiImage> imagePtr(SEXP(imageObject.attr(".nifti_image_ptr")));
        NiftiImage *ptr = imagePtr;
        if (ptr != NULL)
        {
            if (MAYBE_SHARED(object) && !readOnly)
                copy(*ptr);
            else
                acquire(*ptr);
            resolved = true;
            
            if (imageObject.hasAttribute("dim"))
                update(imageObject);
        }
        else if (Rf_isString(object))
            throw std::runtime_error("Internal image is not valid");
        else
            Rf_warning("Ignoring invalid internal pointer");
    }
    
    if (!resolved)
    {
        if (Rf_isNull(object))
            acquire(NULL);
        else if (Rf_isString(object))
        {
            const std::string path = Rcpp::as<std::string>(object);
            acquire(nifti_image_read(internal::stringToPath(path), readData));
            if (this->image == NULL)
                throw std::runtime_error("Failed to read image from path " + path);
        }
        else if (imageObject.inherits("nifti"))
            initFromNiftiS4(imageObject, readData);
        else if (imageObject.inherits("anlz"))
            throw std::runtime_error("Cannot currently convert objects of class \"anlz\"");
        else if (imageObject.inherits("MriImage"))
            initFromMriImage(imageObject, readData);
        else if (Rf_isVectorList(object))
            initFromList(imageObject);
        else if (imageObject.hasAttribute("dim"))
            initFromArray(imageObject, readData);
        else if (imageObject.hasAttribute("class"))
            throw std::runtime_error("Cannot convert object of class \"" + Rcpp::as<std::string>(imageObject.attr("class")) + "\" to a nifti_image");
        else
            throw std::runtime_error("Cannot convert unclassed non-array object");
    }
    
    if (this->image != NULL)
        nifti_update_dims_from_array(this->image);
    
#ifndef NDEBUG
    Rc_printf("Creating NiftiImage with pointer %p (from SEXP)\n", this->image);
#endif
}

#endif // USING_R

inline NiftiImage::NiftiImage (const std::vector<int> &dim, const int datatype)
    : image(NULL), refCount(NULL)
{
    initFromDims(dim, datatype);
#ifndef NDEBUG
    Rc_printf("Creating NiftiImage with pointer %p (from dims)\n", this->image);
#endif
}

inline NiftiImage::NiftiImage (const std::vector<int> &dim, const std::string &datatype)
    : image(NULL), refCount(NULL)
{
    initFromDims(dim, internal::stringToDatatype(datatype));
#ifndef NDEBUG
    Rc_printf("Creating NiftiImage with pointer %p (from dims)\n", this->image);
#endif
}

inline NiftiImage::NiftiImage (const std::string &path, const bool readData)
    : image(NULL), refCount(NULL)
{
    acquire(nifti_image_read(internal::stringToPath(path), readData));
    
    if (image == NULL)
        throw std::runtime_error("Failed to read image from path " + path);
    
#ifndef NDEBUG
    Rc_printf("Creating NiftiImage with pointer %p (from string)\n", this->image);
#endif
}

inline NiftiImage::NiftiImage (const std::string &path, const std::vector<int> &volumes)
    : image(NULL), refCount(NULL)
{
    if (volumes.empty())
        throw std::runtime_error("The vector of volumes is empty");
    
    nifti_brick_list brickList;
    acquire(nifti_image_read_bricks(internal::stringToPath(path), volumes.size(), &volumes[0], &brickList));
    if (image == NULL)
        throw std::runtime_error("Failed to read image from path " + path);
    
    size_t brickSize = image->nbyper * image->nx * image->ny * image->nz;
    image->data = calloc(1, nifti_get_volsize(image));
    for (int i=0; i<brickList.nbricks; i++)
        memcpy((char *) image->data + i * brickSize, brickList.bricks[i], brickSize);
    nifti_free_NBL(&brickList);
    
#ifndef NDEBUG
    Rc_printf("Creating NiftiImage with pointer %p (from string and volume vector)\n", this->image);
#endif
}

inline void NiftiImage::updatePixdim (const std::vector<float> &pixdim)
{
    const int nDims = image->dim[0];
    const std::vector<float> origPixdim(image->pixdim+1, image->pixdim+4);
    
    for (int i=1; i<8; i++)
        image->pixdim[i] = 0.0;
    
    const int pixdimLength = pixdim.size();
    for (int i=0; i<std::min(pixdimLength,nDims); i++)
        image->pixdim[i+1] = pixdim[i];
    
    if (!std::equal(origPixdim.begin(), origPixdim.begin() + std::min(3,nDims), pixdim.begin()))
    {
        mat33 scaleMatrix;
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                if (i != j)
                    scaleMatrix.m[i][j] = 0.0;
                else if (i >= pixdimLength)
                    scaleMatrix.m[i][j] = 1.0;
                else
                    scaleMatrix.m[i][j] = pixdim[i] / origPixdim[i];
            }
        }
        
        if (image->qform_code > 0)
        {
            mat33 prod = nifti_mat33_mul(scaleMatrix, internal::topLeftCorner(image->qto_xyz));
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                    image->qto_xyz.m[i][j] = prod.m[i][j];
            }
            image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
            nifti_mat44_to_quatern(image->qto_xyz, &image->quatern_b, &image->quatern_c, &image->quatern_d, &image->qoffset_x, &image->qoffset_y, &image->qoffset_z, NULL, NULL, NULL, &image->qfac);
        }
        
        if (image->sform_code > 0)
        {
            mat33 prod = nifti_mat33_mul(scaleMatrix, internal::topLeftCorner(image->sto_xyz));
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                    image->sto_xyz.m[i][j] = prod.m[i][j];
            }
            image->sto_ijk = nifti_mat44_inverse(image->sto_xyz);
        }
    }
}

inline void NiftiImage::setPixunits (const std::vector<std::string> &pixunits)
{
    for (size_t i=0; i<pixunits.size(); i++)
    {
        if (pixunits[i] == "m")
            image->xyz_units = NIFTI_UNITS_METER;
        else if (pixunits[i] == "mm")
            image->xyz_units = NIFTI_UNITS_MM;
        else if (pixunits[i] == "um")
            image->xyz_units = NIFTI_UNITS_MICRON;
        else if (pixunits[i] == "s")
            image->time_units = NIFTI_UNITS_SEC;
        else if (pixunits[i] == "ms")
            image->time_units = NIFTI_UNITS_MSEC;
        else if (pixunits[i] == "us")
            image->time_units = NIFTI_UNITS_USEC;
        else if (pixunits[i] == "Hz")
            image->time_units = NIFTI_UNITS_HZ;
        else if (pixunits[i] == "ppm")
            image->time_units = NIFTI_UNITS_PPM;
        else if (pixunits[i] == "rad/s")
            image->time_units = NIFTI_UNITS_RADS;
    }
}

inline NiftiImage & NiftiImage::rescale (const std::vector<float> &scales)
{
    std::vector<float> pixdim(image->pixdim+1, image->pixdim+4);
    
    for (int i=0; i<std::min(3, int(scales.size())); i++)
    {
        if (scales[i] != 1.0)
        {
            pixdim[i] /= scales[i];
            image->dim[i+1] = static_cast<int>(std::floor(image->dim[i+1] * scales[i]));
        }
    }
    
    updatePixdim(pixdim);
    nifti_update_dims_from_array(image);
    
    // Data vector is now the wrong size, so drop it
    nifti_image_unload(image);
    
    image->scl_slope = 0.0;
    image->scl_inter = 0.0;
    
    return *this;
}

inline NiftiImage & NiftiImage::reorient (const int icode, const int jcode, const int kcode)
{
    if (this->isNull())
        return *this;
    if (image->qform_code == 0 && image->sform_code == 0)
    {
        Rf_warning("Image qform and sform codes are both zero, so it cannot be reoriented");
        return *this;
    }
    
    int used[6] = { 0, 0, 0, 0, 0, 0 };
    used[icode-1] = 1;
    used[jcode-1] = 1;
    used[kcode-1] = 1;
    if (used[0]+used[1] != 1 || used[2]+used[3] != 1 || used[4]+used[5] != 1)
        throw std::runtime_error("Each canonical axis should be used exactly once");
    
    const int codes[3] = { icode, jcode, kcode };
    const mat44 native = this->xform();
    
    // Calculate the origin, which requires inverting the current xform
    // Here we use a simplified formula that exploits blockwise inversion and the nature of xforms
    internal::vec3 origin;
    for (int i=0; i<3; i++)
        origin.v[i] = native.m[i][3];
    origin = -internal::matrixVectorProduct(nifti_mat33_inverse(internal::topLeftCorner(native)), origin);
    
    // Create a target xform (rotation matrix only)
    mat33 target;
    for (int j=0; j<3; j++)
    {
        for (int i=0; i<3; i++)
            target.m[i][j] = 0.0;
        
        switch (codes[j])
        {
            case NIFTI_L2R: target.m[0][j] =  1.0; break;
            case NIFTI_R2L: target.m[0][j] = -1.0; break;
            case NIFTI_P2A: target.m[1][j] =  1.0; break;
            case NIFTI_A2P: target.m[1][j] = -1.0; break;
            case NIFTI_I2S: target.m[2][j] =  1.0; break;
            case NIFTI_S2I: target.m[2][j] = -1.0; break;
        }
    }
    
    // Extract (inverse of) canonical axis matrix from native xform
    int nicode, njcode, nkcode;
    nifti_mat44_to_orientation(native, &nicode, &njcode, &nkcode);
    int ncodes[3] = { nicode, njcode, nkcode };
    mat33 nativeAxesTransposed;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
            nativeAxesTransposed.m[i][j] = 0.0;

        switch (ncodes[i])
        {
            case NIFTI_L2R: nativeAxesTransposed.m[i][0] =  1.0; break;
            case NIFTI_R2L: nativeAxesTransposed.m[i][0] = -1.0; break;
            case NIFTI_P2A: nativeAxesTransposed.m[i][1] =  1.0; break;
            case NIFTI_A2P: nativeAxesTransposed.m[i][1] = -1.0; break;
            case NIFTI_I2S: nativeAxesTransposed.m[i][2] =  1.0; break;
            case NIFTI_S2I: nativeAxesTransposed.m[i][2] = -1.0; break;
        }
    }
    
    // Check for no-op case
    if (icode == nicode && jcode == njcode && kcode == nkcode)
        return *this;
    
    // The transform is t(approx_old_xform) %*% target_xform
    // The new xform is old_xform %*% transform
    // NB: "transform" is really 4x4, but the last row is simple and the last column is filled below
    mat33 transform = nifti_mat33_mul(nativeAxesTransposed, target);
    mat44 result;
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<3; j++)
            result.m[i][j] = native.m[i][0] * transform.m[0][j] + native.m[i][1] * transform.m[1][j] + native.m[i][2] * transform.m[2][j];
        
        result.m[3][i] = i == 3 ? 1.0 : 0.0;
    }
    
    // Extract the mapping between dimensions and the signs
    // These vectors are all indexed in the target space, except "revsigns"
    int locs[3], signs[3], newdim[3], revsigns[3];
    float newpixdim[3];
    double maxes[3] = { R_NegInf, R_NegInf, R_NegInf };
    internal::vec3 offset;
    for (int j=0; j<3; j++)
    {
        // Find the largest absolute value in each column, which gives the old dimension corresponding to each new dimension
        for (int i=0; i<3; i++)
        {
            const double value = static_cast<double>(transform.m[i][j]);
            if (fabs(value) > maxes[j])
            {
                maxes[j] = fabs(value);
                signs[j] = value > 0.0 ? 1 : -1;
                locs[j] = i;
            }
        }
        
        // Obtain the sign for the reverse mapping
        revsigns[locs[j]] = signs[j];
        
        // Permute dim and pixdim
        newdim[j] = image->dim[locs[j]+1];
        newpixdim[j] = image->pixdim[locs[j]+1];
        
        // Flip and/or permute the origin
        if (signs[j] < 0)
            offset.v[j] = image->dim[locs[j]+1] - origin.v[locs[j]] - 1.0;
        else
            offset.v[j] = origin.v[locs[j]];
    }
    
    // Convert the origin back to an xform offset and insert it
    offset = -internal::matrixVectorProduct(internal::topLeftCorner(result), offset);
    for (int i=0; i<3; i++)
        result.m[i][3] = offset.v[i];
    
    // Update the xforms with nonzero codes
    if (image->qform_code > 0)
    {
        image->qto_xyz = result;
        image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
        nifti_mat44_to_quatern(image->qto_xyz, &image->quatern_b, &image->quatern_c, &image->quatern_d, &image->qoffset_x, &image->qoffset_y, &image->qoffset_z, NULL, NULL, NULL, &image->qfac);
    }
    if (image->sform_code > 0)
    {
        image->sto_xyz = result;
        image->sto_ijk = nifti_mat44_inverse(image->sto_xyz);
    }
    
    // Calculate strides: the step in target space associated with each dimension in source space
    ptrdiff_t strides[3];
    strides[locs[0]] = 1;
    strides[locs[1]] = strides[locs[0]] * image->dim[locs[0]+1];
    strides[locs[2]] = strides[locs[1]] * image->dim[locs[1]+1];
    
    // Permute the data (if present)
    if (image->data != NULL)
    {    
        size_t volSize = size_t(image->nx * image->ny * image->nz);
        size_t nVolumes = std::max(size_t(1), image->nvox / volSize);
        
        const NiftiImageData oldData = this->data();
        NiftiImageData newData(oldData);
        
        // Where the sign is negative we need to start at the end of the dimension
        size_t volStart = 0;
        for (int i=0; i<3; i++)
        {
            if (revsigns[i] < 0)
                volStart += (image->dim[i+1] - 1) * strides[i];
        }
        
        // Iterate over the data and place it into a new vector
        NiftiImageData::Iterator it = oldData.begin();
        for (size_t v=0; v<nVolumes; v++)
        {
            for (int k=0; k<image->nz; k++)
            {
                ptrdiff_t offset = k * strides[2] * revsigns[2];
                for (int j=0; j<image->ny; j++)
                {
                    for (int i=0; i<image->nx; i++)
                    {
                        newData[volStart + offset] = *it++;
                        offset += strides[0] * revsigns[0];
                    }
                    offset += strides[1] * revsigns[1] - image->nx * strides[0] * revsigns[0];
                }
            }
            volStart += volSize;
        }
        
        // Vector data needs to be reoriented to match the xform
        if (image->intent_code == NIFTI_INTENT_VECTOR && image->dim[5] == 3)
        {
            internal::vec3 oldVec;
            const size_t supervolSize = volSize * image->nt;
            NiftiImageData::Iterator it = newData.begin();
            for (size_t i=0; i<supervolSize; i++, ++it)
            {
                for (int j=0; j<3; j++)
                    oldVec.v[j] = double(*(it + j*supervolSize));
                const internal::vec3 newVec = internal::matrixVectorProduct(transform, oldVec);
                for (int j=0; j<3; j++)
                    *(it + j*supervolSize) = newVec.v[j];
            }
        }
        
        // Replace the existing data in the image
        this->replaceData(newData);
    }
    
    // Copy new dims and pixdims in
    // NB: Old dims are used above, so this must happen last
    std::copy(newdim, newdim+3, image->dim+1);
    std::copy(newpixdim, newpixdim+3, image->pixdim+1);
    nifti_update_dims_from_array(image);
    
    return *this;
}

inline NiftiImage & NiftiImage::reorient (const std::string &orientation)
{
    if (orientation.length() != 3)
        throw std::runtime_error("Orientation string should have exactly three characters");
    
    int codes[3];
    for (int i=0; i<3; i++)
    {
        switch (orientation[i])
        {
            case 'r': case 'R': codes[i] = NIFTI_L2R; break;
            case 'l': case 'L': codes[i] = NIFTI_R2L; break;
            case 'a': case 'A': codes[i] = NIFTI_P2A; break;
            case 'p': case 'P': codes[i] = NIFTI_A2P; break;
            case 's': case 'S': codes[i] = NIFTI_I2S; break;
            case 'i': case 'I': codes[i] = NIFTI_S2I; break;
            
            default:
            throw std::runtime_error("Orientation string is invalid");
        }
    }
    
    return reorient(codes[0], codes[1], codes[2]);
}

#ifdef USING_R

inline NiftiImage & NiftiImage::update (const Rcpp::RObject &object)
{
    if (Rf_isVectorList(object))
    {
        Rcpp::List list(object);
        nifti_1_header *header = NULL;
        if (this->isNull())
        {
            header = nifti_make_new_header(NULL, DT_FLOAT64);
            internal::updateHeader(header, list, true);
        }
        else
        {
            header = (nifti_1_header *) calloc(1, sizeof(nifti_1_header));
            *header = nifti_convert_nim2nhdr(image);
            internal::updateHeader(header, list, true);
        }
        
        if (header != NULL)
        {
            // Retain the data pointer, but otherwise overwrite the stored object with one created from the header
            // The file names can't be preserved through the round-trip, so free them
            void *dataPtr = image->data;
            nifti_image *tempImage = nifti_convert_nhdr2nim(*header, NULL);
            
            if (image->fname != NULL)
                free(image->fname);
            if (image->iname != NULL)
                free(image->iname);
            
            memcpy(image, tempImage, sizeof(nifti_image));
            image->num_ext = 0;
            image->ext_list = NULL;
            image->data = dataPtr;
            
            nifti_image_free(tempImage);
            free(header);
        }
    }
    else if (object.hasAttribute("dim"))
    {
        for (int i=0; i<8; i++)
            image->dim[i] = 0;
        const std::vector<int> dimVector = object.attr("dim");
    
        const int nDims = std::min(7, int(dimVector.size()));
        image->dim[0] = nDims;
        for (int i=0; i<nDims; i++)
            image->dim[i+1] = dimVector[i];
    
        if (object.hasAttribute("pixdim"))
        {
            const std::vector<float> pixdimVector = object.attr("pixdim");
            updatePixdim(pixdimVector);
        }
    
        if (object.hasAttribute("pixunits"))
        {
            const std::vector<std::string> pixunitsVector = object.attr("pixunits");
            setPixunits(pixunitsVector);
        }
    
        // This NIfTI-1 library function clobbers dim[0] if the last dimension is unitary; we undo that here
        nifti_update_dims_from_array(image);
        image->dim[0] = image->ndim = nDims;
    
        image->datatype = NiftiImage::sexpTypeToNiftiType(object.sexp_type());
        if (object.inherits("rgbArray"))
        {
            const int channels = object.attr("channels");
            image->datatype = (channels == 4 ? DT_RGBA32 : DT_RGB24);
        }
        nifti_datatype_sizes(image->datatype, &image->nbyper, NULL);
    
        nifti_image_unload(image);
    
        const size_t dataSize = nifti_get_volsize(image);
        image->data = calloc(1, dataSize);
        if (image->datatype == DT_INT32 || image->datatype == DT_RGBA32)
            memcpy(image->data, INTEGER(object), dataSize);
        else if (image->datatype == DT_RGB24)
            std::copy(INTEGER(object), INTEGER(object)+image->nvox, this->data().begin());
        else if (image->datatype == DT_COMPLEX128)
            memcpy(image->data, COMPLEX(object), dataSize);
        else
            memcpy(image->data, REAL(object), dataSize);
    
        image->scl_slope = 0.0;
        image->scl_inter = 0.0;
    }
    
    return *this;
}

#endif // USING_R

inline mat44 NiftiImage::xform (const bool preferQuaternion) const
{
    if (image == NULL)
    {
        mat44 matrix;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
                matrix.m[i][j] = 0.0;
        }
        return matrix;
    }
    else if (image->qform_code <= 0 && image->sform_code <= 0)
    {
        // No qform or sform so use pixdim (NB: other software may assume differently)
        mat44 matrix;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
            {
                if (i != j)
                    matrix.m[i][j] = 0.0;
                else if (i == 3)
                    matrix.m[3][3] = 1.0;
                else
                    matrix.m[i][j] = (image->pixdim[i+1]==0.0 ? 1.0 : image->pixdim[i+1]);
            }
        }
        return matrix;
    }
    else if ((preferQuaternion && image->qform_code > 0) || image->sform_code <= 0)
        return image->qto_xyz;
    else
        return image->sto_xyz;
}

template <typename TargetType>
inline std::vector<TargetType> NiftiImage::Block::getData (const bool useSlope) const
{
    NiftiImageData data = this->data();
    if (!useSlope)
        data = data.unscaled();
    
    if (image.isNull() || data.isEmpty())
        return std::vector<TargetType>();
    else
    {
        std::vector<TargetType> result(data.size());
        std::copy(data.begin(), data.end(), result.begin());
        return result;
    }
}

template <typename TargetType>
inline std::vector<TargetType> NiftiImage::getData (const bool useSlope) const
{
    NiftiImageData data = this->data();
    if (!useSlope)
        data = data.unscaled();
    
    if (this->isNull() || data.isEmpty())
        return std::vector<TargetType>();
    else
    {
        std::vector<TargetType> result(data.size());
        std::copy(data.begin(), data.end(), result.begin());
        return result;
    }
}

inline NiftiImage & NiftiImage::changeDatatype (const int datatype, const bool useSlope)
{
    if (this->isNull() || image->datatype == datatype)
        return *this;
    
    if (useSlope && this->isDataScaled())
        throw std::runtime_error("Resetting the slope and intercept for an image with them already set is not supported");
    
    const NiftiImageData data(useSlope ? this->data() : this->data().unscaled(), datatype);
    return replaceData(data);
}

inline NiftiImage & NiftiImage::changeDatatype (const std::string &datatype, const bool useSlope)
{
    return changeDatatype(internal::stringToDatatype(datatype), useSlope);
}

template <typename SourceType>
inline NiftiImage & NiftiImage::replaceData (const std::vector<SourceType> &data, const int datatype)
{
    replaceData(NiftiImageData(data.begin(), data.end(), datatype));
    return *this;
}

inline NiftiImage & NiftiImage::replaceData (const NiftiImageData &data)
{
    if (this->isNull())
        return *this;
    else if (data.isEmpty())
    {
        nifti_image_unload(image);
        return *this;
    }
    else if (data.length() != image->nvox)
        throw std::runtime_error("New data length does not match the number of voxels in the image");
    
    // Copy the data
    NiftiImageData copy = data;
    nifti_image_unload(image);
    image->data = copy.blob();
    image->datatype = copy.datatype();
    image->scl_slope = static_cast<float>(copy.slope);
    image->scl_inter = static_cast<float>(copy.intercept);
    nifti_datatype_sizes(image->datatype, &image->nbyper, &image->swapsize);
    
    double min, max;
    copy.minmax(&min, &max);
    image->cal_min = static_cast<float>(min);
    image->cal_max = static_cast<float>(max);
    
    copy.disown();
    
    return *this;
}

inline std::pair<std::string,std::string> NiftiImage::toFile (const std::string fileName, const int datatype) const
{
    const bool changingDatatype = (datatype != DT_NONE && !this->isNull() && datatype != image->datatype);
    
    // Copy the source image only if the datatype will be changed
    NiftiImage imageToWrite(*this, changingDatatype);
    
    if (changingDatatype)
        imageToWrite.changeDatatype(datatype, true);
    
    const int status = nifti_set_filenames(imageToWrite, internal::stringToPath(fileName), false, true);
    if (status != 0)
        throw std::runtime_error("Failed to set filenames for NIfTI object");
    nifti_image_write(imageToWrite);
    
    return std::pair<std::string,std::string>(std::string(imageToWrite->fname), std::string(imageToWrite->iname));
}

inline std::pair<std::string,std::string> NiftiImage::toFile (const std::string fileName, const std::string &datatype) const
{
    return toFile(fileName, internal::stringToDatatype(datatype));
}

#ifdef USING_R

inline Rcpp::RObject NiftiImage::toArray () const
{
    Rcpp::RObject array;
    
    if (this->isNull())
        return array;
    else
    {
        NiftiImageData data = this->data();
        if (data.isEmpty())
        {
            Rf_warning("Internal image contains no data - filling array with NAs");
            array = Rcpp::LogicalVector(image->nvox, NA_LOGICAL);
        }
        else if (data.isComplex())
            array = Rcpp::ComplexVector(data.begin(), data.end());
        else if (data.isFloatingPoint() || data.isScaled())
            array = Rcpp::NumericVector(data.begin(), data.end());
        else
            array = Rcpp::IntegerVector(data.begin(), data.end());
    
        internal::addAttributes(array, *this, true, true, false);
        if (data.isRgb())
        {
            array.attr("class") = Rcpp::CharacterVector::create("niftiImage", "rgbArray", "array");
            array.attr("channels") = (data.datatype() == DT_RGBA32 ? 4 : 3);
        }
        else
            array.attr("class") = Rcpp::CharacterVector::create("niftiImage", "array");
        return array;
    }
}

inline Rcpp::RObject NiftiImage::toPointer (const std::string label) const
{
    if (this->isNull())
        return Rcpp::RObject();
    else
    {
        Rcpp::RObject string = Rcpp::wrap(label);
        internal::addAttributes(string, *this, false);
        string.attr("class") = Rcpp::CharacterVector::create("internalImage", "niftiImage");
        return string;
    }
}

inline Rcpp::RObject NiftiImage::toArrayOrPointer (const bool internal, const std::string label) const
{
    return (internal ? toPointer(label) : toArray());
}

inline Rcpp::RObject NiftiImage::headerToList () const
{
    if (this->image == NULL)
        return Rcpp::RObject();
    
    nifti_1_header header = nifti_convert_nim2nhdr(this->image);
    Rcpp::List result;
    
    result["sizeof_hdr"] = header.sizeof_hdr;
    
    result["dim_info"] = int(header.dim_info);
    result["dim"] = std::vector<short>(header.dim, header.dim+8);
    
    result["intent_p1"] = header.intent_p1;
    result["intent_p2"] = header.intent_p2;
    result["intent_p3"] = header.intent_p3;
    result["intent_code"] = header.intent_code;
    
    result["datatype"] = header.datatype;
    result["bitpix"] = header.bitpix;
    
    result["slice_start"] = header.slice_start;
    result["pixdim"] = std::vector<float>(header.pixdim, header.pixdim+8);
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
    
    Rcpp::List strings;
    strings["datatype"] = nifti_datatype_string(header.datatype);
    strings["intent_code"] = nifti_intent_string(header.intent_code);
    strings["qform_code"] = nifti_xform_string(header.qform_code);
    strings["sform_code"] = nifti_xform_string(header.sform_code);
    strings["slice_code"] = nifti_slice_string(header.slice_code);
    
    internal::addAttributes(result, *this, false, false);
    result.attr("class") = Rcpp::CharacterVector::create("niftiHeader");
    result.attr("strings") = strings;
    
    return result;
}

#endif // USING_R

#endif
