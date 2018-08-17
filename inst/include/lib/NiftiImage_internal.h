#ifndef _NIFTI_IMAGE_INTERNAL_H_
#define _NIFTI_IMAGE_INTERNAL_H_

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

double roundEven (const double value)
{
    if (ISNAN(value))
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

template <typename TargetType>
class DataConverter
{
public:
    enum ConversionMode { CastMode, ScaleMode, IndexMode };
    
protected:
    double slope, intercept;
    ConversionMode mode;
    
public:
    DataConverter ()
        : mode(CastMode), slope(1.0), intercept(0.0) {}
    
    DataConverter (const float slope, const float intercept)
        : mode(ScaleMode), slope(static_cast<double>(slope)), intercept(static_cast<double>(intercept))
    {
        if (slope == 0.0 || (slope == 1.0 && intercept == 0.0))
            mode = CastMode;
    }
    
    double getSlope () const { return slope; }
    double getIntercept () const { return intercept; }
    ConversionMode getMode () const { return mode; }
    
    DataConverter & setMode (const ConversionMode newMode)
    {
        mode = newMode;
        return *this;
    }
    
    template <typename SourceType>
    DataConverter & calibrate (const SourceType *source, const size_t length)
    {
        if (std::numeric_limits<TargetType>::is_integer)
        {
            const double typeMin = static_cast<double>(std::numeric_limits<TargetType>::min());
            const double typeMax = static_cast<double>(std::numeric_limits<TargetType>::max());
            const double dataMin = static_cast<double>(*std::min_element(source, source + length));
            const double dataMax = static_cast<double>(*std::max_element(source, source + length));
            if (!std::numeric_limits<SourceType>::is_integer || dataMin < typeMin || dataMax > typeMax)
            {
                slope = (dataMax - dataMin) / (typeMax - typeMin);
                intercept = dataMin - slope * typeMin;
            }
        }
        return *this;
    }
    
    template <typename SourceType>
    TargetType indexValue (const SourceType value) const
    {
        // If NaNs aren't going to make it across, treat them as zero
        if (std::numeric_limits<SourceType>::has_quiet_NaN && !std::numeric_limits<TargetType>::has_quiet_NaN && ISNAN(value))
            return static_cast<TargetType>(roundEven(-intercept / slope));
        else
            return static_cast<TargetType>(roundEven((static_cast<double>(value) - intercept) / slope));
    }
    
    template <typename SourceType>
    TargetType scaleValue (const SourceType value) const
    {
        // If NaNs aren't going to make it across, treat them as zero
        if (std::numeric_limits<SourceType>::has_quiet_NaN && !std::numeric_limits<TargetType>::has_quiet_NaN && ISNAN(value))
            return TargetType(0);
        else
            return static_cast<TargetType>(static_cast<double>(value) * slope + intercept);
    }
    
    template <typename SourceType>
    TargetType convertValue (const SourceType value) const
    {
        // If NaNs aren't going to make it across, treat them as zero
        if (std::numeric_limits<SourceType>::has_quiet_NaN && !std::numeric_limits<TargetType>::has_quiet_NaN && ISNAN(value))
            return (mode == IndexMode ? static_cast<TargetType>(roundEven(-intercept / slope)) : TargetType(0));
        else if (mode == CastMode)
            return static_cast<TargetType>(value);
        else if (mode == ScaleMode)
            return static_cast<TargetType>(static_cast<double>(value) * slope + intercept);
        else if (mode == IndexMode)
            return static_cast<TargetType>(roundEven((static_cast<double>(value) - intercept) / slope));
        else
            return TargetType(0);
    }
    
    template <typename SourceType>
    void run (const SourceType *source, const size_t length, TargetType *target) const
    {
        for (size_t i=0; i<length; i++)
            target[i] = convertValue(source[i]);
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

template <typename TargetType>
inline void convertData (void *source, const short datatype, const size_t length, TargetType *target, const ptrdiff_t offset = 0, DataConverter<TargetType> *converter = NULL)
{
    if (source == NULL)
        return;
    
    const bool willFreeConverter = (converter == NULL);
    if (willFreeConverter)
        converter = new DataConverter<TargetType>;
    
    switch (datatype)
    {
        case DT_UINT8:
        {
            uint8_t *castSource = static_cast<uint8_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            converter->run(castSource, length, target);
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
    
    if (willFreeConverter)
        delete converter;
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

#ifdef USING_R

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
        for (size_t i=0; i<std::min(dim.size(),size_t(8)); i++)
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
        for (size_t i=0; i<std::min(pixdim.size(),size_t(8)); i++)
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
        for (size_t i=0; i<std::min(srow_x.size(),size_t(4)); i++)
            header->srow_x[i] = srow_x[i];
    }
    if (names.count("srow_y") == 1)
    {
        std::vector<float> srow_y = list["srow_y"];
        for (size_t i=0; i<std::min(srow_y.size(),size_t(4)); i++)
            header->srow_y[i] = srow_y[i];
    }
    if (names.count("srow_z") == 1)
    {
        std::vector<float> srow_z = list["srow_z"];
        for (size_t i=0; i<std::min(srow_z.size(),size_t(4)); i++)
            header->srow_z[i] = srow_z[i];
    }
    
    if (names.count("intent_name") == 1)
        strcpy(header->intent_name, Rcpp::as<std::string>(list["intent_name"]).substr(0,15).c_str());
    if (names.count("magic") == 1)
        strcpy(header->magic, Rcpp::as<std::string>(list["magic"]).substr(0,3).c_str());
}

#endif // USING_R

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

#ifdef USING_R

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
            convertData(source->data, source->datatype, source->nvox, INTEGER(array));
        else if (SexpType == REALSXP)
            convertData(source->data, source->datatype, source->nvox, REAL(array));
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
    
    NiftiImage *wrappedSource = new NiftiImage(source, true);
    wrappedSource->setPersistence(true);
    Rcpp::XPtr<NiftiImage> xptr(wrappedSource);
    R_RegisterCFinalizerEx(SEXP(xptr), &finaliseNiftiImage, FALSE);
    object.attr(".nifti_image_ptr") = xptr;
}

#endif // USING_R

} // namespace

#endif
