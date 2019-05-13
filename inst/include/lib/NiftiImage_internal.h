#ifndef _NIFTI_IMAGE_INTERNAL_H_
#define _NIFTI_IMAGE_INTERNAL_H_

namespace internal {

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

template <typename TargetType>
class DataConverter
{
public:
    enum ConversionMode { CastMode, ScaleMode, IndexMode };
    
protected:
    ConversionMode mode;
    double slope, intercept;
    
    template <typename Type>
    static bool lessThan (Type a, Type b)
    {
        return (!isNaN(a) && !isNaN(b) && a < b);
    }
    
public:
    DataConverter (const ConversionMode mode = CastMode)
        : mode(mode), slope(0.0), intercept(0.0) {}
    
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
            const double dataMin = static_cast<double>(*std::min_element(source, source + length, DataConverter::lessThan<SourceType>));
            const double dataMax = static_cast<double>(*std::max_element(source, source + length, DataConverter::lessThan<SourceType>));
            
            // If the source type is floating-point but values are in range, we will just round them
            if (dataMin < typeMin || dataMax > typeMax)
            {
                slope = (dataMax - dataMin) / (typeMax - typeMin);
                intercept = dataMin - slope * typeMin;
            }
            else
            {
                slope = 1.0;
                intercept = 0.0;
            }
        }
        return *this;
    }
    
    template <typename SourceType>
    TargetType convertValue (const SourceType value) const
    {
        // If NaNs aren't going to make it across, treat them as zero
        if (isNaN(value) && !std::numeric_limits<TargetType>::has_quiet_NaN)
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
    TargetType operator() (const SourceType value) const
    {
        return convertValue(value);
    }
};

template <typename TargetType, class OutputIterator>
inline void convertData (void *source, const short datatype, const size_t length, OutputIterator target, const ptrdiff_t offset = 0, DataConverter<TargetType> *converter = NULL)
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
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_INT16:
        {
            int16_t *castSource = static_cast<int16_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_INT32:
        {
            int32_t *castSource = static_cast<int32_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_FLOAT32:
        {
            float *castSource = static_cast<float *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_FLOAT64:
        {
            double *castSource = static_cast<double *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_INT8:
        {
            int8_t *castSource = static_cast<int8_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_UINT16:
        {
            uint16_t *castSource = static_cast<uint16_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_UINT32:
        {
            uint32_t *castSource = static_cast<uint32_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_INT64:
        {
            int64_t *castSource = static_cast<int64_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
            break;
        }
        
        case DT_UINT64:
        {
            uint64_t *castSource = static_cast<uint64_t *>(source) + offset;
            if (converter->getMode() == DataConverter<TargetType>::IndexMode)
                converter->calibrate(castSource, length);
            std::transform(castSource, castSource + length, target, *converter);
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
        std::transform(begin, end, static_cast<uint8_t *>(target), DataConverter<uint8_t>());
        break;
        
        case DT_INT16:
        std::transform(begin, end, static_cast<int16_t *>(target), DataConverter<int16_t>());
        break;
        
        case DT_INT32:
        std::transform(begin, end, static_cast<int32_t *>(target), DataConverter<int32_t>());
        break;
        
        case DT_FLOAT32:
        std::transform(begin, end, static_cast<float *>(target), DataConverter<float>());
        break;
        
        case DT_FLOAT64:
        std::transform(begin, end, static_cast<double *>(target), DataConverter<double>());
        break;
        
        case DT_INT8:
        std::transform(begin, end, static_cast<int8_t *>(target), DataConverter<int8_t>());
        break;
        
        case DT_UINT16:
        std::transform(begin, end, static_cast<uint16_t *>(target), DataConverter<uint16_t>());
        break;
        
        case DT_UINT32:
        std::transform(begin, end, static_cast<uint32_t *>(target), DataConverter<uint32_t>());
        break;
        
        case DT_INT64:
        std::transform(begin, end, static_cast<int64_t *>(target), DataConverter<int64_t>());
        break;
        
        case DT_UINT64:
        std::transform(begin, end, static_cast<uint64_t *>(target), DataConverter<uint64_t>());
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

inline void addAttributes (const SEXP pointer, const NiftiImage &source, const bool realDim = true, const bool includeXptr = true)
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
        Rcpp::XPtr<NiftiImage> xptr(new NiftiImage(source,false));
        object.attr(".nifti_image_ptr") = xptr;
    }
}

#endif // USING_R

} // namespace

#endif
