#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_


#ifdef USING_R

#include <Rcpp.h>

// Defined since R 3.1.0, according to Tomas Kalibera, but there's no reason to break compatibility with 3.0.x
#ifndef MAYBE_SHARED
#define MAYBE_SHARED(x) (NAMED(x) > 1)
#endif

#else

#define R_NegInf -INFINITY

#include <stdint.h>
#include <cstddef>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <locale>
#include <limits>

#endif


#include "niftilib/nifti1_io.h"

/**
 * @mainpage RNifti: Fast R and C++ Access to NIfTI Images
 * A more extensive overview of the \c RNifti package, and its usage from R, is provided on the
 * package's GitHub page at \c https://github.com/jonclayden/RNifti. The primary role of these
 * pages is to document the \ref RNifti::NiftiImage and \ref RNifti::NiftiImageData C++ classes
 * for package developers linking to \c RNifti.
**/

namespace RNifti {

/**
 * Wrapper class encapsulating a NIfTI data blob, with responsibility for handling data scaling
 * and polymorphism. This class provides read/write data access, iterators, etc., which internally
 * handle conversion to and from the data's native type. It can be linked to the data in a
 * \c nifti_image or used independently.
 * @author Jon Clayden (<code@clayden.org>)
**/
class NiftiImageData
{
public:
    double slope;                       /**< The slope term used to scale data values. Ignored if zero. */
    double intercept;                   /**< The intercept term used to scale data values */
    
protected:
    /**
     * Abstract inner class defining the type-specific functions required in concrete subclasses
    **/
    struct TypeHandler
    {
        virtual ~TypeHandler() {}
        virtual size_t size () const { return 0; }
        virtual bool hasNaN () const { return false; }
        virtual double getDouble (void *ptr) const = 0;
        virtual int getInt (void *ptr) const = 0;
        virtual void setDouble (void *ptr, const double value) const = 0;
        virtual void setInt (void *ptr, const int value) const = 0;
        virtual void minmax (void *ptr, const size_t length, double *min, double *max) const = 0;
    };
    
    /**
     * Concrete inner class template defining behaviour specific to individual data types
    **/
    template <typename Type>
    struct ConcreteTypeHandler : public TypeHandler
    {
        size_t size () const { return (sizeof(Type)); }
        bool hasNaN () const { return std::numeric_limits<Type>::has_quiet_NaN; }
        double getDouble (void *ptr) const { return static_cast<double>(*static_cast<Type*>(ptr)); }
        int getInt (void *ptr) const { return static_cast<int>(*static_cast<Type*>(ptr)); }
        void setDouble (void *ptr, const double value) const { *(static_cast<Type*>(ptr)) = static_cast<Type>(value); }
        void setInt (void *ptr, const int value) const { *(static_cast<Type*>(ptr)) = static_cast<Type>(value); }
        void minmax (void *ptr, const size_t length, double *min, double *max) const;
    };
    
    /**
     * Create a concrete type handler appropriate to the datatype code stored with the data
     * @return The newly allocated type handler, or \c NULL
     * @exception runtime_error If the current datatype is unsupported
    **/
    TypeHandler * createHandler ()
    {
        if (_datatype == DT_NONE)
            return NULL;
        
        switch (_datatype)
        {
            case DT_UINT8:   return new ConcreteTypeHandler<uint8_t>();  break;
            case DT_INT16:   return new ConcreteTypeHandler<int16_t>();  break;
            case DT_INT32:   return new ConcreteTypeHandler<int32_t>();  break;
            case DT_FLOAT32: return new ConcreteTypeHandler<float>();    break;
            case DT_FLOAT64: return new ConcreteTypeHandler<double>();   break;
            case DT_INT8:    return new ConcreteTypeHandler<int8_t>();   break;
            case DT_UINT16:  return new ConcreteTypeHandler<uint16_t>(); break;
            case DT_UINT32:  return new ConcreteTypeHandler<uint32_t>(); break;
            case DT_INT64:   return new ConcreteTypeHandler<int64_t>();  break;
            case DT_UINT64:  return new ConcreteTypeHandler<uint64_t>(); break;
            
            default:
            throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(_datatype)) + ")");
        }
    }
    
    void *dataPtr;                      /**< Opaque pointer to the underlying data blob */
    int _datatype;                      /**< Datatype code indicating the actual type of the elements */
    TypeHandler *handler;               /**< Type handler, which is created to match the datatype */
    size_t _length;                     /**< The number of data elements in the blob */
    bool owner;                         /**< An indicator of whether this object is responsible for cleaning up the data */
    
    /**
     * Initialiser method, used by constructors
     * @param data Pointer to a preallocated data blob, or \c NULL
     * @param length Number of elements in the blob
     * @param datatype NIfTI datatype code appropriate to the blob
     * @param slope Slope parameter for scaling values
     * @param intercept Intercept parameter for scaling values
     * @param alloc If \c true, the default, and \c data is \c NULL, memory will be allocated for
     *   the blob. If \c false, the blob will be \c NULL in this case
    **/
    void init (void *data, const size_t length, const int datatype, const double slope, const double intercept, const bool alloc = true)
    {
        this->_length = length;
        this->_datatype = datatype;
        this->slope = slope;
        this->intercept = intercept;
        
        owner = false;
        handler = createHandler();
        if (handler == NULL)
            dataPtr = NULL;
        else if (alloc && data == NULL)
        {
            dataPtr = calloc(length, handler->size());
            owner = true;
        }
        else
            dataPtr = data;
    }
    
    /**
     * Update the slope and intercept to cover the range of another data object. If the current
     * object's datatype can capture the required range without scaling, the slope and intercept
     * are simply reset
     * @param data Another data object
    **/
    void calibrateFrom (const NiftiImageData &data)
    {
        slope = 1.0;
        intercept = 0.0;
        
        if (this->isInteger())
        {
            double dataMin, dataMax, typeMin, typeMax;
            data.minmax(&dataMin, &dataMax);
            handler->minmax(NULL, 0, &typeMin, &typeMax);
            
            // If the source type is floating-point but values are in range, we will just round them
            if (dataMin < typeMin || dataMax > typeMax)
            {
                slope = (dataMax - dataMin) / (typeMax - typeMin);
                intercept = dataMin - (slope) * typeMin;
            }
        }
    }
    
public:
    /**
     * Inner class representing a single element in the data blob
    **/
    struct Element
    {
    private:
        const NiftiImageData &parent;
        void *ptr;
        
    public:
        /**
         * Primary constructor
         * @param parent A reference to the parent object
         * @param ptr An opaque pointer to the element. If \c NULL, the start of the data blob
         *   encapsulated by the parent will be used
        **/
        Element (const NiftiImageData &parent, void *ptr = NULL)
            : parent(parent)
        {
            this->ptr = (ptr == NULL ? parent.dataPtr : ptr);
        }
        
        /**
         * Copy assignment operator
         * @param value The value to assign. Any basic numeric type supported by NIfTI-1 is
         *   allowed, but \c int is used as an intermediate type for all integers, so values
         *   unrepresentable in a signed 32-bit integer may overflow
         * @return A reference to the callee
        **/
        template <typename SourceType>
        Element & operator= (const SourceType &value);
        
        /**
         * Copy assignment operator
         * @param other Another data element
         * @return A reference to the callee
        **/
        Element & operator= (const Element &other);
        
        /**
         * Implicit type-cast operator, suitable for implicit conversion to basic numeric types
        **/
        template <typename TargetType>
        operator TargetType() const
        {
            if (parent.isScaled())
                return TargetType(parent.handler->getDouble(ptr) * parent.slope + parent.intercept);
            else if (std::numeric_limits<TargetType>::is_integer)
                return TargetType(parent.handler->getInt(ptr));
            else
                return TargetType(parent.handler->getDouble(ptr));
        }
    };
    
    /**
     * Iterator type for \c NiftiImageData, with \c Element as its value type
    **/
    class Iterator : public std::iterator<std::random_access_iterator_tag, Element>
    {
    private:
        const NiftiImageData &parent;
        void *ptr;
        size_t step;
        
    public:
        /**
         * Primary constructor
         * @param parent A reference to the parent object
         * @param ptr An opaque pointer to the memory underpinning the iterator
         * @param step The increment between elements within the blob, in bytes. If zero, the
         *   default, the width associated with the stored datatype will be used.
        **/
        Iterator (const NiftiImageData &parent, void *ptr = NULL, const size_t step = 0)
            : parent(parent)
        {
            this->ptr = (ptr == NULL ? parent.dataPtr : ptr);
            this->step = (step == 0 ? parent.handler->size() : step);
        }
        
        /**
         * Copy constructor
         * @param other Another iterator
        **/
        Iterator (const Iterator &other)
            : parent(other.parent), ptr(other.ptr), step(other.step) {}
        
        Iterator & operator++ () { ptr = static_cast<char*>(ptr) + step; return *this; }
        Iterator operator++ (int) { Iterator copy(*this); ptr = static_cast<char*>(ptr) + step; return copy; }
        Iterator operator+ (ptrdiff_t n) const
        {
            void *newptr = static_cast<char*>(ptr) + (n * step);
            return Iterator(parent, newptr, step);
        }
        Iterator & operator-- () { ptr = static_cast<char*>(ptr) - step; return *this; }
        Iterator operator-- (int) { Iterator copy(*this); ptr = static_cast<char*>(ptr) - step; return copy; }
        Iterator operator- (ptrdiff_t n) const
        {
            void *newptr = static_cast<char*>(ptr) - (n * step);
            return Iterator(parent, newptr, step);
        }
        
        ptrdiff_t operator- (const Iterator &other) const
        {
            const ptrdiff_t difference = static_cast<char*>(ptr) - static_cast<char*>(other.ptr);
            return difference / step;
        }
        
        bool operator== (const Iterator &other) const { return (ptr==other.ptr && step==other.step); }
        bool operator!= (const Iterator &other) const { return (ptr!=other.ptr || step!=other.step); }
        bool operator> (const Iterator &other) const { return (ptr > other.ptr); }
        bool operator< (const Iterator &other) const { return (ptr < other.ptr); }
        
        const Element operator* () const { return Element(parent, ptr); }
        Element operator* () { return Element(parent, ptr); }
        const Element operator[] (const size_t i) const { return Element(parent, static_cast<char*>(ptr) + (i * step)); }
        Element operator[] (const size_t i) { return Element(parent, static_cast<char*>(ptr) + (i * step)); }
    };
    
    /**
     * Default constructor, creating an empty data object
    **/
    NiftiImageData ()
        : slope(1.0), intercept(0.0), dataPtr(NULL), _datatype(DT_NONE), handler(NULL), _length(0), owner(false) {}
    
    /**
     * Primary constructor
     * @param data A pointer to a pre-allocated data blob, or \c NULL. In the latter case, memory
     *   will be allocated by the object, and cleaned up at destruction unless it is disowned
     * @param length The number of elements in the blob
     * @param datatype The NIfTI datatype code corresponding to the type of the data elements
     * @param slope The slope parameter to use for data scaling, if any
     * @param intercept The intercept parameter to use for data scaling, if any
    **/
    NiftiImageData (void *data, const size_t length, const int datatype, const double slope = 1.0, const double intercept = 0.0)
    {
        init(data, length, datatype, slope, intercept);
    }
    
    /**
     * Convenience constructor for a \c nifti_image
     * @param image The image struct whose data the object will wrap
    **/
    NiftiImageData (nifti_image *image)
    {
        if (image == NULL)
            init(NULL, 0, DT_NONE, 0.0, 0.0, false);
        else
            init(image->data, image->nvox, image->datatype, static_cast<double>(image->scl_slope), static_cast<double>(image->scl_inter), false);
    }
    
    /**
     * Copy constructor with optional type conversion
     * @param source Another \c NiftiImageData object to copy data from
     * @param datatype The datatype to convert to, or \c DT_NONE, the default, for no conversion.
     *   If the range of the source data cannot be represented by the chosen type, the slope and
     *   intercept parameters will be set to adjust the range
    **/
    NiftiImageData (const NiftiImageData &source, const int datatype = DT_NONE)
    {
        init(NULL, source.length(), datatype == DT_NONE ? source.datatype() : datatype, source.slope, source.intercept);
        
        if (datatype == DT_NONE || datatype == source.datatype())
            memcpy(dataPtr, source.dataPtr, source.totalBytes());
        else
        {
            calibrateFrom(source);
            std::copy(source.begin(), source.end(), this->begin());
        }
    }
    
    /**
     * Iterator-based constructor
     * @param from Iterator type representing the start of the source data to be copied
     * @param to Iterator type representing the end of the source data to be copied
     * @param datatype The NIfTI datatype to use within the data blob
    **/
    template <class InputIterator>
    NiftiImageData (InputIterator from, InputIterator to, const int datatype)
    {
        const size_t length = static_cast<size_t>(std::distance(from, to));
        init(NULL, length, datatype, 1.0, 0.0);
        std::copy(from, to, this->begin());
    }
    
    /**
     * Destructor which frees the type handler, and the data blob if it is owned by this object
    **/
    virtual ~NiftiImageData ()
    {
        delete handler;
        if (owner)
            free(dataPtr);
    }
    
    /**
     * Copy assignment operator
     * @param source Another \c NiftiImageData object, from which the data and metadata are copied
     * @return A reference to the callee
    **/
    NiftiImageData & operator= (const NiftiImageData &source)
    {
        if (source.dataPtr != NULL)
        {
            // Free the old data, if we allocated it
            if (owner)
                free(dataPtr);
            init(NULL, source.length(), source.datatype(), source.slope, source.intercept);
            memcpy(dataPtr, source.dataPtr, source.totalBytes());
        }
        return *this;
    }
    
    void * blob () const             { return dataPtr; }                /**< Return an opaque pointer to the blob */
    int datatype () const            { return _datatype; }              /**< Return stored datatype code */
    size_t length () const           { return _length; }                /**< Return the number of elements in the data */
    size_t size () const             { return _length; }                /**< Return the number of elements in the data */
    
    /** Return the number of bytes used per element, or zero if the datatype is undefined or the blob is \c NULL */
    size_t bytesPerPixel () const    { return (handler == NULL ? 0 : handler->size()); }
    
    /** Return the total size of the data blob, in bytes */
    size_t totalBytes () const       { return _length * bytesPerPixel(); }
    
    /**
     * Determine whether or not the object is empty
     * @return \c true if the data pointer is \c NULL; \c false otherwise
    **/
    bool isEmpty () const            { return (dataPtr == NULL); }
    
    /**
     * Determine whether the object uses data scaling
     * @return \c true if the slope and intercept parameters are set to nontrivial values;
         \c false otherwise
    **/
    bool isScaled () const           { return (slope != 0.0 && (slope != 1.0 || intercept != 0.0)); }
    
    /**
     * Determine whether the datatype is complex
     * @return Currently \c false, always
    **/
    bool isComplex () const          { return false; }
    
    /**
     * Determine whether the datatype is floating point
     * @return \c true if the data represents 32-bit or 64-bit floating point values; \c false
     *   otherwise
    **/
    bool isFloatingPoint () const    { return (_datatype == DT_FLOAT32 || _datatype == DT_FLOAT64); }
    
    /**
     * Determine whether the datatype is an integer type
     * @return \c true if the data represents integers; \c false otherwise
    **/
    bool isInteger () const          { return nifti_is_inttype(_datatype); }
    
    /**
     * Determine whether the datatype corresponds to an RGB value
     * @return Currently \c false, always
    **/
    bool isRgb () const              { return false; }
    
    /**
     * Return a similar object to the callee, but with the slope and intercept values reset
     * @return A new \c NiftiImageData object, pointing to the same memory as the callee
    **/
    NiftiImageData unscaled () const { return NiftiImageData(dataPtr, _length, _datatype); }
    
    /**
     * Disown the data blob, removing responsibility for freeing it upon destruction
     * @return A reference to the modified callee
    **/
    NiftiImageData & disown ()       { this->owner = false; return *this; }
    
    /** Obtain a constant iterator corresponding to the start of the blob */
    const Iterator begin () const { return Iterator(*this); }
    
    /** Obtain a constant iterator corresponding to the end of the blob */
    const Iterator end () const { return Iterator(*this, static_cast<char*>(dataPtr) + totalBytes()); }
    
    /** Obtain a mutable iterator corresponding to the start of the blob */
    Iterator begin () { return Iterator(*this); }
    
    /** Obtain a mutable iterator corresponding to the end of the blob */
    Iterator end () { return Iterator(*this, static_cast<char*>(dataPtr) + totalBytes()); }
    
    /**
     * Indexing operator, returning a constant element
     * @param i Index value, where the first dimension moves fastest
     * @return Constant element proxy type
    **/
    const Element operator[] (const size_t i) const { return Element(*this, static_cast<char*>(dataPtr) + (i * bytesPerPixel())); }
    
    /**
     * Indexing operator, returning a mutable element
     * @param i Index value, where the first dimension moves fastest
     * @return Mutable element proxy type
    **/
    Element operator[] (const size_t i) { return Element(*this, static_cast<char*>(dataPtr) + (i * bytesPerPixel())); }
    
    /**
     * Calculate the minimum and maximum values in the blob, as doubles
     * @param min Pointer to the minimum value (output parameter). Will be set to zero if the
     *   datatype is unknown or the data is empty
     * @param max Pointer to the maximum value (output parameter). Will be set to zero if the
     *   datatype is unknown or the data is empty
    **/
    void minmax (double *min, double *max) const
    {
        if (handler == NULL)
        {
            *min = 0.0;
            *max = 0.0;
        }
        else
            handler->minmax(dataPtr, _length, min, max);
    }
};


#ifdef USING_R
template <>
inline bool NiftiImageData::ConcreteTypeHandler<int>::hasNaN () const { return true; }
#endif


/**
 * Thin wrapper around a C-style \c nifti_image struct that allows C++-style destruction. Reference
 * counting is used to allow multiple \c NiftiImage objects to wrap the same \c nifti_image
 * pointer, akin to a \c std::shared_ptr (but without requiring C++11).
 * @author Jon Clayden (<code@clayden.org>)
**/
class NiftiImage
{
public:
    /**
     * Inner class referring to a subset of an image. Currently must refer to the last
     * dimension in the image, i.e., a volume in a 4D parent image, or a slice in a 3D image
    **/
    struct Block
    {
        const NiftiImage &image;        /**< The parent image */
        const int dimension;            /**< The dimension along which the block applies (which should be the last) */
        const int index;                /**< The location along \c dimension */
        
        /**
         * Standard constructor for this class
         * @param image The parent image
         * @param dimension The dimension along which the block applies (which should be the last)
         * @param index The location along \c dimension
         * @exception runtime_error If \c dimension is not the last dimension in the image
        **/
        Block (const NiftiImage &image, const int dimension, const int index)
            : image(image), dimension(dimension), index(index)
        {
            if (dimension != image->ndim)
                throw std::runtime_error("Blocks must be along the last dimension in the image");
        }
        
        /**
         * Copy assignment operator, which allows a block in one image to be replaced with
         * the contents of another image
         * @param source A \ref NiftiImage, containing the data to replace the block with
         * @return A reference to the block
         * @exception runtime_error If the \c source is incompatible with the block in size or
         * datatype
        **/
        Block & operator= (const NiftiImage &source)
        {
            if (source->datatype != image->datatype)
                throw std::runtime_error("New data does not have the same datatype as the target block");
            if (source->scl_slope != image->scl_slope || source->scl_inter != image->scl_inter)
                throw std::runtime_error("New data does not have the same scale parameters as the target block");
            
            size_t blockSize = 1;
            for (int i=1; i<dimension; i++)
                blockSize *= image->dim[i];
            
            if (blockSize != source->nvox)
                throw std::runtime_error("New data does not have the same size as the target block");
            
            blockSize *= image->nbyper;
            memcpy(static_cast<char*>(image->data) + blockSize*index, source->data, blockSize);
            return *this;
        }
        
        /**
         * Obtain the data within the block
         * @return A \c NiftiImageData object encapsulating the data
        **/
        NiftiImageData data () const
        {
            if (image.isNull())
                return NiftiImageData();
            else
            {
                size_t blockSize = 1;
                for (int i=1; i<dimension; i++)
                    blockSize *= image->dim[i];
                return NiftiImageData(static_cast<char*>(image->data) + blockSize * index * image->nbyper, blockSize, image->datatype, static_cast<double>(image->scl_slope), static_cast<double>(image->scl_inter));
            }
        }
        
        /**
         * Extract a vector of data from a block, casting it to any required element type
         * @param useSlope If \c true, the default, then the data will be adjusted for the slope
         * and intercept stored with the image, if any
         * @note If the slope and intercept are applied, there is no guarantee that the adjusted
         * values will fit within the requested type. No check is made for this
        **/
        template <typename TargetType>
        std::vector<TargetType> getData (const bool useSlope = true) const;
    };
    
#ifdef USING_R
    /**
     * Convert between R \c SEXP object type and \c nifti_image datatype codes
     * @param sexpType A numeric R \c SEXP type code
     * @return A \c nifti_image datatype code
     * @exception runtime_error If a non-numeric type is passed
    **/
    static int sexpTypeToNiftiType (const int sexpType)
    {
        if (sexpType == INTSXP || sexpType == LGLSXP)
            return DT_INT32;
        else if (sexpType == REALSXP)
            return DT_FLOAT64;
        else
            throw std::runtime_error("Array elements must be numeric");
    }
#endif
    
    /**
     * Extract the pure rotation part of a 4x4 xform matrix
     * @param matrix An xform matrix
     * @return A 3x3 rotation matrix
    **/
    static mat33 xformToRotation (const mat44 matrix);
    
    /**
     * Convert a 4x4 xform matrix to a string describing its canonical axes
     * @param matrix An xform matrix
     * @return A string containing three characters
    **/
    static std::string xformToString (const mat44 matrix);
    
    /**
     * Get the NIfTI format version used by the file at the specified path
     * @param path A string specifying a file path
     * @return An integer: -1 if the file is not present or not valid, 0 for ANALYZE-7.5, or
     *         a value greater than 0 for NIfTI
    **/
    static int fileVersion (const std::string &path);
    

protected:
    nifti_image *image;         /**< The wrapped \c nifti_image pointer */
    int *refCount;              /**< A reference counter, shared with other objects wrapping the same pointer */
    
    /**
     * Acquire the specified pointer to a \c nifti_image \c struct, taking (possibly shared)
     * responsibility for freeing the associated memory. If the object currently wraps another
     * pointer, it will be released
     * @param image The pointer to wrap
    **/
    void acquire (nifti_image * const image);
    
    /**
     * Acquire the same pointer as another \c NiftiImage, incrementing the shared reference count
     * @param source A reference to a \c NiftiImage
    **/
    void acquire (const NiftiImage &source)
    {
        refCount = source.refCount;
        acquire(source.image);
    }
    
    /**
     * Release the currently wrapped pointer, if it is not \c NULL, decrementing the reference
     * count and releasing memory if there are no remaining references to the pointer
    **/
    void release ();
    
    /**
     * Copy the contents of a \c nifti_image to create a new image, acquiring the new pointer
     * @param source A pointer to a \c nifti_image
    **/
    void copy (const nifti_image *source);
    
    /**
     * Copy the contents of another \c NiftiImage to create a new image, acquiring a new pointer
     * @param source A reference to a \c NiftiImage
    **/
    void copy (const NiftiImage &source);
    
    /**
     * Copy the contents of a \ref Block to create a new image, acquiring a new pointer
     * @param source A reference to a \ref Block
    **/
    void copy (const Block &source);


#ifdef USING_R

    /**
     * Initialise the object from an S4 object of class \c "nifti"
     * @param object The source object
     * @param copyData If \c true, the data are copied in; otherwise just the metadata is extracted
    **/
    void initFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);
    
    /**
     * Initialise the object from a reference object of class \c "MriImage"
     * @param object The source object
     * @param copyData If \c true, the data are copied in; otherwise just the metadata is extracted
    **/
    void initFromMriImage (const Rcpp::RObject &object, const bool copyData = true);
    
    /**
     * Initialise the object from an R list with named elements, which can only contain metadata
     * @param object The source object
    **/
    void initFromList (const Rcpp::RObject &object);
    
    /**
     * Initialise the object from an R array
     * @param object The source object
     * @param copyData If \c true, the data are copied in; otherwise just the metadata is extracted
    **/
    void initFromArray (const Rcpp::RObject &object, const bool copyData = true);
   
#endif

    /**
     * Modify the pixel dimensions, and potentially the xform matrices to match
     * @param pixdim Vector of new pixel dimensions
    **/
    void updatePixdim (const std::vector<float> &pixdim);
    
    /**
     * Modify the pixel dimension units
     * @param pixunits Vector of new pixel units, specified using their standard abbreviations
    **/
    void setPixunits (const std::vector<std::string> &pixunits);
    
public:
    /**
     * Default constructor
    **/
    NiftiImage ()
        : image(NULL), refCount(NULL) {}
    
    /**
     * Copy constructor
     * @param source Another \c NiftiImage object
     * @param copy If \c true, the underlying \c nifti_image will be copied; otherwise the new
     * object wraps the same \c nifti_image and increments the shared reference count
    **/
    NiftiImage (const NiftiImage &source, const bool copy = true)
        : image(NULL), refCount(NULL)
    {
        if (copy)
            this->copy(source);
        else
            acquire(source);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from NiftiImage)\n", this->image);
#endif
    }
    
    /**
     * Initialise from a block, copying in the data
     * @param source A \c Block object, referring to part of another \c NiftiImage
    **/
    NiftiImage (const Block &source)
        : image(NULL), refCount(NULL)
    {
        this->copy(source);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from Block)\n", this->image);
#endif
    }
    
    /**
     * Initialise using an existing \c nifti_image pointer
     * @param image An existing \c nifti_image pointer, possibly \c NULL
     * @param copy If \c true, the image data will be copied; otherwise this object just wraps
     * the pointer passed to it
    **/
    NiftiImage (nifti_image * const image, const bool copy = false)
        : image(NULL), refCount(NULL)
    {
        if (copy)
            this->copy(image);
        else
            acquire(image);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from pointer)\n", this->image);
#endif
    }
    
    /**
     * Initialise using a path string
     * @param path A string specifying a path to a valid NIfTI-1 file, possibly gzipped
     * @param readData If \c true, the data will be read as well as the metadata
     * @exception runtime_error If reading from the file fails
    **/
    NiftiImage (const std::string &path, const bool readData = true)
        : image(NULL), refCount(NULL)
    {
        acquire(nifti_image_read(path.c_str(), readData));
        if (image == NULL)
            throw std::runtime_error("Failed to read image from path " + path);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from string)\n", this->image);
#endif
    }
    
    /**
     * Initialise using a path string and sequence of required volumes
     * @param path A string specifying a path to a valid NIfTI-1 file, possibly gzipped
     * @param volumes The volumes to read in (squashing all dimensions above the third together)
     * @exception runtime_error If reading from the file fails, or \c volumes is empty
    **/
    NiftiImage (const std::string &path, const std::vector<int> &volumes);
    
#ifdef USING_R
    /**
     * Initialise from an R object, retrieving an existing image from an external pointer attribute
     * if available; otherwise constructing a new one from the R object itself
     * @param object The source object
     * @param readData If \c true, the data will be retrieved as well as the metadata
     * @param readOnly If \c true, the caller asserts that its intent is read-only. Otherwise, if
     * the \c SEXP may have multiple names at the R level (according to the \c MAYBE_SHARED R
     * macro), an image retrieved from an external pointer will be duplicated to preserve R's usual
     * semantics
    **/
    NiftiImage (const SEXP object, const bool readData = true, const bool readOnly = false);
#endif
    
    /**
     * Destructor which decrements the reference counter, and releases the wrapped pointer if the
     * counter drops to zero
    **/
    virtual ~NiftiImage () { release(); }
    
    /**
     * Allows a \c NiftiImage object to be treated as a pointer to a \c const \c nifti_image
    **/
    operator const nifti_image* () const { return image; }
    
    /**
     * Allows a \c NiftiImage object to be treated as a pointer to a \c nifti_image
    **/
    operator nifti_image* () { return image; }
   
    /**
     * Allows a \c NiftiImage object to be treated as a pointer to a \c const \c nifti_image
    **/
    const nifti_image * operator-> () const { return image; }
    
    /**
     * Allows a \c NiftiImage object to be treated as a pointer to a \c nifti_image
    **/
    nifti_image * operator-> () { return image; }
    
    /**
     * Copy assignment operator, which copies from its argument
     * @param source Another \c NiftiImage
    **/
    NiftiImage & operator= (const NiftiImage &source)
    {
        copy(source);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from NiftiImage)\n", this->image);
#endif
        return *this;
    }
    
    /**
     * Copy assignment operator, which allows a block to be used to replace the contents of a
     * suitably sized image
     * @param source A reference to a suitable \ref Block object
    **/
    NiftiImage & operator= (const Block &source)
    {
        copy(source);
#ifndef NDEBUG
        Rc_printf("Creating NiftiImage with pointer %p (from Block)\n", this->image);
#endif
        return *this;
    }
    
    /**
     * Mark the image as persistent, so that it can be passed back to R
     * @param persistent The new persistence state of the object
     * @return A reference to the callee.
     * @deprecated The persistence mechanism has been replaced with reference counting, so this
     * function no longer has any effect. Instead it returns \c *this, unmodified.
    **/
    NiftiImage & setPersistence (const bool persistent) { return *this; }
    
    /**
     * Determine whether or not the wrapped pointer is \c NULL
     * @return \c true if the wrapped pointer is \c NULL; \c false otherwise
    **/
    bool isNull () const { return (image == NULL); }
    
    /**
     * Determine whether the wrapped pointer is shared with another \c NiftiImage
     * @return \c true if the reference count is greater than 1; \c false otherwise
    **/
    bool isShared () const { return (refCount != NULL && *refCount > 1); }
    
    /**
     * Determine whether or not the image is marked as persistent
     * @return \c false, always
     * @deprecated The persistence mechanism has been replaced with reference counting, so this
     * function will always return \c false. Use \ref isShared instead.
    **/
    bool isPersistent () const { return false; }
    
    /**
     * Determine whether nontrivial scale and slope parameters are set
     * @return \c true if the object wraps an image pointer, its slope is not zero and the slope
     *         and intercept are not exactly one and zero; \c false otherwise
    **/
    bool isDataScaled () const { return (image != NULL && image->scl_slope != 0.0 && (image->scl_slope != 1.0 || image->scl_inter != 0.0)); }
    
    /**
     * Return the number of dimensions in the image
     * @return An integer giving the image dimensionality
    **/
    int nDims () const
    {
        if (image == NULL)
            return 0;
        else
            return image->ndim;
    }
    
    /**
     * Return the dimensions of the image
     * @return A vector of integers giving the width in each dimension
    **/
    std::vector<int> dim () const
    {
        if (image == NULL)
            return std::vector<int>();
        else
            return std::vector<int>(image->dim+1, image->dim+image->ndim+1);
    }
    
    /**
     * Return the dimensions of the pixels or voxels in the image
     * @return A vector of floating-point values giving the pixel width in each dimension
    **/
    std::vector<float> pixdim () const
    {
        if (image == NULL)
            return std::vector<float>();
        else
            return std::vector<float>(image->pixdim+1, image->pixdim+image->ndim+1);
    }
    
    /**
     * Drop unitary dimensions
     * @return Self, after possibly reducing the dimensionality of the image
     * @note This function differs from its R equivalent in only dropping unitary dimensions after
     * the last nonunitary one
    **/
    NiftiImage & drop ()
    {
        int ndim = image->ndim;
        while (image->dim[ndim] < 2)
            ndim--;
        image->dim[0] = image->ndim = ndim;
        
        return *this;
    }
    
    /**
     * Obtain the pixel data within the image
     * @return A constant \c NiftiImageData object encapsulating the data
    **/
    const NiftiImageData data () const { return NiftiImageData(image); }
    
    /**
     * Obtain the pixel data within the image
     * @return A mutable \c NiftiImageData object encapsulating the data
    **/
    NiftiImageData data () { return NiftiImageData(image); }
    
    /**
     * Extract a vector of data from the image, casting it to any required element type
     * @param useSlope If \c true, the default, then the data will be adjusted for the slope and
     * intercept stored with the image, if any
     * @return A vector of data values, cast to the required type
     * @note If the slope and intercept are applied, there is no guarantee that the adjusted values
     * will fit within the requested type. No check is made for this
     * @deprecated Use of the (ultimately more flexible) \ref data methods is now preferred
    **/
    template <typename TargetType>
    std::vector<TargetType> getData (const bool useSlope = true) const;
    
    /**
     * Change the datatype of the image, casting the pixel data if present
     * @param datatype A NIfTI datatype code
     * @param useSlope If \c true, and conversion is to an integer type, the data will be rescaled
     * and the image's slope and intercept set to capture the full range of original values
     * @return Self, after changing the datatype
    **/
    NiftiImage & changeDatatype (const int datatype, const bool useSlope = false);
    
    /**
     * Change the datatype of the image, casting the pixel data if present
     * @param datatype A string specifying the new datatype
     * @param useSlope If \c true, and conversion is to an integer type, the data will be rescaled
     * and the image's slope and intercept set to capture the full range of original values
     * @return Self, after changing the datatype
    **/
    NiftiImage & changeDatatype (const std::string &datatype, const bool useSlope = false);
    
    /**
     * Replace the pixel data in the image with the contents of a vector
     * @param data A data vector, whose elements will be used to replace the image data
     * @param datatype The final datatype required. By default the existing datatype of the image
     * is used
     * @exception runtime_error If the length of the new data does not match the image
     * @return Self, after replacing the data
    **/
    template <typename SourceType>
    NiftiImage & replaceData (const std::vector<SourceType> &data, const int datatype = DT_NONE);
    
    /**
     * Replace the pixel data in the image with the contents of a \c NiftiImageData object
     * @param data A data object, whose elements will be case to match the datatype of the image
     * @exception runtime_error If the length of the new data does not match the image
     * @return Self, after replacing the data
    **/
    NiftiImage & replaceData (const NiftiImageData &data);
    
    /**
     * Drop the data from the image, retaining only the metadata
     * @return Self, after dropping the data
    **/
    NiftiImage & dropData ()
    {
        nifti_image_unload(image);
        return *this;
    }
    
    /**
     * Rescale the image, changing its image dimensions and pixel dimensions
     * @param scales Vector of scale factors along each dimension
     * @return Self, after rescaling the metadata
     * @note No interpolation is performed on the pixel data, which is simply dropped
    **/
    NiftiImage & rescale (const std::vector<float> &scales);
    
    /**
     * Reorient the image by permuting dimensions and potentially reversing some
     * @param i,j,k Constants such as \c NIFTI_L2R, \c NIFTI_P2A and \c NIFTI_I2S, giving the
     * canonical axes to reorient to
     * @return Self, after reorientation
     * @note The pixel data is reordered, but not resampled. The xform matrices will also be
     * adjusted in line with the transformation
    **/
    NiftiImage & reorient (const int i, const int j, const int k);
    
    /**
     * Reorient the image by permuting dimensions and potentially reversing some
     * @param orientation A string containing some permutation of the letters \c L or \c R,
     * \c P or \c A, \c I or \c S, giving the canonical axes to reorient to
     * @return Self, after reorientation
     * @note The pixel data is reordered, but not resampled. The xform matrices will also be
     * adjusted in line with the transformation
    **/
    NiftiImage & reorient (const std::string &orientation);
    
#ifdef USING_R
    /**
     * Update the image from an R array
     * @param array An R array or list object
     * @return Self, after updating data and/or metadata
    **/
    NiftiImage & update (const Rcpp::RObject &object);
#endif
    
    /**
     * Obtain an xform matrix, indicating the orientation of the image
     * @param preferQuaternion If \c true, use the qform matrix in preference to the sform
     * @return A 4x4 matrix
    **/
    mat44 xform (const bool preferQuaternion = true) const;
    
    /**
     * Return the number of blocks in the image
     * @return An integer giving the number of blocks in the image
    **/
    int nBlocks () const
    {
        if (image == NULL)
            return 0;
        else
            return image->dim[image->ndim];
    }
    
    /**
     * Extract a block from the image
     * @param i The block number required
     * @return A \ref Block object
     * @note \ref slice and \ref volume are variants of this function specific to 3D and 4D images,
     * respectively, which may be preferred in some cases for clarity
    **/
    const Block block (const int i) const { return Block(*this, nDims(), i); }
    
    /**
     * Extract a block from the image
     * @param i The block number required
     * @return A \ref Block object
     * @note \ref slice and \ref volume are variants of this function specific to 3D and 4D images,
     * respectively, which may be preferred in some cases for clarity
    **/
    Block block (const int i) { return Block(*this, nDims(), i); }
    
    /**
     * Extract a slice block from a 3D image
     * @param i The slice number required
     * @return A \ref Block object
    **/
    const Block slice (const int i) const { return Block(*this, 3, i); }
    
    /**
     * Extract a slice block from a 3D image
     * @param i The slice number required
     * @return A \ref Block object
    **/
    Block slice (const int i) { return Block(*this, 3, i); }
    
    /**
     * Extract a volume block from a 4D image
     * @param i The volume number required
     * @return A \ref Block object
    **/
    const Block volume (const int i) const { return Block(*this, 4, i); }
    
    /**
     * Extract a volume block from a 4D image
     * @param i The volume number required
     * @return A \ref Block object
    **/
    Block volume (const int i) { return Block(*this, 4, i); }
    
    /**
     * Write the image to a NIfTI-1 file
     * @param fileName The file name to write to, with appropriate suffix (e.g. ".nii.gz")
     * @param datatype The datatype to use when writing the file
    **/
    void toFile (const std::string fileName, const int datatype = DT_NONE) const;
    
    /**
     * Write the image to a NIfTI-1 file
     * @param fileName The file name to write to, with appropriate suffix (e.g. ".nii.gz")
     * @param datatype The datatype to use when writing the file, or "auto"
    **/
    void toFile (const std::string fileName, const std::string &datatype) const;
    
#ifdef USING_R
    
    /**
     * Create an R array from the image
     * @return A numeric array object with an external pointer attribute
    **/
    Rcpp::RObject toArray () const;
    
    /**
     * Create an internal image to pass back to R
     * @param label A string labelling the image
     * @return An R character string with additional attributes
    **/
    Rcpp::RObject toPointer (const std::string label) const;
    
    /**
     * A conditional method that calls either \ref toArray or \ref toPointer
     * @param internal If \c true, \ref toPointer will be called; otherwise \ref toArray
     * @param label A string labelling the image
     * @return An R object
    **/
    Rcpp::RObject toArrayOrPointer (const bool internal, const std::string label) const;
    
    /**
     * Create an R list containing raw image metadata
     * @return An R list
    **/
    Rcpp::RObject headerToList () const;
    
#endif

};

// Include implementations
#include "lib/NiftiImage_impl.h"

} // main namespace

#endif
