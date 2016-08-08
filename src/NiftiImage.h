#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

#include "nifti1_io.h"

// Thin wrapper around a C-style nifti_image struct that allows C++-style destruction
class NiftiImage
{
public:
    struct Block
    {
        const NiftiImage &image;
        const int dimension;
        const int index;
        
        Block (const NiftiImage &image, const int dimension, const int index)
            : image(image), dimension(dimension), index(index)
        {
            if (dimension != image->ndim)
                throw std::runtime_error("Blocks must be along the last dimension in the image");
        }
        
        Block & operator= (const NiftiImage &source)
        {
            if (source->datatype != image->datatype)
                throw std::runtime_error("New data does not have the same datatype as the target block");
            
            size_t blockSize = 1;
            for (int i=1; i<dimension; i++)
                blockSize *= image->dim[i];
            
            if (blockSize != source->nvox)
                throw std::runtime_error("New data does not have the same size as the target block");
            
            blockSize *= image->nbyper;
            memcpy(static_cast<char*>(image->data) + blockSize*index, source->data, blockSize);
            return *this;
        }
    };
    
    static std::map<std::string,short> DatatypeCodes;
    static std::map<std::string,short> buildDatatypeCodes ()
    {
        std::map<std::string,short> codes;
        codes["auto"] = DT_NONE;
        codes["none"] = DT_NONE;
        codes["unknown"] = DT_NONE;
        codes["uint8"] = DT_UINT8;
        codes["char"] = DT_UINT8;
        codes["int16"] = DT_INT16;
        codes["short"] = DT_INT16;
        codes["int32"] = DT_INT32;
        codes["int"] = DT_INT32;
        codes["float32"] = DT_FLOAT32;
        codes["float"] = DT_FLOAT32;
        codes["float64"] = DT_FLOAT64;
        codes["double"] = DT_FLOAT64;
        codes["int8"] = DT_INT8;
        codes["uint16"] = DT_UINT16;
        codes["uint32"] = DT_UINT32;
        codes["int64"] = DT_INT64;
        codes["uint64"] = DT_UINT64;
        return codes;
    }
    
    static short sexpTypeToNiftiType (const int sexpType)
    {
        if (sexpType == INTSXP || sexpType == LGLSXP)
            return DT_INT32;
        else if (sexpType == REALSXP)
            return DT_FLOAT64;
        else
            throw std::runtime_error("Array elements must be numeric");
    }
    
protected:
    nifti_image *image;
    bool persistent;
    
    void copy (nifti_image * const source);
    void copy (const NiftiImage &source);
    void copy (const Block &source);
    
    void initFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);
    void initFromMriImage (const Rcpp::RObject &object, const bool copyData = true);
    void initFromList (const Rcpp::RObject &object);
    void initFromArray (const Rcpp::RObject &object, const bool copyData = true);
    
    void updatePixdim (const std::vector<float> &pixdim);
    void setPixunits (const std::vector<std::string> &pixunits);
    
public:
    NiftiImage ()
        : image(NULL), persistent(false) {}
    
    NiftiImage (const NiftiImage &source)
        : persistent(false)
    {
        this->copy(source);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p\n", this->image);
#endif
    }
    
    NiftiImage (nifti_image * const image, const bool copy = false)
        : image(image), persistent(false)
    {
        if (copy)
            this->copy(image);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p\n", this->image);
#endif
    }
    
    NiftiImage (const SEXP object, const bool readData = true);
    
    ~NiftiImage ()
    {
        if (!persistent)
        {
#ifndef NDEBUG
            Rprintf("Freeing NiftiImage with pointer %p\n", this->image);
#endif
            nifti_image_free(image);
        }
    }
    
    operator nifti_image* () const { return image; }
    
    nifti_image * operator-> () const { return image; }
    
    NiftiImage & operator= (const NiftiImage &source)
    {
        copy(source);
        return *this;
    }
    
    NiftiImage & operator= (const Block &source)
    {
        copy(source);
        return *this;
    }
    
    void setPersistence (const bool persistent)
    {
        this->persistent = persistent;
#ifndef NDEBUG
        if (persistent)
            Rprintf("Setting NiftiImage with pointer %p to be persistent\n", this->image);
#endif
    }
    
    bool isNull () const { return (image == NULL); }
    bool isPersistent () const { return persistent; }
    int nDims () const
    {
        if (image == NULL)
            return 0;
        else
            return image->ndim;
    }
    
    // Note that this function differs from its R equivalent in only dropping unitary dimensions after the last nonunitary one
    NiftiImage & drop ()
    {
        int ndim = image->ndim;
        while (image->dim[ndim] < 2)
            ndim--;
        image->dim[0] = image->ndim = ndim;
        
        return *this;
    }
    
    void rescale (const std::vector<float> &scales);
    void update (const SEXP array);
    
    mat44 xform (const bool preferQuaternion = true) const;
    
    const Block slice (const int i) const { return Block(*this, 3, i); }
    const Block volume (const int i) const { return Block(*this, 4, i); }
    
    Block slice (const int i) { return Block(*this, 3, i); }
    Block volume (const int i) { return Block(*this, 4, i); }
    
    void toFile (const std::string fileName, const short datatype) const;
    
    Rcpp::RObject toArray () const;
    Rcpp::RObject toPointer (const std::string label) const;
    Rcpp::RObject toArrayOrPointer (const bool internal, const std::string label) const;
    Rcpp::RObject headerToList () const;
};

NiftiImage allocateMultiregResult (const NiftiImage &source, const NiftiImage &target, const bool forceDouble);

#endif
