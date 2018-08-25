#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_


#ifdef USING_R

#include <Rcpp.h>

#else

#define R_NegInf -INFINITY
#define ISNAN(x) (x != x)

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
 * pages is to document the \ref RNifti::NiftiImage C++ class for package developers linking to
 * \c RNifti.
**/

namespace RNifti {

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

} // internal namespace

/**
 * Thin wrapper around a C-style \c nifti_image struct that allows C++-style destruction
 * @author Jon Clayden
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
    static short sexpTypeToNiftiType (const int sexpType)
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
     * Convert a 4x4 xform matrix to a string describing its canonical axes
     * @param matrix An xform matrix
     * @return A string containing three characters
    **/
    static std::string xformToString (const mat44 matrix)
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
    
    /**
     * Get the NIfTI format version used by the file at the specified path
     * @param path A string specifying a file path
     * @return An integer: -1 if the file is not present or not valid, 0 for ANALYZE-7.5, or
     *         a value greater than 0 for NIfTI
    **/
    static int fileVersion (const std::string &path)
    {
        nifti_1_header *header = nifti_read_header(path.c_str(), NULL, false);
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
            return version;
        }
    }
    

protected:
    nifti_image *image;         /**< The wrapped \c nifti_image pointer */
    bool persistent;            /**< Marker of persistence, which determines whether the nifti_image should be freed on destruction */
    
    /**
     * Copy the contents of a \c nifti_image to create a new image
     * @param source A pointer to a \c nifti_image
    **/
    void copy (const nifti_image *source);
    
    /**
     * Copy the contents of another \c NiftiImage to create a new image
     * @param source A reference to a \c NiftiImage
    **/
    void copy (const NiftiImage &source);
    
    /**
     * Copy the contents of a \ref Block to create a new image
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
        : image(NULL), persistent(false) {}
    
    /**
     * Copy constructor
     * @param source Another \c NiftiImage object
    **/
    NiftiImage (const NiftiImage &source)
        : image(NULL), persistent(false)
    {
        this->copy(source);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p (from NiftiImage)\n", this->image);
#endif
    }
    
    /**
     * Initialise from a block, copying in the data
     * @param source A \c Block object, referring to part of another \c NiftiImage
    **/
    NiftiImage (const Block &source)
        : image(NULL), persistent(false)
    {
        this->copy(source);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p (from Block)\n", this->image);
#endif
    }
    
    /**
     * Initialise using an existing \c nifti_image pointer
     * @param image An existing \c nifti_image pointer, possibly \c NULL
     * @param copy If \c true, the image data will be copied; otherwise this object just wraps
     * the pointer passed to it
    **/
    NiftiImage (nifti_image * const image, const bool copy = false)
        : image(NULL), persistent(false)
    {
        if (copy)
            this->copy(image);
        else
            this->image = image;
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p (from pointer)\n", this->image);
#endif
    }
    
    /**
     * Initialise using a path string
     * @param path A string specifying a path to a valid NIfTI-1 file, possibly gzipped
     * @param readData If \c true, the data will be read as well as the metadata
     * @exception runtime_error If reading from the file fails
    **/
    NiftiImage (const std::string &path, const bool readData = true)
        : persistent(false)
    {
        this->image = nifti_image_read(path.c_str(), readData);
        if (this->image == NULL)
            throw std::runtime_error("Failed to read image from path " + path);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p (from string)\n", this->image);
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
     * Initialise from an R object
     * @param object The source object
     * @param readData If \c true, the data will be copied as well as the metadata
    **/
    NiftiImage (const SEXP object, const bool readData = true);
#endif
    
    /**
     * Destructor which frees the wrapped pointer, unless the object is marked as persistent
    **/
    virtual ~NiftiImage ()
    {
        if (!persistent)
        {
#ifndef NDEBUG
            Rprintf("Freeing NiftiImage with pointer %p\n", this->image);
#endif
            nifti_image_free(image);
        }
    }
    
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
        Rprintf("Creating NiftiImage with pointer %p (from NiftiImage)\n", this->image);
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
        Rprintf("Creating NiftiImage with pointer %p (from Block)\n", this->image);
#endif
        return *this;
    }
    
    /**
     * Marked the image as persistent, so that it can be passed back to R
     * @param persistent The new persistence state of the object
    **/
    NiftiImage & setPersistence (const bool persistent)
    {
        this->persistent = persistent;
#ifndef NDEBUG
        if (persistent)
            Rprintf("Setting NiftiImage with pointer %p to be persistent\n", this->image);
#endif
        return *this;
    }
    
    /**
     * Determine whether or not the internal pointer is \c NULL
    **/
    bool isNull () const { return (image == NULL); }
    
    /**
     * Determine whether or not the image is marked as persistent
    **/
    bool isPersistent () const { return persistent; }
    
    /**
     * Determine whether nontrivial scale and slope parameters are set
    **/
    bool isDataScaled () const { return (image != NULL && image->scl_slope != 0.0 && (image->scl_slope != 1.0 || image->scl_inter != 0.0)); }
    
    /**
     * Return the number of dimensions in the image
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
     * Extract a vector of data from the image, casting it to any required element type
     * @param useSlope If \c true, the default, then the data will be adjusted for the slope and
     * intercept stored with the image, if any
     * @note If the slope and intercept are applied, there is no guarantee that the adjusted values
     * will fit within the requested type. No check is made for this
    **/
    template <typename TargetType>
    std::vector<TargetType> getData (const bool useSlope = true) const;
    
    /**
     * Change the datatype of the image, casting the pixel data if present
     * @param datatype A NIfTI datatype code
     * @param useSlope If \c true, and conversion is to an integer type, the data will be rescaled
     * and the image's slope and intercept set to capture the full range of original values
    **/
    NiftiImage & changeDatatype (const short datatype, const bool useSlope = false);
    
    /**
     * Change the datatype of the image, casting the pixel data if present
     * @param datatype A string specifying the new datatype
     * @param useSlope If \c true, and conversion is to an integer type, the data will be rescaled
     * and the image's slope and intercept set to capture the full range of original values
    **/
    NiftiImage & changeDatatype (const std::string &datatype, const bool useSlope = false);
    
    /**
     * Replace the pixel data in the image with the contents of a vector
     * @param data A data vector, whose elements will be cast to match the datatype of the image.
     * An exception will be raised if this does not have a length matching the image
     * @param datatype The final datatype required. By default the existing datatype of the image
     * is used
    **/
    template <typename SourceType>
    NiftiImage & replaceData (const std::vector<SourceType> &data, const short datatype = DT_NONE);
    
    /**
     * Drop the data from the image, retaining only the metadata
    **/
    NiftiImage & dropData ()
    {
        nifti_image_unload(image);
        return *this;
    }
    
    /**
     * Rescale the image, changing its image dimensions and pixel dimensions
     * @param scales Vector of scale factors along each dimension
     * @note No interpolation is performed on the pixel data, which is simply dropped
    **/
    NiftiImage & rescale (const std::vector<float> &scales);
    
    /**
     * Reorient the image by permuting dimensions and potentially reversing some
     * @param i,j,k Constants such as \c NIFTI_L2R, \c NIFTI_P2A and \c NIFTI_I2S, giving the
     * canonical axes to reorient to
     * @note The pixel data is reordered, but not resampled. The xform matrices will also be
     * adjusted in line with the transformation
    **/
    NiftiImage & reorient (const int i, const int j, const int k);
    
    /**
     * Reorient the image by permuting dimensions and potentially reversing some
     * @param orientation A string containing some permutation of the letters \c L or \c R,
     * \c P or \c A, \c I or \c S, giving the canonical axes to reorient to
     * @note The pixel data is reordered, but not resampled. The xform matrices will also be
     * adjusted in line with the transformation
    **/
    NiftiImage & reorient (const std::string &orientation);
    
#ifdef USING_R
    /**
     * Update the image from an R array
     * @param array An R array or list object
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
    void toFile (const std::string fileName, const short datatype = DT_NONE) const;
    
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

// Include helper functions
#include "lib/NiftiImage_internal.h"

inline void NiftiImage::copy (const nifti_image *source)
{
    if (image != NULL && !persistent)
        nifti_image_free(image);
        
    if (source == NULL)
        image = NULL;
    else
    {
        image = nifti_copy_nim_info(source);
        if (source->data != NULL)
        {
            size_t dataSize = nifti_get_volsize(source);
            image->data = calloc(1, dataSize);
            memcpy(image->data, source->data, dataSize);
        }
    }
    
    persistent = false;
}

inline void NiftiImage::copy (const NiftiImage &source)
{
    const nifti_image *sourceStruct = source;
    copy(sourceStruct);
}

inline void NiftiImage::copy (const Block &source)
{
    if (image != NULL && !persistent)
        nifti_image_free(image);
    
    const nifti_image *sourceStruct = source.image;
    if (sourceStruct == NULL)
        image = NULL;
    else
    {
        image = nifti_copy_nim_info(sourceStruct);
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
    
    persistent = false;
}

#ifdef USING_R

// Convert an S4 "nifti" object, as defined in the oro.nifti package, to a "nifti_image" struct
inline void NiftiImage::initFromNiftiS4 (const Rcpp::RObject &object, const bool copyData)
{
    nifti_1_header header;
    header.sizeof_hdr = 348;
    
    const std::vector<short> dims = object.slot("dim_");
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
    
    header.scl_slope = object.slot("scl_slope");
    header.scl_inter = object.slot("scl_inter");
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
    
    if (header.datatype == DT_UINT8 || header.datatype == DT_INT16 || header.datatype == DT_INT32 || header.datatype == DT_INT8 || header.datatype == DT_UINT16 || header.datatype == DT_UINT32)
        header.datatype = DT_INT32;
    else if (header.datatype == DT_FLOAT32 || header.datatype == DT_FLOAT64)
        header.datatype = DT_FLOAT64;  // This assumes that sizeof(double) == 8
    else
        throw std::runtime_error("Data type is not supported");
    
    this->image = nifti_convert_nhdr2nim(header, NULL);
    
    const SEXP data = PROTECT(object.slot(".Data"));
    if (!copyData || Rf_length(data) <= 1)
        this->image->data = NULL;
    else
    {
        const size_t dataSize = nifti_get_volsize(this->image);
        this->image->data = calloc(1, dataSize);
        if (header.datatype == DT_INT32)
        {
            Rcpp::IntegerVector intData(data);
            std::copy(intData.begin(), intData.end(), static_cast<int32_t*>(this->image->data));
        }
        else
        {
            Rcpp::DoubleVector doubleData(data);
            std::copy(doubleData.begin(), doubleData.end(), static_cast<double*>(this->image->data));
        }
    }
    UNPROTECT(1);
}

inline void NiftiImage::initFromMriImage (const Rcpp::RObject &object, const bool copyData)
{
    Rcpp::Reference mriImage(object);
    Rcpp::Function getXform = mriImage.field("getXform");
    Rcpp::NumericMatrix xform = getXform();
    
    this->image = NULL;
    
    if (Rf_length(mriImage.field("tags")) > 0)
        initFromList(mriImage.field("tags"));
    
    Rcpp::RObject data = mriImage.field("data");
    if (data.inherits("SparseArray"))
    {
        Rcpp::Language call("as.array", data);
        data = call.eval();
    }
    
    const short datatype = (Rf_isNull(data) ? DT_INT32 : sexpTypeToNiftiType(data.sexp_type()));
    
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
        this->image = nifti_make_new_nim(dims, datatype, FALSE);
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
    
    this->image = nifti_convert_nhdr2nim(*header, NULL);
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
    
    const short datatype = sexpTypeToNiftiType(object.sexp_type());
    this->image = nifti_make_new_nim(dims, datatype, int(copyData));
    
    if (copyData)
    {
        const size_t dataSize = nifti_get_volsize(image);
        if (datatype == DT_INT32)
            memcpy(this->image->data, INTEGER(object), dataSize);
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

inline NiftiImage::NiftiImage (const SEXP object, const bool readData)
    : persistent(false)
{
    Rcpp::RObject imageObject(object);
    bool resolved = false;
    
    if (imageObject.hasAttribute(".nifti_image_ptr"))
    {
        Rcpp::XPtr<NiftiImage> imagePtr(SEXP(imageObject.attr(".nifti_image_ptr")));
        NiftiImage *ptr = imagePtr;
        if (ptr != NULL)
        {
            this->image = ptr->image;
            this->persistent = true;
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
            this->image = NULL;
        else if (Rf_isString(object))
        {
            const std::string path = Rcpp::as<std::string>(object);
            this->image = nifti_image_read(path.c_str(), readData);
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
        else
            throw std::runtime_error("Cannot convert object of class \"" + Rcpp::as<std::string>(imageObject.attr("class")) + "\" to a nifti_image");
    }
    
    if (this->image != NULL)
        nifti_update_dims_from_array(this->image);
    
#ifndef NDEBUG
    Rprintf("Creating NiftiImage with pointer %p (from SEXP)\n", this->image);
#endif
}

#endif // USING_R

inline NiftiImage::NiftiImage (const std::string &path, const std::vector<int> &volumes)
    : persistent(false)
{
    if (volumes.empty())
        throw std::runtime_error("The vector of volumes is empty");
    
    nifti_brick_list brickList;
    this->image = nifti_image_read_bricks(path.c_str(), volumes.size(), &volumes[0], &brickList);
    if (this->image == NULL)
        throw std::runtime_error("Failed to read image from path " + path);
    
    size_t brickSize = image->nbyper * image->nx * image->ny * image->nz;
    image->data = calloc(1, nifti_get_volsize(image));
    for (int i=0; i<brickList.nbricks; i++)
        memcpy((char *) image->data + i * brickSize, brickList.bricks[i], brickSize);
    nifti_free_NBL(&brickList);
    
#ifndef NDEBUG
    Rprintf("Creating NiftiImage with pointer %p (from string and volume vector)\n", this->image);
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
                else if (i >= nDims)
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
    if (!persistent)
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
        
        const std::vector<double> oldData = this->getData<double>();
        std::vector<double> newData(image->nvox);
        
        // Where the sign is negative we need to start at the end of the dimension
        size_t volStart = 0;
        for (int i=0; i<3; i++)
        {
            if (revsigns[i] < 0)
                volStart += (image->dim[i+1] - 1) * strides[i];
        }
        
        // Iterate over the data and place it into a new vector
        std::vector<double>::const_iterator it = oldData.begin();
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
        switch(orientation[i])
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
            nifti_image *newImage = nifti_convert_nhdr2nim(*header, NULL);
            if (this->image->data != NULL)
            {
                size_t dataSize = nifti_get_volsize(image);
                newImage->data = calloc(1, dataSize);
                memcpy(newImage->data, image->data, dataSize);
            }
            if (!persistent)
                nifti_image_free(image);
            this->image = newImage;
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
        nifti_datatype_sizes(image->datatype, &image->nbyper, NULL);
    
        if (!persistent)
            nifti_image_unload(image);
    
        const size_t dataSize = nifti_get_volsize(image);
        image->data = calloc(1, dataSize);
        if (image->datatype == DT_INT32)
            memcpy(image->data, INTEGER(object), dataSize);
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
    if (image.isNull())
        return std::vector<TargetType>();
    
    size_t blockSize = 1;
    for (int i=1; i<dimension; i++)
        blockSize *= image->dim[i];

    std::vector<TargetType> data(blockSize);
    internal::DataConverter<TargetType> *converter = NULL;
    if (useSlope && image.isDataScaled())
        converter = new internal::DataConverter<TargetType>(image->scl_slope, image->scl_inter);
    
    internal::convertData<TargetType>(image->data, image->datatype, blockSize, data.begin(), blockSize*index, converter);
    
    delete converter;
    return data;
}

template <typename TargetType>
inline std::vector<TargetType> NiftiImage::getData (const bool useSlope) const
{
    if (this->isNull())
        return std::vector<TargetType>();
    
    std::vector<TargetType> data(image->nvox);
    internal::DataConverter<TargetType> *converter = NULL;
    if (useSlope && this->isDataScaled())
        converter = new internal::DataConverter<TargetType>(image->scl_slope, image->scl_inter);
    
    internal::convertData<TargetType>(image->data, image->datatype, image->nvox, data.begin(), 0, converter);
    
    delete converter;
    return data;
}

inline NiftiImage & NiftiImage::changeDatatype (const short datatype, const bool useSlope)
{
    if (this->isNull() || image->datatype == datatype)
        return *this;
    
    if (useSlope && this->isDataScaled())
        throw std::runtime_error("Resetting the slope and intercept for an image with them already set is not supported");
    
    if (image->data != NULL)
    {
        int bytesPerPixel;
        nifti_datatype_sizes(datatype, &bytesPerPixel, NULL);
        void *data = calloc(image->nvox, bytesPerPixel);
        
        switch (datatype)
        {
            case DT_UINT8:
            {
                internal::DataConverter<uint8_t> converter(useSlope ? internal::DataConverter<uint8_t>::IndexMode : internal::DataConverter<uint8_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<uint8_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_INT16:
            {
                internal::DataConverter<int16_t> converter(useSlope ? internal::DataConverter<int16_t>::IndexMode : internal::DataConverter<int16_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<int16_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_INT32:
            {
                internal::DataConverter<int32_t> converter(useSlope ? internal::DataConverter<int32_t>::IndexMode : internal::DataConverter<int32_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<int32_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_FLOAT32:
            {
                internal::DataConverter<float> converter(useSlope ? internal::DataConverter<float>::IndexMode : internal::DataConverter<float>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<float *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_FLOAT64:
            {
                internal::DataConverter<double> converter(useSlope ? internal::DataConverter<double>::IndexMode : internal::DataConverter<double>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<double *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_INT8:
            {
                internal::DataConverter<int8_t> converter(useSlope ? internal::DataConverter<int8_t>::IndexMode : internal::DataConverter<int8_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<int8_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_UINT16:
            {
                internal::DataConverter<uint16_t> converter(useSlope ? internal::DataConverter<uint16_t>::IndexMode : internal::DataConverter<uint16_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<uint16_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_UINT32:
            {
                internal::DataConverter<uint32_t> converter(useSlope ? internal::DataConverter<uint32_t>::IndexMode : internal::DataConverter<uint32_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<uint32_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_INT64:
            {
                internal::DataConverter<int64_t> converter(useSlope ? internal::DataConverter<int64_t>::IndexMode : internal::DataConverter<int64_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<int64_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            case DT_UINT64:
            {
                internal::DataConverter<uint64_t> converter(useSlope ? internal::DataConverter<uint64_t>::IndexMode : internal::DataConverter<uint64_t>::CastMode);
                internal::convertData(image->data, image->datatype, image->nvox, static_cast<uint64_t *>(data), 0, &converter);
                image->scl_slope = static_cast<float>(converter.getSlope());
                image->scl_inter = static_cast<float>(converter.getIntercept());
                break;
            }
        
            default:
            throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(datatype)) + ")");
        }
    
        nifti_image_unload(image);
        image->data = data;
    }
    
    image->datatype = datatype;
    nifti_datatype_sizes(datatype, &image->nbyper, &image->swapsize);
    
    return *this;
}

inline NiftiImage & NiftiImage::changeDatatype (const std::string &datatype, const bool useSlope)
{
    return changeDatatype(internal::stringToDatatype(datatype), useSlope);
}

template <typename SourceType>
inline NiftiImage & NiftiImage::replaceData (const std::vector<SourceType> &data, const short datatype)
{
    if (this->isNull())
        return *this;
    else if (data.size() != image->nvox)
        throw std::runtime_error("New data length does not match the number of voxels in the image");
    
    if (datatype != DT_NONE)
    {
        nifti_image_unload(image);
        image->datatype = datatype;
        nifti_datatype_sizes(datatype, &image->nbyper, &image->swapsize);
    }
    
    if (image->data == NULL)
        image->data = calloc(image->nvox, image->nbyper);
    internal::replaceData<SourceType>(data.begin(), data.end(), image->data, image->datatype);
    
    image->scl_slope = 0.0;
    image->scl_inter = 0.0;
    image->cal_min = static_cast<float>(*std::min_element(data.begin(), data.end()));
    image->cal_max = static_cast<float>(*std::max_element(data.begin(), data.end()));
    
    return *this;
}

inline void NiftiImage::toFile (const std::string fileName, const short datatype) const
{
    const bool changingDatatype = (datatype != DT_NONE && !this->isNull() && datatype != image->datatype);
    
    // Copy the source image only if the datatype will be changed
    NiftiImage imageToWrite(image, changingDatatype);
    
    if (changingDatatype)
        imageToWrite.changeDatatype(datatype, true);
    else
        imageToWrite.setPersistence(true);
    
    const int status = nifti_set_filenames(imageToWrite, fileName.c_str(), false, true);
    if (status != 0)
        throw std::runtime_error("Failed to set filenames for NIfTI object");
    nifti_image_write(imageToWrite);
}

inline void NiftiImage::toFile (const std::string fileName, const std::string &datatype) const
{
    toFile(fileName, internal::stringToDatatype(datatype));
}

#ifdef USING_R

inline Rcpp::RObject NiftiImage::toArray () const
{
    Rcpp::RObject array;
    
    if (this->isNull())
        return array;
    else if (this->isDataScaled())
    {
        array = internal::imageDataToArray<REALSXP>(image);
        std::transform(REAL(array), REAL(array)+Rf_length(array), REAL(array), internal::DataConverter<double>(image->scl_slope,image->scl_inter));
    }
    else
    {
        switch (image->datatype)
        {
            case DT_UINT8:
            case DT_INT16:
            case DT_INT32:
            case DT_INT8:
            case DT_UINT16:
            case DT_UINT32:
            case DT_INT64:
            case DT_UINT64:
            array = internal::imageDataToArray<INTSXP>(image);
            break;
        
            case DT_FLOAT32:
            case DT_FLOAT64:
            array = internal::imageDataToArray<REALSXP>(image);
            break;
        
            default:
            throw std::runtime_error("Unsupported data type (" + std::string(nifti_datatype_string(image->datatype)) + ")");
        }
    }
    
    internal::addAttributes(array, image);
    array.attr("class") = Rcpp::CharacterVector::create("niftiImage", "array");
    return array;
}

inline Rcpp::RObject NiftiImage::toPointer (const std::string label) const
{
    if (this->isNull())
        return Rcpp::RObject();
    else
    {
        Rcpp::RObject string = Rcpp::wrap(label);
        internal::addAttributes(string, image, false);
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
    
    result.attr("class") = Rcpp::CharacterVector::create("niftiHeader");
    
    return result;
}

#endif // USING_R

} // namespace

#endif
