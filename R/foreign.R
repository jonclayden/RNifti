#' asNifti() methods for foreign image classes
#'
#' These \code{\link{asNifti}} methods extend it to classes from other packages: \code{"nifti"}
#' objects from \code{oro.nifti}, and \code{"MriImage"} objects from \code{tractor.base}. They
#' are defined here in exactly the way that any third-party package could extend \code{asNifti}
#' for its own classes, as a proof of concept for that general mechanism. The conversion work
#' itself is done in R, where extracting S4 slots and reference-class fields is a more natural
#' fit than the equivalent C++; only genuinely delicate numerical work (the xform matrix, for
#' \code{"MriImage"}) is left to the existing, already-exposed C++ machinery, via the ordinary
#' \code{srow_x}/\code{srow_y}/\code{srow_z} header fields rather than a fresh mutation step.
#' Historically, this logic was instead hardcoded into the underlying C++ constructor, from
#' before \code{asNifti} became an S3 generic.
#'
#' Unlike \code{\link{asNifti}} itself, these methods do not accept a \code{reference} argument:
#' the source object is already a complete description of the image, so there is no external
#' information to merge in, and folding one in would mean a second, redundant pass through
#' \code{asNifti} to apply it.
#'
#' @param x An object of the corresponding class.
#' @param datatype,internal,... See \code{\link{asNifti}}.
#' @return An array or internal image, with class \code{"niftiImage"}.
#'
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{asNifti}}
#' @name asNifti-foreign
#' @rdname asNifti-foreign
NULL

#' @rdname asNifti-foreign
#' @export
asNifti.nifti <- function (x, datatype = "auto", internal = NA, ...)
{
    # oro.nifti upscales data to one of these two types on read, and applies its own scaling,
    # so the original on-disk datatype and slope/intercept are not relevant here
    if (!(x@datatype %in% c(2L,4L,8L,256L,512L,768L,16L,64L)))
        stop("Data type is not supported")
    data <- x@.Data
    storage.mode(data) <- if (x@datatype %in% c(16L,64L)) "double" else "integer"
    
    # The S4 slots use NIfTI-native field names almost throughout, so can be extracted wholesale
    # rather than one at a time, aside from a couple of small adjustments. "dim_info" is dropped
    # because oro.nifti stores it as an unused empty string, not the numeric value expected here
    header <- attributes(x)
    header$class <- NULL
    header$dim <- header$dim_
    header$dim_info <- NULL
    header$scl_slope <- 0
    header$scl_inter <- 0
    
    asNifti(data, header, datatype=datatype, internal=internal, ...)
}

#' @rdname asNifti-foreign
#' @export
asNifti.MriImage <- function (x, datatype = "auto", internal = NA, ...)
{
    data <- x$data
    if (inherits(data, "SparseArray"))
        data <- as.array(data)
    
    header <- niftiHeader(x$tags)
    attribNames <- setdiff(names(x$tags), names(header))
    
    ndim <- length(x$imageDims)
    
    # This method is available since TractoR v3.0.0, and returns an implicit xform if none is stored
    xform <- structure(x$getXform(), code=2L)
    
    if (is.null(data))
    {
        # dim and pixdim should be set before assigning an xform
        qfac <- determinant(xform)$sign
        header$dim <- as.integer(c(ndim, x$imageDims, rep(1L,7-ndim)))
        header$pixdim <- c(qfac, abs(x$voxelDims), rep(0,7-ndim))
        header$xyzt_units <- sum(c(unknown=0L,m=1L,mm=2L,um=3L,s=8L,ms=16L,us=24L)[x$voxelDimUnits])
        
        sform(header) <- xform
        qform(header) <- xform
        
        # Creating an internal image makes sense when there's no data
        result <- asNifti(header, datatype=datatype, internal=TRUE, ...)
    }
    else
    {
        # Setting this metadata here means it will be properly interpreted by the C++ code in one
        # shot, without having to finagle it into NIfTI's required form
        dim(data) <- x$imageDims
        attr(data, "pixdim") <- abs(x$voxelDims)
        attr(data, "pixunits") <- x$voxelDimUnits
        
        sform(header) <- xform
        qform(header) <- xform
        
        result <- asNifti(data, header, datatype=datatype, internal=internal, ...)
    }
    
    if (length(attribNames) > 0L)
        imageAttributes(result) <- x$tags[attribNames]
    
    return (result)
}
