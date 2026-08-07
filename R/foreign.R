#' @rdname asNifti
#' @order 3
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

#' @rdname asNifti
#' @order 4
#' @export
asNifti.MriImage <- function (x, datatype = "auto", internal = NA, ...)
{
    data <- x$data
    if (inherits(data, "SparseArray"))
        data <- as.array(data)
    
    # MriImage tags may include other NIfTI fields, but also unrelated metadata that becomes attributes
    header <- niftiHeader(x$tags)
    attribNames <- setdiff(names(x$tags), names(header))
    
    ndim <- length(x$imageDims)
    
    # This method is available since TractoR v3.0.0, and returns an implicit xform if none is stored
    xform <- structure(x$getXform(), code=2L)
    sform(header) <- xform
    qform(header) <- xform
    
    if (is.null(data))
    {
        qfac <- determinant(xform)$sign
        header$dim <- as.integer(c(ndim, x$imageDims, rep(1L,7-ndim)))
        header$pixdim <- c(qfac, abs(x$voxelDims), rep(0,7-ndim))
        header$xyzt_units <- sum(c(unknown=0L,m=1L,mm=2L,um=3L,s=8L,ms=16L,us=24L)[x$voxelDimUnits])
        
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
        
        # The header is used first below, so its own pixdim (still at niftiHeader()'s arbitrary
        # default here) needs to already agree with the array's, or asNifti() will interpret the
        # "change" as a request to rescale the xform we just set, rather than a benign no-op
        header$pixdim[1 + seq_len(ndim)] <- abs(x$voxelDims)
        
        result <- asNifti(data, header, datatype=datatype, internal=internal, ...)
    }
    
    if (length(attribNames) > 0L)
        imageAttributes(result) <- x$tags[attribNames]
    
    return (result)
}
