#' Internal images
#' 
#' An internal image is a simple R object with a few attributes including a
#' pointer to an internal C structure, which contains the full image data. They
#' are used in the package for efficiency, but can be converted to a normal
#' R array using the \code{as.array} method. Attributes of these objects should
#' not be changed.
#' 
#' @param x An \code{"internalImage"} object.
#' @param value Not used. Changing the dimensions of (or data in) an internal
#'   image is invalid, and will produce an error. Convert to an array first.
#' @param i,j Index vectors. May be missing, which indicates that the whole of
#'   the relevant dimension should be obtained.
#' @param ... Additional parameters to methods. Only used for additional
#'   indices.
#' @param drop If \code{TRUE} (the default), unitary indices in the result will
#'   be dropped. This mirrors the behaviour of standard array indexing.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @aliases internalImage
#' @rdname internalImage
#' @export
dim.internalImage <- function (x)
{
    return (attr(x, "imagedim"))
}

#' @rdname internalImage
#' @export
"dim<-.internalImage" <- function (x, value)
{
    stop("Dimensions of an internal image cannot be changed")
}

#' @rdname internalImage
#' @export
as.array.internalImage <- function (x, ...)
{
    return (.Call("pointerToArray", x, PACKAGE="RNifti"))
}

#' @rdname internalImage
#' @export
"[.internalImage" <- function (x, i, j, ..., drop = TRUE)
{
    nArgs <- nargs() - as.integer(!missing(drop))
    if (nArgs < 2)
        return (as.array(x))
    
    # Evaluate the indices, replacing missing values with -1
    indices <- substitute(list(i,j,...))
    present <- (sapply(indices, as.character)[-1] != "")
    if (any(!present))
        indices[which(!present)+1] <- -1
    indices <- eval(indices, parent.frame())
    lengths <- rep(-1L, nArgs - 1)
    lengths[present] <- sapply(indices[present], length)
    
    dims <- dim(x)
    data <- NULL
    
    if (all(lengths == -1))
        return (as.array(x))
    else if (any(lengths == 0))
        return (numeric(0))
    # else if (nArgs != length(dims) + 1)          # TODO: vector and matrix indexing, which only involve i
        # stop("Number of indices (", nArgs-1, ") not equal to the dimensionality of the image (", length(dims), ")")
    else if (nArgs == length(dims) + 1)
    {
        if (all(lengths %in% c(-1L,1L)))
        {
            data <- .Call("indexCollapsed", x, as.integer(indices), PACKAGE="RNifti")
            dim(data) <- ifelse(present, 1L, dims)
        }
        else if (all(sapply(indices, function(n) all(diff(as.integer(n))==1L))))
        {
            starts <- ifelse(present, sapply(indices,min,na.rm=TRUE), 1L)
            sizes <- ifelse(present, lengths, dims)
            data <- .Call("indexBlock", x, starts, sizes, PACKAGE="RNifti")
            dim(data) <- sizes
        }
    }
    
    if (drop)
        data <- drop(data)
    return (data)
}

#' @rdname internalImage
#' @export
"[<-.internalImage" <- function (x, i, j, ..., value)
{
    stop("The data in an internal image cannot be changed - convert to array first")
}

#' @export
print.niftiImage <- function (x, ...)
{
    dim <- dim(x)
    ndim <- length(dim)
    pixdim <- attr(x, "pixdim")
    pixunits <- attr(x, "pixunits")
    
    if (inherits(x, "internalImage"))
        cat(paste0("Internal image: \"", x, "\"\n"))
    else
        cat(paste0("Image array of mode \"", storage.mode(x), "\" (", format(object.size(x),"auto"), ")\n"))
    
    cat(paste("-", paste(dim,collapse=" x "), ifelse(ndim>2,"voxels\n","pixels\n")))
    
    if (!is.null(pixdim))
    {
        spaceUnit <- grep("m$", pixunits, perl=TRUE, value=TRUE)
        cat(paste("-", paste(signif(pixdim[1:min(3,ndim)],4),collapse=" x ")))
        if (length(spaceUnit) > 0)
            cat(paste0(" ", spaceUnit[1]))
        
        if (ndim > 3)
        {
            timeUnit <- grep("s$", pixunits, perl=TRUE, value=TRUE)
            cat(paste(" x", signif(pixdim[4],4)))
            if (length(timeUnit) > 0)
                cat(paste0(" ", timeUnit[1]))
        }
        if (ndim > 4)
            cat(paste(" x", paste(signif(pixdim[5:ndim],4),collapse=" x ")))
        
        cat(paste(" per", ifelse(ndim>2,"voxel\n","pixel\n")))
    }
}

#' Number of dimensions
#' 
#' This function is shorthand for \code{length(dim(object))}.
#' 
#' @param object An R object.
#' @return The dimensionality of the object. Objects without a \code{dim}
#'   attribute will produce zero.
#' 
#' @examples
#' ndim(array(0L, dim=c(10,10)))
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
ndim <- function (object)
{
    length(dim(object))
}

#' Pixel dimensions and units
#' 
#' By default, these generic functions return or replace the \code{"pixdim"}
#' and \code{"pixunits"} attributes of their arguments. These represent the
#' physical step size between pixel or voxel centre points, and the spatial and
#' temporal units that they are given in. The former defaults to 1 in each
#' dimension, if there is no attribute.
#' 
#' @param object An R object, generally an image.
#' @param value Numeric vector of pixel dimensions along each axis, or
#'   character vector of abbreviated units. For dimensions, a scalar
#'   \code{value} will be recycled if necessary.
#' @return \code{pixdim} returns a numeric vector of pixel dimensions.
#'   \code{pixunits} returns a character vector of length up to two, giving the
#'   spatial and temporal unit names.
#' 
#' @examples
#' im <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' pixdim(im)
#' pixunits(im)
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
pixdim <- function (object)
{
    UseMethod("pixdim")
}

#' @rdname pixdim
#' @export
pixdim.default <- function (object)
{
    if (!is.null(attr(object, "pixdim")))
        return (attr(object, "pixdim"))
    else
    {
        value <- try(attr(niftiHeader(object), "pixdim"), silent=TRUE)
        if (inherits(value, "try-error"))
            return (1)
        else
            return (value)
    }
}

#' @rdname pixdim
#' @export
"pixdim<-" <- function (object, value)
{
    UseMethod("pixdim<-")
}

#' @rdname pixdim
#' @export
"pixdim<-.default" <- function (object, value)
{
    if (inherits(object, "internalImage"))
        stop("Pixel dimensions of an internal image cannot be changed")
    
    if (is.numeric(value))
    {
        if (length(value) == ndim(object))
            attr(object, "pixdim") <- value
        else if (length(value) == 1)
            attr(object, "pixdim") <- rep(value, ndim(object))
    }
    return (object)
}

#' @rdname pixdim
#' @export
pixunits <- function (object)
{
    UseMethod("pixunits")
}

#' @rdname pixdim
#' @export
pixunits.default <- function (object)
{
    if (!is.null(attr(object, "pixunits")))
        return (attr(object, "pixunits"))
    else
    {
        value <- try(attr(niftiHeader(object), "pixunits"), silent=TRUE)
        if (inherits(value, "try-error"))
            return ("Unknown")
        else
            return (value)
    }
}

#' @rdname pixdim
#' @export
"pixunits<-" <- function (object, value)
{
    UseMethod("pixunits<-")
}

#' @rdname pixdim
#' @export
"pixunits<-.default" <- function (object, value)
{
    if (inherits(object, "internalImage"))
        stop("Pixel units of an internal image cannot be changed")
    
    if (is.character(value))
        attr(object, "pixunits") <- value
    
    return (object)
}

#' Access to metadata elements
#' 
#' These methods provide shorthand access to metadata elements from the NIfTI
#' header corresponding to an image. The extraction version returns the
#' corresponding element from the result of \code{niftiHeader}, while the
#' replacement version calls \code{updateNifti} to replace it.
#' 
#' @param x A \code{"niftiImage"} object, internal or otherwise.
#' @param name A string naming the field required.
#' @param value A new value for the field.
#' 
#' @examples
#' im <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' print(im$descrip)
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftiHeader}}, \code{\link{updateNifti}}
#' @rdname indexing
#' @export
"$.niftiImage" <- function (x, name)
{
    return (niftiHeader(x)[[name]])
}

#' @rdname indexing
#' @export
"$<-.niftiImage" <- function (x, name, value)
{
    template <- list()
    template[[name]] <- value
    return (updateNifti(x, template))
}
