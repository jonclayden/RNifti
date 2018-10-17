#' Internal images
#' 
#' An internal image is a simple R object with a few attributes including a
#' pointer to an internal C structure, which contains the full image data. They
#' are used in the package for efficiency, but can be converted to a normal
#' R array using the \code{as.array} method. Attributes of these objects should
#' not be changed.
#' 
#' @param x An \code{"internalImage"} object.
#' @param value Not used. Changing the dimensions of an internal image is
#'   invalid, and will produce an error.
#' @param ... Additional parameters to methods. Currently unused.
#' 
#' @author Jon Clayden <code@@clayden.org>
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

#' @export
print.niftiImage <- function (x, ...)
{
    dim <- dim(x)
    ndim <- length(dim)
    pixdim <- attr(x, "pixdim")
    pixunits <- attr(x, "pixunits")
    
    if ("internalImage" %in% class(x))
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
#' fname = system.file("extdata", "example.nii.gz", package="RNifti")
#' im <- readNifti(fname)
#' pixdim(im)
#' pixdim(fname)
#' pixunits(im)
#' pixunits(fname) # should be unknown
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
    else if (!is.null(dim(object)))
        return (rep(1, length(dim(object))))
    else {
        pdim = tryCatch({
            hdr = niftiHeader(object)
            d1 = abs(hdr$dim[1])
            hdr$pixdim[seq(2, 2 + d1 - 1)]
        })
        if (!inherits(pdim, "try-error")) {
            return(pdim)
        }
        return (1)
        
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
    if ("internalImage" %in% class(object))
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
        return ("Unknown")
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
    if ("internalImage" %in% class(object))
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
