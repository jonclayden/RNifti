#' RGB arrays
#' 
#' The \code{rgbArray} function constructs an integer array whose values are
#' byte-packed representations of 8-bit RGBA colour values. The \code{channels}
#' attribute (with value 3 or 4) indicates how many channels are being used.
#' The resulting array can be used to construct an RGB(A) NIfTI image, or
#' converted to standard R colour strings using the \code{as.character} method.
#' 
#' @param red A numeric vector (or array) of red channel values. If this is the
#'   only channel argument, it can also be a character vector of colour values
#'   (including alpha, if required), or a numeric array whose last dimension is
#'   2 (for grey + alpha), 3 (for RGB) or 4 (for RGBA).
#' @param green,blue,alpha Numeric vectors (or arrays) containing values for
#'   the appropriate channel. These will be combined with the \code{red} values
#'   using \code{cbind}, and hence recycled as necessary. Alpha, or green and
#'   blue, can be missing.
#' @param max The maximum possible value for any channel. The default
#'   is 255 when the data is of integer mode, and 1 otherwise.
#' @param dim An integer vector of dimensions for the final array. The
#'   dimensions of \code{red} are used if this is \code{NULL}.
#' @param ... For \code{rgbArray}, additional attributes to set on the result,
#'   such as \code{pixdim}. These are passed directly to
#'   \code{\link{structure}}. For the \code{as.character} method, this argument
#'   is ignored.
#' @param x An \code{rgbArray} object.
#' @param flatten Logical value. If \code{FALSE}, the dimensions of \code{x}
#'   will be retained in the result. The default is \code{TRUE}, for
#'   consistency with the usual behaviour of \code{as.character}, which strips
#'   all attributes.
#' @return \code{rgbArray} returns an integer-mode array of class
#'   \code{"rgbArray"}. The \code{as.character} method returns a character-mode
#'   vector of colour strings.
#' 
#' @note The values of an \code{"rgbArray"} are not easily interpreted, and
#'   may depend on the endianness of the platform. For manipulation or use as
#'   colours they should generally be converted to character mode, or the
#'   channels extracted using the \code{\link{channels}} function.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
rgbArray <- function (red, green, blue, alpha, max = NULL, dim = NULL, ...)
{
    source <- NULL
    channels <- 0L
    
    if (!missing(green) && !missing(blue) && !missing(alpha))
    {
        source <- cbind(red, green, blue, alpha)
        channels <- 4L
    }
    else if (!missing(green) && !missing(blue))
    {
        source <- cbind(red, green, blue)
        channels <- 3L
    }
    else if (!missing(alpha))
    {
        source <- cbind(red, alpha)
        channels <- 2L
    }
    else if (is.character(red))
    {
        channels <- ifelse(max(nchar(red),na.rm=TRUE) > 7L, 4L, 3L)
        source <- t(col2rgb(red, alpha=(channels==4L)))
    }
    else if (is.numeric(red) && is.array(red))
    {
        source <- red
        channels <- dim(red)[ndim(red)]
        if (channels < 2L || channels > 4L)
            stop("If only one numeric argument is supplied, its last dimension must be 2, 3 or 4")
        if (is.null(dim))
            dim <- dim(red)[-ndim(red)]
    }
    else
        stop("The combination of channels provided is not supported")
    
    if (is.null(dim))
        dim <- dim(red)
    if (is.null(max))
        max <- switch(storage.mode(source), integer=255, 1)
    
    result <- .Call("packRgb", source, channels, max, PACKAGE="RNifti")
    return (structure(result, ..., channels=channels, dim=dim, class="rgbArray"))
}

#' @rdname rgbArray
#' @export
as.character.rgbArray <- function (x, flatten = TRUE, ...)
{
    result <- .Call("rgbToStrings", x, PACKAGE="RNifti")
    if (!flatten)
        dim(result) <- dim(x)
    return (result)
}

#' Extract channels from RGB data
#' 
#' Extract one or more channels from an RGB data array that was obtained from
#' an RGB NIfTI image or created by the \code{\link{rgbArray}} function. The
#' result is more amenable to numeric manipulation.
#' 
#' @param array An image, an \code{rgbArray}, or another array that can be
#'   converted to the latter.
#' @param channels A character vector of channels to extract.
#' @param raw Boolean value: if \code{TRUE}, return a raw array, which is the
#'   most compact representation; otherwise return an integer array.
#' @return A raw-mode or integer-mode array with one more dimension than the
#'   first argument, corresponding to channels.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
channels <- function (array, channels = c("red","green","blue","alpha"), raw = FALSE)
{
    if (!inherits(array,"niftiImage") && !inherits(array,"rgbArray"))
        array <- rgbArray(array)
    
    channels <- match.arg(channels, several.ok=TRUE)
    channelNumbers <- c(red=1L, green=2L, blue=3L, alpha=4L)[channels]
    
    result <- .Call("unpackRgb", array, channelNumbers, PACKAGE="RNifti")
    if (!raw)
        storage.mode(result) <- "integer"
    dimnames(result) <- c(rep(list(NULL),ndim(array)), list(channels))
    return (result)
}
