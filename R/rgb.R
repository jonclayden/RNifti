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
            stop("If only one argument is supplied, its last dimension must be 2, 3 or 4")
        if (is.null(dim))
            dim <- dim(red)[-ndim(red)]
    }
    else
        stop("At least two channels must be specified, or the first argument should be multidimensional")
    
    if (is.null(dim))
        dim <- dim(red)
    if (is.null(max))
        max <- switch(storage.mode(source), integer=255, 1)
    
    result <- .Call("packRgb", source, channels, max, PACKAGE="RNifti")
    return (structure(result, ..., channels=channels, dim=dim, class="rgbArray"))
}

#' @export
as.character.rgbArray <- function (x, ...)
{
    .Call("rgbToStrings", x, PACKAGE="RNifti")
}

#' @export
extractChannels <- function (array, channels = c("red","green","blue","alpha"), raw = FALSE)
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
