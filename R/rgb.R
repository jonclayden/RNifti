rgbArray <- function (r, g, b, a, max = NULL, dim = NULL, ...)
{
    if (is.null(dim))
        dim <- dim(r)
    
    source <- NULL
    channels <- 0L
    
    if (!missing(g) && !missing(b) && !missing(a))
    {
        source <- c(r, g, b, a)
        channels <- 4L
    }
    else if (!missing(g) && !missing(b))
    {
        source <- c(r, g, b)
        channels <- 4L
    }
    else if (is.matrix(r) || is.array(r))
    {
        source <- r
        channels <- dim(r)[ndim(r)]
    }
    
    if (is.null(max))
        max <- switch(storage.mode(source), integer=255, 1)
    
    result <- .Call("packRgb", source, channels, max, PACKAGE="RNifti")
    return (structure(result, ..., dim=dim, class="rgbArray"))
}
