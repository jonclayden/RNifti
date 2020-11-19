#' NIfTI extensions
#' 
#' The NIfTI-1 and NIfTI-2 formats have a simple extension mechanism that
#' allows additional metadata to be stored with their headers. The format of
#' this extension data is unspecified by the NIfTI standard, but extension
#' codes indicate what type of information is present. These functions provide
#' access to this extension metadata.
#' 
#' @param image An image, in any acceptable form (see \code{\link{asNifti}}).
#' @param code Integer value specifying which extension code is required.
#' @param mode The required mode of the extracted data.
#' @param ... Additional arguments to \code{\link{readBin}}.
#' @param simplify Logical value. If \code{TRUE}, the default, a single
#'   extension will be returned as a vector; otherwise a list is always
#'   returned.
#' @return For \code{extensions}, a list of raw vectors containing the bytes
#'   stored in each available header. For \code{extension}, a list of vector
#'   of values, converted to the required mode, for the extension code
#'   specified. If the extension code is not used in the image, the return
#'   value is \code{NULL}.
#' 
#' @note A list of registered extension codes is available in the
#'   \code{nifti2_io.h} header file, which is distributed with this package.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
extensions <- function (image)
{
    .Call("getExtensions", asNifti(image,internal=TRUE), 0L, PACKAGE="RNifti")
}

#' @rdname extensions
#' @export
extension <- function (image, code, mode = c("raw","character","numeric","double","integer","logical","complex"), ..., simplify = TRUE)
{
    mode <- match.arg(mode)
    raw <- .Call("getExtensions", asNifti(image,internal=TRUE), code, PACKAGE="RNifti")
    
    convert <- function(x) switch(mode, raw=x, character=rawToChar(x), readBin(rawConnection(x), mode, length(x), ...))
    
    if (length(raw) == 0L)
        return (NULL)
    else if (simplify && length(raw) == 1L)
        return (convert(raw[[1]]))
    else
        return (lapply(raw, convert))
}
