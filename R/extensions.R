#' @export
ExtensionCodes <- list(ignore=0L, DICOM=2L, AFNI=4L, comment=6L, XCEDE=8L, jimdiminfo=10L, workflow_fwds=12L, FreeSurfer=14L, pypickle=16L, MiND_ident=18L, b_value=20L, spherical_direction=22L, DT_component=24L, SHC_degreeorder=26L, voxbo=28L, Caret=30L, CIFTI=32L, variable_frame_timing=34L, eval=38L, MATLAB=40L, Quantiphyse=42L, MRS=44L)

.resolveCode <- function (code)
{
    name <- as.character(substitute(code, parent.frame()))
    loc <- match(tolower(name), tolower(names(ExtensionCodes)))
    if (length(name) == 1L && !is.na(loc))
        return (ExtensionCodes[[loc]])
    else if (is.numeric(code) && code < 0)
        stop("Extension code must not be negative")
    else
        return (code)
}

#' NIfTI extensions
#' 
#' The NIfTI-1 and NIfTI-2 formats have a simple extension mechanism that
#' allows additional metadata to be stored with their headers. The format of
#' this extension data is unspecified by the NIfTI standard, but extension
#' codes indicate what type of information is present. These functions provide
#' access to this extension metadata.
#' 
#' The plural version, \code{extensions}, extracts or replaces all extensions
#' at once. The retrieval form returns a list of raw vectors, each with the
#' corresponding code in an attribute, and the replacement form accepts a list
#' of atomic vectors with code attributes, or \code{NULL}, which removes all
#' extensions. The singular version, \code{extension}, gets all extensions with
#' the specified code, or appends an extension with that code.
#' 
#' NIfTI extensions are stored as a simple, unstructured byte stream, which is
#' naturally represented in R as a vector of mode \code{"raw"}. However, these
#' functions will perform some conversion to and from other atomic types for
#' convenience. The NIfTI standard makes no guarantees about byte order within
#' the data stream, but the \code{endian} argument to \code{\link{readBin}} can
#' be passed through when converting to a non-raw type.
#' 
#' @param image An image, in any acceptable form (see \code{\link{asNifti}}).
#' @param code Integer value specifying which extension code is required.
#' @param mode The required mode of the extracted data.
#' @param ... Additional arguments to \code{\link{readBin}}.
#' @param simplify Logical value. If \code{TRUE}, the default, a single
#'   extension will be returned as a vector; otherwise a list is always
#'   returned.
#' @param value New value for the extension(s).
#' @return For \code{extensions}, a list of raw vectors containing the bytes
#'   stored in each available header. For \code{extension}, a list of vector
#'   of values, converted to the required mode, for the extension code
#'   specified. If the extension code is not used in the image, the return
#'   value is \code{NULL}. The replacement forms return the modified image.
#' 
#' @note A list of registered extension codes is available in the
#'   \code{nifti2_io.h} header file, which is distributed with this package.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
extensions <- function (image)
{
    .Call("getExtensions", asNifti(image,internal=TRUE), -1L, PACKAGE="RNifti")
}

#' @rdname extensions
#' @export
extension <- function (image, code, mode = c("raw","character","numeric","double","integer","logical","complex"), ..., simplify = TRUE)
{
    mode <- match.arg(mode)
    raw <- .Call("getExtensions", asNifti(image,internal=TRUE), .resolveCode(code), PACKAGE="RNifti")
    
    convert <- function(x) switch(mode, raw=x, character=rawToChar(x), readBin(rawConnection(x), mode, length(x), ...))
    
    if (length(raw) == 0L)
        return (NULL)
    else if (simplify && length(raw) == 1L)
        return (convert(raw[[1]]))
    else
        return (lapply(raw, convert))
}

#' @rdname extensions
#' @export
`extensions<-` <- function (image, value)
{
    if (!is.list(value))
        value <- list(value)
    .Call("setExtensions", asNifti(image), value, -1L, PACKAGE="RNifti")
}

#' @rdname extensions
#' @export
`extension<-` <- function (image, code, value)
{
    if (is.list(value))
    {
        if (length(value) > 1L)
            warning("List elements after the first are ignored when adding an image extension")
        value <- value[[1]]
    }
    .Call("setExtensions", asNifti(image), value, .resolveCode(code), PACKAGE="RNifti")
}
