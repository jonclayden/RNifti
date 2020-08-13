#' Read NIfTI or ANALYZE format files
#' 
#' This function reads one or more NIfTI-1, NIfTI-2 or ANALYZE-7.5 files into
#' R, using the standard NIfTI C library.
#' 
#' @param file A character vector of file names.
#' @param internal Logical value. If \code{FALSE} (the default), an array
#'   of class \code{"niftiImage"}, containing the image pixel or voxel values,
#'   will be returned. If \code{TRUE}, the return value will be an object of
#'   class \code{"internalImage"}, which contains only minimal metadata about
#'   the image. Either way, the return value has an attribute which points to a
#'   C data structure containing the full image.
#' @param volumes An integer vector giving the volumes to read (counting along
#'   all dimensions beyond the third jointly), or \code{NULL}, the default, in
#'   which case every volume is read. This cannot currently be set differently
#'   for each file read.
#' @return An array or internal image, with class \code{"niftiImage"} (and
#'   possibly also \code{"internalImage"}), or a list of such objects if
#'   \code{file} has length greater than one.
#' 
#' @note If the \code{internal} argument is \code{FALSE} (the default), the
#'   data type of the image pointer will be set to match one of R's native
#'   numeric data types, i.e., 32-bit signed integer or 64-bit double-precision
#'   floating-point. In these circumstances the data type reported by the
#'   \code{\link{niftiHeader}} function will therefore not, in general, match
#'   the storage type used in the file. See also the \code{datatype} argument
#'   to \code{\link{writeNifti}}.
#' 
#' @examples
#' path <- system.file("extdata", "example.nii.gz", package="RNifti")
#' readNifti(path)
#' readNifti(path, internal=TRUE)
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{writeNifti}}
#' @references The NIfTI-1 standard (\url{http://www.nitrc.org/docman/view.php/26/64/nifti1.h}).
#' @aliases readAnalyze
#' @export readNifti readAnalyze
readNifti <- readAnalyze <- function (file, internal = FALSE, volumes = NULL)
{
    if (!is.character(file))
        stop("File name(s) must be specified in a character vector")
    if (length(file) == 0)
        stop("File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", f, internal, volumes, PACKAGE="RNifti"))
    else
        .Call("readNifti", file, internal, volumes, PACKAGE="RNifti")
}

#' Write a NIfTI or ANALYZE format file
#' 
#' These functions write an image to NIfTI-1, NIfTI-2 or ANALYZE-7.5 format,
#' using the standard NIfTI C library.
#' 
#' @param image An image, in any acceptable form (see \code{\link{asNifti}}).
#' @param file A character string containing a file name.
#' @param template An optional template object to derive NIfTI-1 properties
#'   from. Passed to \code{\link{asNifti}} if \code{image} is an array.
#' @param datatype The NIfTI datatype to use when writing the data out. The
#'   default, \code{"auto"} uses the R type or, for internal images, the
#'   original datatype. Other possibilities are \code{"float"}, \code{"int16"},
#'   etc., which may be preferred to reduce file size. However, no checks are
#'   done to ensure that the coercion maintains precision.
#' @param version An integer (1 or 2) giving the NIfTI file format version to
#'   use. Version 2 is usually only needed for very large images or where
#'   metadata needs to be stored with high precision. The types available for
#'   storing the pixel data are the same in both cases.
#' @return An invisible, named character vector giving the image and header
#'   file names written to.
#' 
#' @note The ANALYZE-7.5 file format is a legacy format and use of it is not
#'   recommended, except for compatibility. In particular, the format does
#'   not reliably specify the spatial orientation of the image.
#' 
#' @examples
#' \dontrun{writeNifti(im, "image.nii.gz", datatype="float")}
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{asNifti}}
#' @references The NIfTI-1 standard (\url{http://www.nitrc.org/docman/view.php/26/64/nifti1.h}).
#' @export
writeNifti <- function (image, file, template = NULL, datatype = "auto", version = 1)
{
    invisible(.Call("writeNifti", asNifti(image,template,internal=TRUE), file, tolower(datatype), switch(version,"nifti1","nifti2"), PACKAGE="RNifti"))
}

#' @rdname writeNifti
#' @export
writeAnalyze <- function (image, file, template = NULL, datatype = "auto")
{
    invisible(.Call("writeNifti", asNifti(image,template,internal=TRUE), file, tolower(datatype), "analyze", PACKAGE="RNifti"))
}

#' Create or modify an NIfTI image object
#' 
#' This function converts a filename, array or other image class into an object
#' of class \code{"niftiImage"}, and optionally updates its metadata from a
#' reference image and/or changes its internal datatype. The dimensions and
#' pixel dimensions from the \code{image} will replace those from the reference
#' object, if they are available.
#' 
#' If the \code{image} has an internal NIfTI pointer, that will be retrieved
#' directly. Otherwise, if it is a string, it will be taken to be a filename.
#' If it looks like a \code{"nifti"} object (from package \code{oro.nifti}),
#' or an \code{"MriImage"} object (from package \code{tractor.base}), a
#' conversion will be performed. A list will be assumed to be of the form
#' produced by \code{\link{niftiHeader}}. Finally, a numeric array or matrix,
#' or RGB array, will be converted using default image parameters.
#' 
#' If \code{reference} is a complete list of NIfTI-1 header fields, like that
#' produced by \code{\link{niftiHeader}}, or an image, then it will be used to
#' create the internal object, and then the data and metadata associated with
#' the \code{image} will overwrite the appropriate parts. If \code{reference}
#' is an incomplete list, the \code{image} will be used to create the internal
#' object, and then the specified fields will be overwritten from the list.
#' This allows users to selectively update certain fields while leaving others
#' alone (but see the note below).
#' 
#' If multiple values are passed for a field that expects a scalar (which is
#' most of them), the first element of the vector will be used, with a warning.
#' An empty vector will be ignored, also with a warning. If a value of the
#' wrong length is passed to a vector-valued field, an error will be generated.
#' 
#' Datatype information in a list \code{reference} is ignored. The datatype can
#' only be changed using the \code{datatype} argument, but in this case the
#' internal object gets out of sync with the R array, so an internal image is
#' returned to avoid the mismatch. Changing the internal datatype in this way
#' is for advanced usage only.
#' 
#' \code{retrieveNifti} and \code{updateNifti} are soft-deprecated alternative
#' interfaces to this function, which behave like the pre-existing functions of
#' the same names. They may be removed in future.
#' 
#' @param x Any suitable object (see Details).
#' @param reference An image, or a named list of NIfTI-1 properties like that
#'   produced by \code{\link{niftiHeader}}. The default of \code{NULL} will
#'   have no effect.
#' @param ... Additional parameters to methods.
#' @param datatype The NIfTI datatype to use within the internal image. The
#'   default, \code{"auto"} uses the R type. Other possibilities are
#'   \code{"float"}, \code{"int16"}, etc., which may be preferred to reduce
#'   object size. However, no checks are done to ensure that the coercion
#'   maintains precision, and this option is for advanced usage only.
#' @param internal Logical value. If \code{FALSE}, the result will be an array
#'   of class \code{"niftiImage"} containing the image pixel or voxel values,
#'   with some metadata in attributes. If \code{TRUE}, the result will be an
#'   object of class \code{"internalImage"}, which exposes some basic metadata
#'   to R but stores the pixel data internally. If \code{NA}, the default, the
#'   result will be an internal image only if the input \code{image} is. If a
#'   new \code{datatype} is set then this value is implicitly \code{TRUE}.
#' @return An array or internal image, with class \code{"niftiImage"} (and
#'   possibly also \code{"internalImage"}).
#' 
#' @note The \code{scl_slope} and \code{scl_inter} fields affect the numerical
#'   interpretation of the pixel data, so it is impossible in general to change
#'   them without also changing the array values on both the C and the R side.
#'   Therefore, to avoid unexpected side-effects, these fields are not affected
#'   by this function. The \code{dim} and \code{pixdim} fields can be changed,
#'   but for most users the accessor functions of the same name are much safer,
#'   and should be used in preference.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{$.niftiImage}},
#'   \code{\link{dim.internalImage}}, \code{\link{pixdim}}, \code{\link{xform}}
#' @aliases retrieveNifti updateNifti
#' @export
asNifti <- function (x, reference = NULL, ...)
{
    UseMethod("asNifti")
}

#' @rdname asNifti
#' @export
asNifti.default <- function (x, reference = NULL, datatype = "auto", internal = NA, ...)
{
    .Call("asNifti", x, reference, datatype, internal, PACKAGE="RNifti")
}

#' @export
retrieveNifti <- function (object)
{
    asNifti(object, internal=TRUE)
}

#' @export
updateNifti <- function (image, template = NULL, datatype = "auto")
{
    asNifti(image, template, datatype, internal=FALSE)
}

#' Dump or construct a raw NIfTI or ANALYZE header
#' 
#' These functions extract the contents of a NIfTI-1 or ANALYZE-7.5 header,
#' closely approximating how it is (or would be) stored on disk. Defaults will
#' be used where information is missing, but no processing is performed on the
#' metadata.
#' 
#' The NIfTI-1 standard was originally formulated as a roughly backwards-
#' compatible improvement on the ANALYZE format. Both formats use a binary
#' header structure of 348 bytes, but the field names and their interpretation
#' is often non-equivalent. These functions dump these fields, without regard
#' to whether or not the result makes proper sense.
#' 
#' \code{dumpNifti} is an alias of \code{niftiHeader}, but the former is now
#' soft-deprecated.
#' 
#' @param image An image, in any acceptable form (see \code{\link{asNifti}}).
#'   A list containing partial header information is acceptable, including an
#'   empty list, which returns defaults for every field.
#' @param x A \code{"niftiHeader"} object.
#' @param ... Ignored.
#' @return For \code{niftiHeader}, a list of class \code{"niftiHeader"}, with
#'   named components corresponding to the elements in a raw NIfTI-1 header.
#'   For \code{analyzeHeader}, the equivalent for ANALYZE-7.5.
#' 
#' @note Several medical image analysis packages, such as SPM and FSL, use the
#'   ANALYZE \code{originator} field to store a coordinate origin. This
#'   interpretation is also returned, in the \code{origin} field.
#' 
#' @examples
#' niftiHeader(system.file("extdata", "example.nii.gz", package="RNifti"))
#' 
#' # Default header for a standard R array
#' niftiHeader(array(0L, dim=c(10,10)))
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftiVersion}}
#' @references The NIfTI-1 standard (\url{http://www.nitrc.org/docman/view.php/26/64/nifti1.h}).
#' @aliases dumpNifti
#' @export niftiHeader dumpNifti
niftiHeader <- dumpNifti <- function (image = list())
{
    .Call("niftiHeader", asNifti(image,internal=TRUE), PACKAGE="RNifti")
}

#' @rdname niftiHeader
#' @export
analyzeHeader <- function (image = list())
{
    .Call("analyzeHeader", asNifti(image,internal=TRUE), PACKAGE="RNifti")
}

#' @rdname niftiHeader
#' @export
print.niftiHeader <- function (x, ...)
{
    cat("NIfTI-1 header\n")
    widths <- nchar(names(x), "width")
    maxWidth <- max(widths)
    names <- names(x)
    strings <- attr(x, "strings")
    
    for (i in seq_along(widths))
    {
        if (!is.null(strings[[names[i]]]))
            cat(paste0(paste(rep(" ",maxWidth-widths[i]),collapse=""), names[i], ": ", paste(format(x[[i]],trim=TRUE),collapse="  "), " (", strings[[names[i]]], ")\n"))
        else
            cat(paste0(paste(rep(" ",maxWidth-widths[i]),collapse=""), names[i], ": ", paste(format(x[[i]],trim=TRUE),collapse="  "), "\n"))
    }
}

#' @rdname niftiHeader
#' @export
print.analyzeHeader <- function (x, ...)
{
    cat("ANALYZE-7.5 header\n")
    widths <- nchar(names(x), "width")
    maxWidth <- max(widths)
    
    for (i in seq_along(widths))
        cat(paste0(paste(rep(" ",maxWidth-widths[i]),collapse=""), names(x)[i], ": ", paste(format(x[[i]],trim=TRUE),collapse="  "), "\n"))
}

#' Check the format version of a file
#' 
#' This function identifies the likely NIfTI format variant used by one or more
#' files on disk.
#' 
#' @param file A character vector of file names.
#' @return A vector of integers, of the same length as \code{file}. Each
#'   element will be 0 for ANALYZE format (the precursor to NIfTI-1), 1 for
#'   NIfTI-1 (which is now most common), 2 for NIfTI-2, or -1 if the file
#'   doesn't exist or doesn't look plausible in any of these formats.
#' 
#' @note NIfTI-2 format, mostly a variant of NIfTI-1 with wider datatypes used
#'   for many fields, is not currently supported for reading, but it is
#'   detected by this function.
#' 
#' @examples
#' path <- system.file("extdata", "example.nii.gz", package="RNifti")
#' niftiVersion(path)       # 1
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{niftiHeader}}
#' @export
niftiVersion <- function (file)
{
    sapply(file, function(f) .Call("niftiVersion", f, PACKAGE="RNifti"))
}

addresses <- function (image)
{
    .Call("getAddresses", image, PACKAGE="RNifti")
}

hasData <- function (image)
{
    .Call("hasData", image, PACKAGE="RNifti")
}

rescaleNifti <- function (image, scales)
{
    .Call("rescaleImage", image, scales, PACKAGE="RNifti")
}

niftiDebug <- function (level = 1L)
{
    invisible(.Call("setDebugLevel", level, PACKAGE="RNifti"))
}
