#' Read a NIfTI-1 format file
#' 
#' This function reads one or more NIfTI-1 files into R, using the standard
#' NIfTI-1 C library.
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
#'   \code{\link{dumpNifti}} function will therefore not, in general, match
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
#' @export
readNifti <- function (file, internal = FALSE, volumes = NULL)
{
    if (!is.character(file))
        stop("File name(s) must be specified in a character vector")
    if (length(file) == 0)
        stop("File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", path.expand(f), internal, volumes, PACKAGE="RNifti"))
    else
        .Call("readNifti", path.expand(file), internal, volumes, PACKAGE="RNifti")
}

#' Write a NIfTI-1 format file
#' 
#' This function writes an image to NIfTI-1 format, using the standard NIfTI-1
#' C library.
#' 
#' @param image An image, in any acceptable form (see
#'   \code{\link{retrieveNifti}}).
#' @param file A character string containing a file name.
#' @param template An optional template object to derive NIfTI-1 properties
#'   from. Passed to \code{\link{updateNifti}} if \code{image} is an array.
#' @param datatype The NIfTI datatype to use when writing the data out. The
#'   default, \code{"auto"} uses the R type or, for internal images, the
#'   original datatype. Other possibilities are \code{"float"}, \code{"int16"},
#'   etc., which may be preferred to reduce file size. However, no checks are
#'   done to ensure that the coercion maintains precision.
#' 
#' @examples
#' \dontrun{writeNifti(im, "image.nii.gz", datatype="float")}
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{updateNifti}}
#' @references The NIfTI-1 standard (\url{http://www.nitrc.org/docman/view.php/26/64/nifti1.h}).
#' @export
writeNifti <- function (image, file, template = NULL, datatype = "auto")
{
    if (is.array(image) && !is.null(template))
        image <- updateNifti(image, template)
    
    invisible(.Call("writeNifti", image, path.expand(file), tolower(datatype), PACKAGE="RNifti"))
}

#' Obtain an internal NIfTI representation of an object
#' 
#' This function converts filenames, arrays and other image classes into an
#' object of class \code{"internalImage"}.
#' 
#' If the \code{object} has an internal NIfTI pointer, that will be retrieved
#' directly. Otherwise, if it is a string, it will be taken to be a filename.
#' If it looks like a \code{"nifti"} object (from package \code{oro.nifti}),
#' or an \code{"MriImage"} object (from package \code{tractor.base}), a
#' conversion will be attempted. A list will be assumed to be of the form
#' produced by \code{\link{dumpNifti}}. Finally, a numeric array or matrix
#' will be converted using default image parameters.
#'
#' @param object Any suitable object (see Details).
#' @return An internal image.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{updateNifti}}
#' @export
retrieveNifti <- function (object)
{
    .Call("readNifti", object, TRUE, NULL, PACKAGE="RNifti")
}

#' Update an internal NIfTI-1 object using a template
#' 
#' This function adds or updates the internal NIfTI-1 object for an array,
#' using metadata from the template. The dimensions and, if available, pixel
#' dimensions, from the \code{image} will replace those from the template.
#' 
#' If \code{template} is a complete list of NIfTI-1 header fields, like that
#' produced by \code{\link{dumpNifti}}, or an image, then it will be used to
#' create the internal object, and then the data and metadata associated with
#' the \code{image} will overwrite the appropriate parts. If \code{template}
#' is an incomplete list, the \code{image} will be used to create the internal
#' object, and then the specified fields will be overwritten from the list.
#' This allows users to selectively update certain fields while leaving others
#' alone (but see the note below).
#' 
#' Datatype information in a list \code{template} is ignored. The datatype can
#' only be changed using the \code{datatype} argument, but in this case the
#' internal object gets out of sync with the R array, so an internal image is
#' returned to avoid the mismatch. Changing the internal datatype in this way
#' is for advanced usage only.
#' 
#' @param image A numeric array.
#' @param template An image, in any acceptable form (see
#'   \code{\link{retrieveNifti}}), or a named list of NIfTI-1 properties like
#'   that produced by \code{\link{dumpNifti}}. The default of \code{NULL} will
#'   have no effect.
#' @param datatype The NIfTI datatype to use within the internal image. The
#'   default, \code{"auto"} uses the R type. Other possibilities are
#'   \code{"float"}, \code{"int16"}, etc., which may be preferred to reduce
#'   object size. However, no checks are done to ensure that the coercion
#'   maintains precision, and this option is for advanced usage only.
#' @return A copy of the original \code{image}, with its internal image
#'   attribute set or updated appropriately. If \code{datatype} is not
#'   \code{"auto"} then the result is an internal image.
#' 
#' @note The \code{scl_slope} and \code{scl_inter} fields affect the numerical
#'   interpretation of the pixel data, so it is impossible in general to change
#'   them without also changing the array values on both the C and the R side.
#'   Therefore, to avoid unexpected side-effects, these fields are not affected
#'   by this function.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
updateNifti <- function (image, template = NULL, datatype = "auto")
{
    .Call("updateNifti", image, template, datatype, PACKAGE="RNifti")
}

#' Dump the contents of an internal NIfTI-1 object
#' 
#' This function extracts the contents of an internal NIfTI-1 object into an R
#' list. No processing is done to the elements.
#' 
#' @param image An image, in any acceptable form (see
#'   \code{\link{retrieveNifti}}).
#' @param x A \code{"niftiHeader"} object.
#' @param ... Ignored.
#' @return For \code{dumpNifti}, a list of class \code{"niftiHeader"}, with
#'   named components corresponding to the elements in a raw NIfTI-1 file.
#' 
#' @examples
#' dumpNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' 
#' # Default header for a standard R array
#' dumpNifti(array(0L, dim=c(10,10)))
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{http://www.nitrc.org/docman/view.php/26/64/nifti1.h}).
#' @export
dumpNifti <- function (image)
{
    .Call("dumpNifti", image, PACKAGE="RNifti")
}

#' @rdname dumpNifti
#' @export
print.niftiHeader <- function (x, ...)
{
    cat("NIfTI-1 header\n")
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
#' @seealso \code{\link{readNifti}}
#' @export
niftiVersion <- function (file)
{
    sapply(file, function(f) .Call("niftiVersion", path.expand(f), PACKAGE="RNifti"))
}

rescaleNifti <- function (image, scales)
{
    .Call("rescaleImage", image, scales, PACKAGE="RNifti")
}
