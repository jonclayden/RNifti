.RNiftiAttribs <- "^\\.|^(image|pix)?dim$|^pixunits$|^class$"

.Bids <- list(
    mappingFromJson=c(MagneticFieldStrength="fieldStrength",
                      ManufacturersModelName="scannerModelName",
                      SpacingBetweenSlices="sliceSpacing",
                      TotalReadoutTime="effectiveReadoutTime",
                      MultibandAccelerationFactor="multibandFactor",
                      ImageComments="comments"),
    mappingToJson=c(fieldStrength="MagneticFieldStrength",
                    scannerModelName="ManufacturersModelName",
                    sliceSpacing="SpacingBetweenSlices",
                    effectiveReadoutTime="TotalReadoutTime",
                    multibandFactor="MultibandAccelerationFactor",
                    comments="ImageComments"),
    toScale="^(Echo|Repetition|Inversion)Time$")

#' Extended image attributes
#' 
#' These functions extract and replace medical image attributes that go beyond
#' the core metadata associated with the NIfTI file formats.
#' 
#' The DICOM (Digital Imaging and Communications in Medicine) format, and BIDS
#' (Brain Imaging Data Structure), which extends NIfTI, can both encapsulate
#' copious amounts of metadata about a scan and the patient. This metadata can
#' be useful for more advanced or research-focussed post-processing methods,
#' and standard R attributes are a natural place to store it. The
#' \code{imageAttributes} function returns a list of just these extended
#' attributes, if they exist, ignoring other attributes used by the package.
#' The replacement form allows this metadata to be modified or removed. These
#' functions currently only act on objects inheriting from the
#' \code{niftiImage} class.
#' 
#' @param x A \code{niftiImage} object.
#' @param value A list of new image attributes to replace any existing ones.
#' @return A list of image attributes, or a modified object with these changed.
#'   These are essentially all attributes except those used for basic
#'   \code{niftiImage} objects by the package.
#' 
#' @examples
#' path <- system.file("extdata", "raw", package="divest")
#' images <- readDicom(path, interactive=FALSE)
#' imageAttributes(images[[1]])
#' @author Jon Clayden <code@@clayden.org>
#' @seealso The \code{divest} package can convert DICOM files to NIfTI formats,
#'   and extract embedded metadata. More information about DICOM is available
#'   at \url{https://www.dicomstandard.org}.
#' @references More information about metadata captured by the BIDS format can
#'   be found at \url{https://bids.neuroimaging.io} or in the paper cited
#'   below.
#' 
#' K.J. Gorgolewski, T. Auer, V.D. Calhoun, et al. The brain imaging data
#' structure, a format for organizing and describing outputs of neuroimaging
#' experiments (2016). Scientific Data 3:160044. \doi{10.1038/sdata.2016.44}.
#' @export
imageAttributes <- function (x)
{
    if (!inherits(x, "niftiImage"))
        return (NULL)
    attribs <- attributes(x)
    if (length(attribs) == 0L || is.null(names(attribs)))
        return (NULL)
    attribs <- attribs[!grepl(.RNiftiAttribs,names(attribs),perl=TRUE) & names(attribs) != ""]
    return (attribs)
}

#' @rdname imageAttributes
#' @export
`imageAttributes<-` <- function (x, value)
{
    if (!inherits(x, "niftiImage"))
        return (x)
    attribs <- attributes(x)
    if (length(attribs) == 0L || is.null(names(attribs)))
        return (x)
    attribs <- attribs[grepl(.RNiftiAttribs,names(attribs),perl=TRUE) | names(attribs) == ""]
    attributes(x) <- c(attribs, as.list(value))
    return (x)
}

renameFromBids <- function (bids)
{
    tractor <- list()
    
    passthrough <- grepl(.RNiftiAttribs, names(bids), perl=TRUE)
    if (any(passthrough))
        tractor <- bids[passthrough]
    
    for (bidsName in names(bids)[!passthrough])
    {
        value <- bids[[bidsName]]
        if (is.character(value) && all(grepl("^\\s*$", value, perl=TRUE)))
            next
        
        if (bidsName == "PhaseEncodingDirection" && grepl("^([ijk])(-)?$", value, perl=TRUE))
        {
            tractor$phaseEncodingDirection <- substring(value,1,1)
            tractor$phaseEncodingSign <- ifelse(nchar(value)==1L, 1L, -1L)
            next
        }
        
        # Attribute names mostly reflect the BIDS ones, but with the first letter downcased
        # A few exceptions are mapped explicitly
        if (bidsName %in% names(.Bids$mappingFromJson))
            name <- .Bids$mappingFromJson[[bidsName]]
        else if (grepl("^[A-Z][a-z]", bidsName, perl=TRUE))
            name <- paste0(tolower(substring(bidsName,1,1)), substring(bidsName,2))
        else
            name <- bidsName
        
        # BIDS always uses seconds for time fields, but we use milliseconds in some places
        if (grepl(.Bids$toScale, bidsName, perl=TRUE))
            value <- value * 1e3
        tractor[[name]] <- value
    }
    
    return (tractor)
}

renameToBids <- function (tractor)
{
    bids <- list()
    
    passthrough <- grepl(.RNiftiAttribs, names(tractor), perl=TRUE)
    if (any(passthrough))
        bids <- tractor[passthrough]
    
    if (all(c("phaseEncodingDirection","phaseEncodingSign") %in% names(tractor)))
    {
        bids$PhaseEncodingDirection <- paste0(tractor$phaseEncodingDirection, ifelse(tractor$phaseEncodingSign < 0,"-",""))
        tractor <- tractor[setdiff(names(tractor), c("phaseEncodingDirection","phaseEncodingSign"))]
    }
    
    for (name in names(tractor)[!passthrough])
    {
        value <- tractor[[name]]
        if (is.character(value) && all(grepl("^\\s*$", value, perl=TRUE)))
            next
        
        if (name %in% names(.Bids$mappingToJson))
            bidsName <- .Bids$mappingToJson[[name]]
        else if (grepl("^[a-z]", name, perl=TRUE))
            bidsName <- paste0(toupper(substring(name,1,1)), substring(name,2))
        else
            bidsName <- name
        
        if (grepl(.Bids$toScale, bidsName, perl=TRUE))
            value <- value / 1e3
        bids[[bidsName]] <- value
    }
    
    return (bids)
}

#' Conversion to and from BIDS JSON
#' 
#' Functions to convert to and from BIDS JSON format for image metadata. They
#' are wrappers around functions from the \code{jsonlite} package, with the
#' additional ability to convert between the tag naming convention used by the
#' \code{divest} and \code{tractor.base} packages and the BIDS equivalent. The
#' differences are mostly in capitalisation, and the units used for magnetic
#' resonance echo, repetition and inversion times.
#' 
#' @param source A list containing metadata (see \code{\link{imageAttributes}})
#'   or, for \code{fromBidsJson}, a string containing literal JSON or the path
#'   to a file containing it.
#' @param rename Logical value. If \code{TRUE}, element names are also
#'   converted to or from the BIDS convention; otherwise this is just a
#'   conversion between an R list and a JSON string.
#' @param path For \code{toBidsJson}, the path to write the JSON output to. If
#'   \code{NULL}, the default, the JSON text is returned in an object.
#' @return \code{fromBidsJson} returns a list of image attributes.
#'   \code{toBidsJson} returns a character vector if \code{path} is
#'   \code{NULL}, otherwise nothing.
#' 
#' @note These functions do not check BIDS metadata for validity in either
#'   direction, either in terms of the names of the fields or their contents.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references More information about metadata captured by the BIDS format can
#'   be found at \url{https://bids.neuroimaging.io} or in the paper cited
#'   below.
#' 
#' K.J. Gorgolewski, T. Auer, V.D. Calhoun, et al. The brain imaging data
#' structure, a format for organizing and describing outputs of neuroimaging
#' experiments (2016). Scientific Data 3:160044. \doi{10.1038/sdata.2016.44}.
#' @rdname bidsJson
#' @export
fromBidsJson <- function (source, rename = FALSE)
{
    if (length(source) == 0L)
        return (list())
    result <- jsonlite::fromJSON(source, simplifyVector=TRUE)
    if (rename)
        result <- renameFromBids(result)
    return (result)
}

#' @rdname bidsJson
#' @export
toBidsJson <- function (source, path = NULL, rename = FALSE)
{
    if (is.list(source) && rename)
        source <- renameToBids(source)
    if (is.null(path))
        return (jsonlite::toJSON(source, auto_unbox=TRUE, digits=NA, pretty=TRUE))
    else
        jsonlite::write_json(source, path, auto_unbox=TRUE, digits=NA, pretty=TRUE)
}
