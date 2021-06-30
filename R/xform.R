.isXformMatrix <- function (x) return (is.numeric(x) && identical(dim(x),c(4L,4L)))

#' Obtain or replace the ``xform'' transforms for an image
#' 
#' These functions convert the ``qform'' or ``sform'' information in a NIfTI
#' header to or from a corresponding affine matrix. These two ``xform''
#' mechanisms are defined by the NIfTI standard, and may both be in use in a
#' particular image header. They define the relationship between the storage
#' order of the image and real space.
#' 
#' Image orientation is indicated using a three-character string, with each
#' character indicating the approximate world-space direction of the positive
#' axes in the first, second and third dimensions, in order. Each character may
#' be `R' for left-to-right, `L' for right-to-left, `A' for posterior-to-
#' anterior, `P' for anterior-to-posterior, `S' for inferior-to-superior, or
#' `I' for superior-to-inferior. The default for NIfTI is RAS, meaning that the
#' first dimension points towards the right, the second towards the front and
#' the third towards the top. An xform matrix is an affine transform relative
#' to that default.
#' 
#' The upper-left 3x3 matrix in a 3D affine transform governs scale, rotation
#' and skew, while the last column is a translation. (The \code{rotation}
#' function extracts the rotation part alone.) The final row is always
#' (0,0,0,1). Reorienting an image involves permuting and possibly reversing
#' some of the axes, both in the data and the metadata. The sense of the
#' translation may also need to be reversed, but this is only possible if the
#' image dimensions are known.
#' 
#' @param image,x An image, in any acceptable form (see \code{\link{asNifti}}),
#'   or a 4x4 numeric xform matrix.
#' @param useQuaternionFirst A single logical value. If \code{TRUE}, the
#'   ``qform'' matrix will be used first, if it is defined; otherwise the
#'   ``sform'' matrix will take priority.
#' @param value A new 4x4 qform or sform matrix, or orientation string. If a
#'   matrix has a \code{"code"} attribute, the appropriate qform or sform code
#'   is also set.
#' @return For \code{xform}, an affine matrix corresponding to the ``qform''
#'   or ``sform'' information in the image header, with an \code{"imagedim"}
#'   attribute giving the original image dimensions and a \code{"code"}
#'   attribute giving the corresponding xform code. For \code{orientation}, a
#'   string with three characters indicating the (approximate) orientation of
#'   the image. The replacement forms return the modified object.
#' 
#' @note The qform and sform replacement functions are for advanced users only.
#'   Modifying the transforms without knowing what you're doing is usually
#'   unwise, as you can make the image object inconsistent.
#' 
#' @examples
#' im <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' xform(im)
#' 
#' # Remove the qform information
#' qform(im) <- structure(diag(4), code=0L)
#' 
#' # The same as above, since the sform is unmodified
#' xform(im)
#' 
#' # The identity matrix corresponds to RAS orientation
#' orientation(diag(4))
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{https://www.nitrc.org/docman/view.php/26/64/nifti1.h})
#'   is the definitive reference on ``xform'' conventions.
#' @export
xform <- function (image, useQuaternionFirst = TRUE)
{
    if (.isXformMatrix(image))
        return (image)
    else
        return (.Call("getXform", asNifti(image,internal=TRUE), isTRUE(useQuaternionFirst), PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
"qform<-" <- function (x, value)
{
    if (!inherits(x, "niftiHeader"))
        x <- asNifti(x)
    
    return (.Call("setXform", x, value, TRUE, PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
"sform<-" <- function (x, value)
{
    if (!inherits(x, "niftiHeader"))
        x <- asNifti(x)
    
    return (.Call("setXform", x, value, FALSE, PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
orientation <- function (x, useQuaternionFirst = TRUE)
{
    if (!.isXformMatrix(x))
        x <- asNifti(x, internal=TRUE)
    
    return (.Call("getOrientation", x, isTRUE(useQuaternionFirst), PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
"orientation<-" <- function (x, value)
{
    if (!.isXformMatrix(x))
        x <- asNifti(x)
    
    return (.Call("setOrientation", x, as.character(value), PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
rotation <- function (x, useQuaternionFirst = TRUE)
{
    if (!.isXformMatrix(x))
        x <- asNifti(x, internal=TRUE)
    
    return (.Call("getRotation", x, isTRUE(useQuaternionFirst), PACKAGE="RNifti"))
}

#' Transform points between voxel and ``world'' coordinates
#' 
#' These functions are used to transform points from dimensionless pixel or
#' voxel coordinates to ``real-world'' coordinates, typically in millimetres,
#' and back. Actual pixel units can be obtained using the
#' \code{\link{pixunits}} function. The \code{origin} function gives the voxel
#' coordinates of the real-world origin.
#' 
#' @param points A vector giving the coordinates of a point, or a matrix with
#'   one point per row.
#' @param image The image in whose space the points are given, or a 4x4 numeric
#'   xform matrix.
#' @param simple A logical value: if \code{TRUE} then the transformation is
#'   performed simply by rescaling the points according to the voxel dimensions
#'   recorded in the \code{image}. Otherwise the full xform matrix is used.
#' @param ... Additional arguments to \code{\link{xform}}.
#' @return A vector or matrix of transformed points.
#' 
#' @note Voxel coordinates are assumed by these functions to use R's indexing
#'   convention, beginning from 1.
#' 
#' @examples
#' im <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
#' 
#' # Find the origin
#' origin(im)
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{xform}}, \code{\link{pixdim}}, \code{\link{pixunits}}
#' @export
voxelToWorld <- function (points, image, simple = FALSE, ...)
{
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- pixdim(image)[seq_len(ncol(points))]
        return (drop(t(apply(points-1, 1, function(x) x*abs(voxelDims)))))
    }
    else
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        affine <- xform(image, ...)
        
        nDims <- ncol(points)
        if (nDims != 2 && nDims != 3)
            stop("Points must be two- or three-dimensional")
        
        if (nDims == 2)
            affine <- matrix(affine[c(1,2,4,5,6,8,13,14,16)], ncol=3, nrow=3)
        
        points <- cbind(points-1, 1)
        newPoints <- affine %*% t(points)
        return (drop(t(newPoints[1:nDims,,drop=FALSE])))
    }
}

#' @rdname voxelToWorld
#' @export
worldToVoxel <- function (points, image, simple = FALSE, ...)
{
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- pixdim(image)[seq_len(ncol(points))]
        return (drop(t(apply(points, 1, function(x) x/abs(voxelDims)) + 1)))
    }
    else
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        affine <- solve(xform(image, ...))
        
        nDims <- ncol(points)
        if (nDims != 2 && nDims != 3)
            stop("Points must be two- or three-dimensional")
        
        if (nDims == 2)
            affine <- matrix(affine[c(1,2,4,5,6,8,13,14,16)], ncol=3, nrow=3)
        
        points <- cbind(points, 1)
        newPoints <- affine %*% t(points) + 1
        return (drop(t(newPoints[1:nDims,,drop=FALSE])))
    }
}

#' @rdname voxelToWorld
#' @export
origin <- function (image, ...)
{
    worldToVoxel(c(0,0,0), image, ...)
}
