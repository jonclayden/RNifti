#' Obtain or replace the ``xform'' transforms for an image
#' 
#' These functions convert the ``qform'' or ``sform'' information in a NIfTI
#' header to or from a corresponding affine matrix. These two ``xform''
#' mechanisms are defined by the NIfTI standard, and may both be in use in a
#' particular image header.
#' 
#' @param image,x An image, in any acceptable form (see
#'   \code{\link{retrieveNifti}}).
#' @param useQuaternionFirst A single logical value. If \code{TRUE}, the
#'   ``qform'' matrix will be used first, if it is defined; otherwise the
#'   ``sform'' matrix will take priority.
#' @param value A new 4x4 qform or sform matrix. If the matrix has a
#'   \code{"code"} attribute, the appropriate qform or sform code is also set.
#' @return A affine matrix corresponding to the ``qform'' or ``sform''
#'   information in the image header. This is a plain matrix, which does not
#'   have the \code{"affine"} class or \code{source} and \code{target}
#'   attributes.
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
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1})
#'   is the definitive reference on ``xform'' conventions.
#' @export
xform <- function (image, useQuaternionFirst = TRUE)
{
    return (.Call("getXform", image, isTRUE(useQuaternionFirst), PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
"qform<-" <- function (x, value)
{
    return (.Call("setXform", x, value, TRUE, PACKAGE="RNifti"))
}

#' @rdname xform
#' @export
"sform<-" <- function (x, value)
{
    return (.Call("setXform", x, value, FALSE, PACKAGE="RNifti"))
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
#' @param image The image in whose space the points are given.
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
