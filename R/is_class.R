#' Logical Checks on Nifti Image classes
#'
#' @param x An object of class \code{niftiImage} or \code{internalImage}
#'
#' @return A logical \code{TRUE} or \code{FALSE}
#' @rdname is_class
#' @export
#' 
#' @examples
#' file = system.file("extdata", "example.nii.gz")
#' is.niftiImage(file)
#' img = readNifti(file)
#' is.niftiImage(img)
#' is.internalImage(img)
#' int_img = readNifti(file, internal = TRUE)
#' is.niftiImage(int_img)
#' is.internalImage(int_img)
is.niftiImage <- function(x) {
  inherits(x, "niftiImage")
}

#' @rdname is_class
#' @export
is.internalImage <- function(x) {
  inherits(x, "internalImage")
}
