#' @title Convert niftiImage to oro.nifti nifti object
#' @description Converts a niftiImage from RNifti to a 
#' nifti object from the oro.nifti package
#' @param image nifitImage object
#'
#' @return Object of class \code{\link{nifti}}
#' @export
#' @importFrom oro.nifti nifti
#' @importFrom methods as slot<- slotNames
nii2oro = function(image) {
  arr = as(image, "array")
  hdr = dumpNifti(image)
  
  hdr_names = sort(names(hdr))
  
  stopifnot( !("dim_" %in% names(hdr)))
  hdr$dim_ = hdr$dim
  hdr$dim = NULL
  
  img = nifti(arr)
  
  #######################################
  # Get the names
  #######################################
  img_names = sort(slotNames(img))
  
  #######################################
  # Stop if these have changed
  # This means one of the standards have changed
  #######################################
  wrong = setdiff(hdr_names, img_names)
  stopifnot(length(wrong) == 0)
  
  #######################################
  # See where they both work
  #######################################
  both = intersect(img_names, hdr_names)
  
  
  #######################################
  # 
  #######################################
  # add = setdiff(img_names, hdr_names)
  
  hdr$dim_info = ""
  for (islot in both) {
    slot(img, islot) = hdr[[islot]]
  }
  
  return(img)
}
