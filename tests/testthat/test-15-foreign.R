context("Working with objects from other packages")

test_that("image objects from oro.nifti can be read", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    
    skip_if_not_installed("oro.nifti")
    
    # The oro.nifti package warns about nonzero slope, which is nothing to worry about
    image <- suppressWarnings(oro.nifti::readNIfTI(imagePath))
    expect_equal(niftiHeader(image)$bitpix, 32L)
})

test_that("image objects from tractor.base can be read", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    
    skip_if_not_installed("tractor.base")
    
    reportr::setOutputLevel(reportr::OL$Warning)
    
    # NB: the $ operator shortcut can't be used since "image" isn't a niftiImage
    image <- tractor.base::readImageFile(imagePath)
    expect_equal(niftiHeader(image)$bitpix, 32L)
    
    image <- tractor.base::readImageFile(imagePath, sparse=TRUE)
    expect_equal(niftiHeader(image)$bitpix, 32L)
})
