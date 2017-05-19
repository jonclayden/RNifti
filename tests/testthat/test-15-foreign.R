context("Working with objects from other packages")

test_that("image objects from other packages can be read", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    
    if (system.file(package="oro.nifti") == "")
        skip("The \"oro.nifti\" package is not available")
    else
    {
        # The oro.nifti package warns about nonzero slope, which is nothing to worry about
        image <- suppressWarnings(oro.nifti::readNIfTI(imagePath))
        expect_that(dumpNifti(image)$bitpix, equals(32L))
    }
    
    if (system.file(package="tractor.base") == "")
        skip("The \"tractor.base\" package is not available")
    else
    {
        reportr::setOutputLevel(reportr::OL$Warning)
        
        image <- tractor.base::readImageFile(imagePath)
        expect_that(dumpNifti(image)$bitpix, equals(32L))
        
        image <- tractor.base::readImageFile(imagePath, sparse=TRUE)
        expect_that(dumpNifti(image)$bitpix, equals(32L))
    }
})