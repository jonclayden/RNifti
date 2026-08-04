imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")

if (requireNamespace("oro.nifti", quietly=TRUE))
{
    # The oro.nifti package warns about nonzero slope, which is nothing to worry about
    image <- suppressWarnings(oro.nifti::readNIfTI(imagePath))
    expect_equal(niftiHeader(image)$bitpix, 32L)
    expect_equal(asNifti(image)[40,40,30], 368)
}
    
if (requireNamespace("tractor.base", quietly=TRUE))
{
    reportr::setOutputLevel(reportr::OL$Warning)
    
    # NB: the $ operator shortcut can't be used since "image" isn't a niftiImage
    image <- tractor.base::readImageFile(imagePath)
    expect_equal(niftiHeader(image)$bitpix, 32L)
    expect_equal(asNifti(image)[40,40,30], 368)
    
    image <- tractor.base::readImageFile(imagePath, sparse=TRUE)
    expect_equal(niftiHeader(image)$bitpix, 32L)
    expect_equal(asNifti(image)[40,40,30], 368)
}

# A package extending RNifti to a class it knows nothing about natively, purely by providing an
# asNifti() S3 method, should still work correctly through any function that ultimately constructs
# a NiftiImage in C++
referenceImage <- readNifti(imagePath)

asNifti.customTestImage <- function (x, ...) asNifti(x$wrapped, ...)
registerS3method("asNifti", "customTestImage", asNifti.customTestImage)

# The underlying representation is a list with an unrelated field name, so if the custom method
# were bypassed, generic list handling would silently misinterpret it rather than erroring
customImage <- structure(list(wrapped=referenceImage), class="customTestImage")
expect_equal(orientation(customImage), orientation(referenceImage))
expect_equal(xform(customImage), xform(referenceImage))
expect_equal(niftiHeader(customImage)$bitpix, niftiHeader(referenceImage)$bitpix)
expect_equal(asNifti(customImage)[40,40,30], asNifti(referenceImage)[40,40,30])

# A classed object with no asNifti() method at all is unaffected: falls back to old behaviour,
# i.e., an error if it isn't otherwise recognisable (not list- or array-shaped)...
expect_error(orientation(structure(1:5, class="noAsNiftiMethodForThis")))

# ...but no error, just a meaningless default result, if it happens to be list-shaped
expect_equal(orientation(structure(list(a=1), class="noAsNiftiMethodForThis")), "RAS")

# A method that fails to actually convert its input should not cause infinite recursion
asNifti.brokenTestImage <- function (x, ...) x
registerS3method("asNifti", "brokenTestImage", asNifti.brokenTestImage)
brokenImage <- structure(list(a=1), class="brokenTestImage")
expect_equal(orientation(brokenImage), "RAS")
