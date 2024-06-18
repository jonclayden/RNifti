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
