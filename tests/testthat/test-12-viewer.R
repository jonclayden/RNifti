context("Displaying images")

test_that("greyscale and RGB images can be viewed", {
    expect_null(view(system.file("extdata", "example.nii.gz", package="RNifti"), interactive=FALSE))
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    expect_message(view(lyr(image,max=NA), interactive=FALSE), "(0, 748)")
    rgbImage <- readNifti(system.file("extdata", "example_rgb.nii.gz", package="RNifti"))
    expect_equal(rgbImage$datatype, 128L)
    expect_null(view(rgbImage, interactive=FALSE))
    grDevices::dev.off()
})
