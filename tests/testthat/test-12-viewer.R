test_that("greyscale and RGB images can be viewed", {
    expect_null(view(system.file("extdata", "example.nii.gz", package="RNifti"), interactive=FALSE))
    
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    expect_message(view(lyr(image,max=NA), interactive=FALSE), "(0, 748)")
    expect_null(view(image[,,34], interactive=FALSE))
    
    rgbImage <- readNifti(system.file("extdata", "example_rgb.nii.gz", package="RNifti"))
    expect_equal(rgbImage$datatype, 128L)
    expect_null(view(rgbImage, interactive=FALSE))
    expect_null(view(lyr(rgbImage, mask=image<400), interactive=FALSE))
    
    # NB: This conversion is flawed, since it does not restore sign, but it's fine for a test
    vecImage <- structure(channels(rgbImage, c("red","green","blue"))/255, dim=c(dim(rgbImage),1L,3L))
    vecImage <- updateNifti(vecImage, rgbImage)
    vecImage$intent_code <- 1007
    expect_null(view(image, vecImage, interactive=FALSE))
    
    volImage <- readNifti(system.file("extdata", "example_4d.nii.gz", package="RNifti"))
    expect_null(view(lyr(volImage,max=3), point=c(48,48,20), infoPanel=timeSeriesPanel, interactive=FALSE))
    
    grDevices::dev.off()
})
