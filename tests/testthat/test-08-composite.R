test_that("complex datatypes are handled properly", {
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    data <- sqrt(array(as.complex(image)-468, dim=dim(image)))
    complexImage <- asNifti(data, image)
    
    expect_equal(complexImage[40,40,30], 0+10i)
    expect_equal(complexImage$datatype, 1792L)
    expect_equal(complexImage$bitpix, 128L)
    expect_equivalent(data, unclass(as.array(complexImage)))
    expect_output(print(complexImage), "complex")
    
    writeNifti(complexImage, tempPath)
    header <- niftiHeader(tempPath)
    expect_equal(header$datatype, 1792L)
    expect_equal(header$bitpix, 128L)

    complexImage <- readNifti(tempPath, internal=TRUE)
    expect_equal(complexImage[40,40,30], 0+10i)
})

test_that("RGB datatypes are handled properly", {
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    expect_equal(as.character(rgbArray(diag(3))), c("#FF0000","#00FF00","#0000FF"))
    expect_equivalent(channels(diag(3),"red")[,1], c(255L,0L,0L))
    
    k <- stats::kmeans(as.vector(image), 3L)
    data <- rgbArray(k$cluster==1, k$cluster==2, k$cluster==3, dim=dim(image))
    rgbImage <- asNifti(data, image)
    cluster <- k$cluster[40 + 39*96 + 29*9216]
    refValue <- rgbArray(diag(3))[cluster]
    
    expect_s3_class(data, "rgbArray")
    expect_equal(rgbImage[40,40,30], refValue)
    expect_equal(rgbImage$datatype, 128L)
    expect_equal(rgbImage$bitpix, 24L)
    expect_equivalent(unclass(data), unclass(as.array(rgbImage)))
    expect_equivalent(channels(rgbImage,c("red","green","blue")[cluster])[40,40,30,1], 255L)
    
    writeNifti(rgbImage, tempPath)
    header <- niftiHeader(tempPath)
    expect_equal(header$datatype, 128L)
    expect_equal(header$bitpix, 24L)
    
    complexImage <- readNifti(tempPath, internal=TRUE)
    expect_equal(complexImage[40,40,30], refValue)
})
