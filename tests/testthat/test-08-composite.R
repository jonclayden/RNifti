context("Handling composite types")

test_that("complex datatypes are handled properly", {
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    data <- sqrt(array(as.complex(image)-468, dim=dim(image)))
    complexImage <- updateNifti(data, image)
    
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
    
    set.seed(1)
    k <- stats::kmeans(as.vector(image), 3L)
    data <- rgbArray(k$cluster==1, k$cluster==2, k$cluster==3, dim=dim(image))
    rgbImage <- updateNifti(data, image)
    
    expect_s3_class(data, "rgbArray")
    expect_equal(rgbImage$datatype, 128L)
})
