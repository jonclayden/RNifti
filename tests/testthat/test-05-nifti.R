context("Reading, writing and manipulating NIfTI objects")

test_that("NIfTI objects can be created from data", {
    data <- array(rnorm(24), dim=c(3L,2L,4L))
    image <- retrieveNifti(data)
    
    expect_equal(ndim(data), 3L)
    expect_s3_class(image, "internalImage")
    expect_equal(dim(image), c(3L,2L,4L))
    expect_equal(pixdim(data), c(1,1,1))
    expect_equal(pixdim(image), c(1,1,1))
    expect_equal(pixunits(data), "Unknown")
    expect_equal(pixunits(image), "Unknown")
    expect_equal(image$datatype, 64L)
})

test_that("NIfTI files can be read and written", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    volumeImagePath <- system.file("extdata", "example_4d.nii.gz", package="RNifti")
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    expect_equal(niftiVersion(imagePath), structure(1L,names=imagePath))
    expect_equal(dim(readNifti(imagePath,internal=FALSE)), c(96L,96L,60L))
    expect_equal(dim(readNifti(imagePath,internal=TRUE)), c(96L,96L,60L))
    expect_output(print(readNifti(imagePath,internal=TRUE)), "2.5 mm per voxel")
    
    image <- readNifti(imagePath)
    expect_s3_class(image, "niftiImage")
    expect_equal(image[40,40,30], 368)
    expect_equal(pixunits(image), c("mm","s"))
    
    expect_equal(pixdim(imagePath), c(2.5,2.5,2.5))
    expect_equal(pixunits(imagePath), c("mm","s"))
    
    data <- as.vector(image)        # strips all attributes
    expect_equal(pixdim(data), 1)
    data <- updateNifti(data, image)
    expect_equal(pixdim(data), c(2.5,2.5,2.5))
    
    writeNifti(image, tempPath)
    expect_equal(pixdim(readNifti(tempPath)), c(2.5,2.5,2.5))
    unlink(tempPath)
    
    expect_output(print(niftiHeader(image)), "NIfTI-1 header")
    expect_equal(image$bitpix, 32L)
    writeNifti(image, tempPath, datatype="short")
    expect_equal(niftiHeader(tempPath)$bitpix, 16L)
    unlink(tempPath)
    
    # Type compression with index mapping
    writeNifti(image, tempPath, datatype="char")
    expect_equal(niftiHeader(tempPath)$datatype, 2L)
    compressedImage <- readNifti(tempPath)
    expect_equal(compressedImage$datatype, 64L)
    expect_equal(min(image), min(compressedImage))
    expect_equal(max(image), max(compressedImage))
    expect_equal(image, compressedImage, tolerance=0.01)
    expect_equal(round(compressedImage[40,40,30]), 363)
    
    image <- readNifti(imagePath, internal=TRUE)
    expect_s3_class(image, "internalImage")
    expect_s3_class(image, "niftiImage")
    expect_equal(as.array(image)[40,40,30], 368)
    expect_error(dim(image) <- c(60L,96L,96L))
    
    image <- readNifti(volumeImagePath, volumes=1:2)
    expect_equal(dim(image), c(96L,96L,60L,2L))
    expect_equal(max(image), 2)
    expect_output(print(image), "x 1 s")    # time units only appear for 4D+ images
    
    analyze <- analyzeHeader()
    expect_s3_class(analyze, "analyzeHeader")
    expect_output(print(analyze), "ANALYZE-7.5")
    expect_equal(analyze$regular, "r")
})

test_that("image objects can be manipulated", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath)
    
    expect_equal(image$dim, c(3L,96L,96L,60L,1L,1L,1L,1L))
    pixdim(image) <- c(5,5,5)
    expect_equal(image$pixdim, c(-1,5,5,5,0,0,0,0))
    pixunits(image) <- c("m","ms")
    expect_equal(image$xyzt_units, 17L)
    
    image$intent_code <- 1000L
    expect_equal(image$intent_code, 1000L)
    image <- updateNifti(image, datatype="float")
    expect_equal(image$datatype, 16L)
    
    image <- readNifti(imagePath, internal=TRUE)
    image <- RNifti:::rescaleNifti(image, c(0.5,0.5,0.5))
    expect_equal(pixdim(image), c(5,5,5))
    expect_warning(as.array(image), "no data")
})

test_that("NIfTI objects have the expected copying semantics", {
    im1 <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    im2 <- im1
    expect_equal(RNifti:::addresses(im1), RNifti:::addresses(im2))
    expect_true(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))
    im1$intent_code <- 1000L
    expect_false(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))
})
