context("Reading, writing and manipulating NIfTI objects")

test_that("NIfTI objects can be created from data", {
    data <- array(rnorm(24), dim=c(3L,2L,4L))
    image <- retrieveNifti(data)
    
    expect_equal(ndim(data), 3L)
    expect_is(image, "internalImage")
    expect_equal(dim(image), c(3L,2L,4L))
    expect_equal(pixdim(data), c(1,1,1))
    expect_equal(pixdim(image), c(1,1,1))
    expect_equal(pixunits(data), "Unknown")
    expect_equal(pixunits(image), "Unknown")
    expect_equal(dumpNifti(image)$datatype, 64L)
})

test_that("NIfTI files can be read and written", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    volumeImagePath <- system.file("extdata", "example_4d.nii.gz", package="RNifti")
    compressedImagePath <- system.file("extdata", "example_compressed.nii.gz", package="RNifti")
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    expect_equal(niftiVersion(imagePath), structure(1L,names=imagePath))
    expect_equal(dim(readNifti(imagePath,internal=FALSE)), c(96L,96L,60L))
    expect_equal(dim(readNifti(imagePath,internal=TRUE)), c(96L,96L,60L))
    expect_output(print(readNifti(imagePath,internal=TRUE)), "2.5 mm per voxel")
    
    image <- readNifti(imagePath)
    expect_equal(image[40,40,30], 368)
    expect_equal(pixunits(image), c("mm","s"))
    
    writeNifti(image, tempPath)
    expect_equal(pixdim(readNifti(tempPath)), c(2.5,2.5,2.5))
    unlink(tempPath)
    
    expect_output(print(dumpNifti(image)), "NIfTI-1 header")
    expect_equal(dumpNifti(image)$bitpix, 32L)
    writeNifti(image, tempPath, datatype="short")
    expect_equal(dumpNifti(tempPath)$bitpix, 16L)
    unlink(tempPath)
    
    expect_equal(dumpNifti(compressedImagePath)$datatype, 2L)
    compressedImage <- readNifti(compressedImagePath)
    expect_equal(dumpNifti(compressedImage)$datatype, 64L)
    expect_equal(round(compressedImage[40,40,30]), 363)
    
    image <- readNifti(imagePath, internal=TRUE)
    expect_equal(as.array(image)[40,40,30], 368)
    expect_error(dim(image) <- c(60L,96L,96L))
    
    image <- readNifti(volumeImagePath, volumes=1:2)
    expect_equal(dim(image), c(96L,96L,60L,2L))
    expect_equal(max(image), 2)
    expect_output(print(image), "x 1 s")    # time units only appear for 4D+ images
})

test_that("image objects can be manipulated", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath)
    
    expect_equal(dumpNifti(image)$dim, c(3L,96L,96L,60L,1L,1L,1L,1L))
    pixdim(image) <- c(5,5,5)
    expect_equal(dumpNifti(image)$pixdim, c(-1,5,5,5,0,0,0,0))
    pixunits(image) <- c("m","ms")
    expect_equal(dumpNifti(image)$xyzt_units, 17L)
    
    image <- updateNifti(image, list(intent_code=1000L))
    expect_equal(dumpNifti(image)$intent_code, 1000L)
    image <- updateNifti(image, datatype="float")
    expect_equal(dumpNifti(image)$datatype, 16L)
    
    image <- readNifti(imagePath, internal=TRUE)
    image <- RNifti:::rescaleNifti(image, c(0.5,0.5,0.5))
    expect_equal(pixdim(image), c(5,5,5))
    expect_warning(as.array(image), "no data")
})
