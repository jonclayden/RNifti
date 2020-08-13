context("Reading, writing and manipulating NIfTI objects")

test_that("NIfTI objects can be created from data", {
    data <- array(rnorm(24), dim=c(3L,2L,4L))
    image <- asNifti(data)
    
    expect_equal(ndim(data), 3L)
    expect_s3_class(image, "niftiImage")
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
    data <- asNifti(data, image)
    expect_equal(pixdim(data), c(2.5,2.5,2.5))
    
    paths <- writeNifti(image, tempPath)
    expect_equal(paths, c(header=tempPath,image=tempPath))
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
    array <- as.array(image)
    expect_s3_class(image, "internalImage")
    expect_s3_class(image, "niftiImage")
    expect_equal(image[], array)
    expect_equal(image[,,], array)
    expect_equal(array[40,40,30], 368)
    expect_equal(image[40,40,30], 368)
    expect_equal(image[271048], 368)
    expect_identical(image[552961], NA_integer_)
    expect_equal(array[40,,30], image[40,,30])
    expect_equal(array[40:42,,30:32], image[40:42,,30:32])
    m <- matrix(c(40,40,30,42,30,32), byrow=TRUE, ncol=3)
    expect_equal(array[m], image[m])
    expect_equal(array[c(40,42),,c(30,32)], image[c(40,42),,c(30,32)])
    expect_error(dim(image) <- c(60L,96L,96L))
    expect_error(pixdim(image) <- c(1,1,1))
    expect_error(image[40,40,30] <- 400)
    
    # Check that internal indexing still works when data scaling is in play
    compressedInternalImage <- readNifti(tempPath, internal=TRUE)
    expect_equal(round(compressedInternalImage[40,40,30]), 363)
    expect_equal(compressedInternalImage[271048], compressedImage[271048])
    expect_equal(compressedInternalImage[40,,30], compressedImage[40,,30])
    expect_equal(compressedInternalImage[c(40,42),,c(30,32)], compressedImage[c(40,42),,c(30,32)])
    
    image <- readNifti(volumeImagePath, volumes=1:2)
    expect_equal(dim(image), c(96L,96L,60L,2L))
    expect_equal(max(image), 2)
    expect_output(print(image), "x 1 s")    # time units only appear for 4D+ images
    
    analyze <- analyzeHeader()
    expect_s3_class(analyze, "analyzeHeader")
    expect_output(print(analyze), "ANALYZE-7.5")
    expect_equal(analyze$regular, "r")
})

test_that("NIfTI-2 and ANALYZE-7.5 format files can be written and read", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    tempStem <- tempfile()
    tempPath <- paste(tempStem, "nii.gz", sep=".")
    analyzePaths <- c(header=paste(tempStem,"hdr",sep="."), image=paste(tempStem,"img",sep="."))
    
    image <- readNifti(imagePath, internal=TRUE)
    paths <- writeNifti(image, tempPath, version=2)
    expect_equal(paths, c(header=tempPath,image=tempPath))
    expect_equal(niftiVersion(tempPath), structure(2L,names=tempPath))
    nifti2Image <- readNifti(tempPath, internal=TRUE)
    expect_equal(nifti2Image$vox_offset, 544L)
    expect_equal(image[40,40,30], nifti2Image[40,40,30])
    unlink(tempPath)
    
    paths <- writeAnalyze(image, analyzePaths[1])
    expect_equal(paths, analyzePaths)
    expect_equal(niftiVersion(tempPath), structure(0L,names=tempPath))
    analyzeImage <- readAnalyze(tempPath, internal=TRUE)
    expect_equal(image[40,40,30], analyzeImage[40,40,30])
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
    image <- asNifti(image, datatype="float")
    expect_equal(image$datatype, 16L)
    
    # Empty values are ignored (with a warning)
    expect_warning(image$intent_code <- integer(0))
    
    # Only the first value will be used (with a warning)
    expect_warning(image$intent_code <- 1002:1005)
    expect_equal(image$intent_code, 1002L)
    
    # Vector-valued fields must have the right length
    expect_error(image$srow_x <- 1.0)
    image$srow_x[1] <- 1.0
    expect_equal(image$srow_x[1], 1.0)
    
    image <- readNifti(imagePath, internal=TRUE)
    image <- RNifti:::rescaleNifti(image, c(0.5,0.5,0.5))
    expect_equal(pixdim(image), c(5,5,5))
    expect_warning(as.array(image), "no data")
})

test_that("NIfTI objects have the expected copying semantics", {
    im1 <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"), internal=TRUE)
    im2 <- im1
    # Only applies for internal images, because otherwise the data will be updated to match the R array
    expect_true(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))
    im1$intent_code <- 1000L
    expect_false(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))
})

test_that("NAs are preserved across datatypes", {
    # Original datatype is int16/short
    image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    image[40,40,30] <- NA
    writeNifti(image, tempPath, datatype="int")
    expect_equal(readNifti(tempPath)[40,40,30], NA_integer_)
    writeNifti(image, tempPath, datatype="double")
    image <- readNifti(tempPath)
    expect_equal(image[40,40,30], NA_real_)
    image[42,42,32] <- NA
    writeNifti(image, tempPath, datatype="int")
    expect_equal(readNifti(tempPath)[42,42,32], NA_integer_)
    writeNifti(image, tempPath, datatype="short")
    expect_equal(readNifti(tempPath)[42,42,32], 0)
})
