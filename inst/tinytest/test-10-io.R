# Reading and writing NIfTI-1 format, and considering its interpretation as ANALYZE-7.5
imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
volumeImagePath <- system.file("extdata", "example_4d.nii.gz", package="RNifti")
tempPath <- paste(tempfile(), "nii.gz", sep=".")

expect_equal(niftiVersion(imagePath), structure(1L,names=imagePath))
expect_equal(dim(readNifti(imagePath,internal=FALSE)), c(96L,96L,60L))
expect_equal(dim(readNifti(imagePath,internal=TRUE)), c(96L,96L,60L))
expect_stdout(print(readNifti(imagePath,internal=TRUE)), "2.5 mm per voxel")

# Non-NIfTI and missing files should return -1
expect_warning(version <- niftiVersion(system.file("COPYRIGHTS",package="RNifti")))
expect_equivalent(version, -1L)
expect_warning(version <- niftiVersion(system.file("no-such-file",package="RNifti")))
expect_equivalent(version, -1L)

image <- readNifti(imagePath)
expect_inherits(image, "niftiImage")
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

expect_stdout(print(niftiHeader(image)), "NIfTI-1 header")
expect_equal(image$bitpix, 32L)
expect_null(image$glmax)
expect_equal(niftiHeader(imagePath,unused=TRUE)$glmax, 0L)
expect_equal(analyzeHeader(imagePath)$datatype, 4L)
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
expect_inherits(image, "internalImage")
expect_inherits(image, "niftiImage")
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

# Direct access to the NIfTI data blob
imageHeader <- niftiHeader(imagePath)
expect_equal(RNifti:::readBlob(imagePath, 1, imageHeader$datatype, imageHeader$vox_offset + imageHeader$bitpix * 271047 / 8), 368)

# Check that internal indexing still works when data scaling is in play
compressedInternalImage <- readNifti(tempPath, internal=TRUE)
expect_equal(round(compressedInternalImage[40,40,30]), 363)
expect_equal(compressedInternalImage[271048], compressedImage[271048])
expect_equal(compressedInternalImage[40,,30], compressedImage[40,,30])
expect_equal(compressedInternalImage[c(40,42),,c(30,32)], compressedImage[c(40,42),,c(30,32)])

image <- readNifti(volumeImagePath, volumes=1:2)
expect_equal(dim(image), c(96L,96L,60L,2L))
expect_equal(max(image), 2)
expect_stdout(print(image), "x 1 s")    # time units only appear for 4D+ images

analyze <- analyzeHeader()
expect_inherits(analyze, "analyzeHeader")
expect_stdout(print(analyze), "ANALYZE-7.5")
expect_equal(analyze$regular, "r")

# Reading and writing NIfTI-2 and ANALYZE-7.5 formats
imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
tempStem <- tempfile()
tempPath <- paste(tempStem, "nii.gz", sep=".")
analyzePaths <- c(header=paste(tempStem,"hdr",sep="."), image=paste(tempStem,"img",sep="."))

image <- readNifti(imagePath, internal=TRUE)
paths <- writeNifti(image, tempPath, version=2)
expect_equal(paths, c(header=tempPath,image=tempPath))
expect_equal(niftiVersion(tempPath), structure(2L,names=tempPath))
nifti2Image <- readNifti(tempPath, internal=TRUE)
expect_stdout(print(niftiHeader(nifti2Image)), "NIfTI-2 header")
expect_equal(nifti2Image$vox_offset, 544L)
expect_equal(image[40,40,30], nifti2Image[40,40,30])
unlink(tempPath)

paths <- writeAnalyze(image, analyzePaths[1])
expect_equal(paths, analyzePaths)
expect_equal(niftiVersion(tempPath), structure(0L,names=tempPath))
analyzeImage <- readAnalyze(tempPath, internal=TRUE)
expect_equal(image[40,40,30], analyzeImage[40,40,30])
