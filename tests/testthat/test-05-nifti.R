context("Reading, writing and manipulating NIfTI objects")

test_that("NIfTI files can be read and written", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    compressedImagePath <- system.file("extdata", "example_compressed.nii.gz", package="RNifti")
    tempPath <- paste(tempfile(), "nii.gz", sep=".")
    
    expect_that(dim(readNifti(imagePath,internal=FALSE)), equals(c(96L,96L,60L)))
    expect_that(dim(readNifti(imagePath,internal=TRUE)), equals(c(96L,96L,60L)))
    expect_that(print(readNifti(imagePath,internal=TRUE)), prints_text("2.5 mm per voxel",fixed=TRUE))
    
    image <- readNifti(imagePath)
    expect_that(image[40,40,30], equals(368))
    expect_that(pixunits(image), equals(c("mm","s")))
    
    writeNifti(image, tempPath)
    expect_that(pixdim(readNifti(tempPath)), equals(c(2.5,2.5,2.5)))
    unlink(tempPath)
    
    expect_that(dumpNifti(image)$bitpix, equals(32L))
    writeNifti(image, tempPath, datatype="short")
    expect_that(dumpNifti(tempPath)$bitpix, equals(16L))
    unlink(tempPath)
    
    expect_that(dumpNifti(compressedImagePath)$datatype, equals(2L))
    compressedImage <- readNifti(compressedImagePath)
    expect_that(dumpNifti(compressedImage)$datatype, equals(64L))
    expect_that(round(compressedImage[40,40,30]), equals(363))
    
    image <- readNifti(imagePath, internal=TRUE)
    expect_that(as.array(image)[40,40,30], equals(368))
})

test_that("image objects can be manipulated", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath)
    
    expect_that(dumpNifti(image)$dim, equals(c(3L,96L,96L,60L,1L,1L,1L,1L)))
    pixdim(image) <- c(5,5,5)
    expect_that(dumpNifti(image)$pixdim, equals(c(-1,5,5,5,0,0,0,0)))
    pixunits(image) <- c("m","ms")
    expect_that(dumpNifti(image)$xyzt_units, equals(17L))
    
    image <- updateNifti(image, list(intent_code=1000L), datatype="float")
    expect_that(dumpNifti(image)$intent_code, equals(1000L))
    expect_that(dumpNifti(image)$datatype, equals(16L))
    
    image <- readNifti(imagePath, internal=TRUE)
    image <- RNifti:::rescaleNifti(image, c(0.5,0.5,0.5))
    expect_that(pixdim(image), equals(c(5,5,5)))
    expect_that(as.array(image), gives_warning("no data"))
})
