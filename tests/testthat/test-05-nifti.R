context("Reading and writing NIfTI files")

test_that("NIfTI files can be read and written", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    tempPath <- tempfile()
    
    expect_that(dim(readNifti(imagePath,internal=FALSE)), equals(c(96L,96L,60L)))
    expect_that(dim(readNifti(imagePath,internal=TRUE)), equals(c(96L,96L,60L)))
    expect_that(print(readNifti(imagePath,internal=TRUE)), prints_text("2.5 mm per voxel",fixed=TRUE))
    
    image <- readNifti(imagePath)
    expect_that(pixunits(image), equals(c("mm","s")))
    
    writeNifti(image, tempPath)
    expect_that(pixdim(readNifti(tempPath)), equals(c(2.5,2.5,2.5)))
    unlink(tempPath)
    
    writeNifti(image, tempPath, datatype="short")
    expect_that(dumpNifti(tempPath)$bitpix, equals(16L))
    unlink(tempPath)
    
    expect_that(dumpNifti(image)$dim, equals(c(3L,96L,96L,60L,1L,1L,1L,1L)))
    pixdim(image) <- c(5,5,5)
    expect_that(dumpNifti(image)$pixdim, equals(c(-1,5,5,5,0,0,0,0)))
    pixunits(image) <- c("m","ms")
    expect_that(dumpNifti(image)$xyzt_units, equals(17L))
    
    image <- updateNifti(image, list(intent_code=1000L))
    expect_that(dumpNifti(image)$intent_code, equals(1000L))
})
