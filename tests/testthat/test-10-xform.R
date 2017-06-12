context("NIfTI sform/qform operations")

test_that("NIfTI sform/qform operations work", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath)
    
    expect_that(diag(xform(image)), equals(c(-2.5,2.5,2.5,1)))
    expect_that(diag(xform(image,useQuaternionFirst=FALSE)), equals(c(-2.5,2.5,2.5,1)))
    
    point <- c(40, 40, 20)
    expect_that(round(voxelToWorld(point,image)), equals(c(25,2,-8)))
    expect_that(voxelToWorld(point,image,simple=TRUE), equals(c(97.5,97.5,47.5)))
    expect_that(worldToVoxel(voxelToWorld(point,image),image), equals(point))
    expect_that(worldToVoxel(voxelToWorld(point,image,simple=TRUE),image,simple=TRUE), equals(point))
    
    # NB: indexing into the image currently discards the metadata, including the xforms
    expect_that(voxelToWorld(point[1:2],image[,,20L]), equals(c(39,39)))
    
    # Copy image object to check copy semantics
    originalImage <- image
    xform <- structure(round(xform(image)), code=4)
    qform(image) <- xform
    sform(image) <- xform
    expect_that(dumpNifti(image)$qform_code, equals(4L))
    expect_that(dumpNifti(image)$sform_code, equals(4L))
    expect_that(dumpNifti(image)$srow_x, equals(c(-2,0,0,122)))
    expect_that(dumpNifti(originalImage)$qform_code, equals(2L))
    
    image <- readNifti(imagePath, internal=FALSE)
    reorientedImage <- image
    orientation(reorientedImage) <- "PIR"
    expect_that(orientation(image), equals("LAS"))
    expect_that(orientation(reorientedImage), equals("PIR"))
    expect_that(image[40,40,30], equals(reorientedImage[57,31,57]))
})
