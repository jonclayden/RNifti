test_that("NIfTI sform/qform operations work", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath)
    
    expect_equal(diag(xform(image)), c(-2.5,2.5,2.5,1))
    expect_equal(diag(xform(image,useQuaternionFirst=FALSE)), c(-2.5,2.5,2.5,1))
    expect_equal(round(origin(image)), c(50,39,23))
    expect_equal(rotation(image), diag(c(-1,1,1)))
    
    xform <- xform(image)
    expect_equal(round(origin(xform)), c(50,39,23))
    expect_equal(rotation(xform), diag(c(-1,1,1)))
    expect_equal(orientation(xform), "LAS")
    rotationPIR <- matrix(0,3,3)
    rotationPIR[c(2,6,7)] <- c(-1,-1,1)
    orientation(xform) <- "PIR"
    expect_equal(rotation(xform), rotationPIR)
    expect_equal(orientation(xform), "PIR")
    
    point <- c(40, 40, 20)
    expect_equal(round(voxelToWorld(point,image)), c(25,2,-8))
    expect_equal(voxelToWorld(point,image,simple=TRUE), c(97.5,97.5,47.5))
    expect_equal(worldToVoxel(voxelToWorld(point,image),image), point)
    expect_equal(worldToVoxel(voxelToWorld(point,image,simple=TRUE),image,simple=TRUE), point)
    
    # NB: indexing into the image currently discards the metadata, including the xforms
    expect_equal(voxelToWorld(point[1:2],image[,,20L]), c(39,39))
    
    # Create a rotation-only affine: 5 degrees about the x-axis
    theta <- 5 / 180 * pi
    affine <- diag(4)
    affine[c(6,7,10,11)] <- c(cos(theta), sin(theta), -sin(theta), cos(theta))
    
    # Copy image object to check copy semantics
    originalImage <- image
    xform <- structure(affine %*% xform(image), code=4)
    qform(image) <- xform
    sform(image) <- xform
    expect_equal(image$qform_code, 4L)
    expect_equal(image$sform_code, 4L)
    expect_equal(round(image$srow_z,1), c(0,0.2,2.5,-63.1))
    expect_equal(originalImage$qform_code, 2L)
    
    # Header-only xform replacement
    header <- niftiHeader(image)
    qform(header) <- xform
    sform(header) <- xform
    expect_equal(round(header$quatern_d,2), 0.04)
    expect_equal(round(header$srow_z,1), c(0,0.2,2.5,-63.1))
})

test_that("image data and metadata can be reoriented", {
    imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
    image <- readNifti(imagePath, internal=FALSE)
    reorientedImage <- image
    
    orientation(reorientedImage) <- "PIR"
    expect_equal(orientation(image), "LAS")
    expect_equal(orientation(reorientedImage), "PIR")
    expect_equal(origin(reorientedImage), (dim(image)-origin(image)+1)[c(2,3,1)], tolerance=1e-4)
    expect_equal(image[40,40,30], reorientedImage[57,31,57])
    orientation(reorientedImage) <- "RAS"
    expect_equal(image[40,40,30], reorientedImage[57,40,30])
    orientation(reorientedImage) <- "LAS"
    expect_equal(xform(image), xform(reorientedImage), tolerance=1e-4)
})
