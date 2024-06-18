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
