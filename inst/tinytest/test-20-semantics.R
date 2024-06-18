# Check for expected copying semantics
im1 <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"), internal=TRUE)
im2 <- im1

# Only applies for internal images, because otherwise the data will be updated to match the R array
expect_true(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))
im1$intent_code <- 1000L
expect_false(all(RNifti:::addresses(im1) == RNifti:::addresses(im2)))

# Unwrapping and wrapping NiftiImage pointers involves copies
im3 <- RNifti:::wrapPointer(RNifti:::unwrapPointer(im1))
expect_equal(niftiHeader(im1), niftiHeader(im3))
expect_inherits(im3, "niftiImage")
expect_false(all(RNifti:::addresses(im1) == RNifti:::addresses(im3)))

# Check that NAs are preserved across datatypes (original datatype is int16/short)
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
