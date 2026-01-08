# Create a small 2D "image" with random data (and 1 NA) to test for equivalence of summary operations
data <- matrix(rnorm(9), 3, 3)
data[sample(1:9,1)] <- NA

# Standard niftiImage cases should be pretty simple because these are arrays
image <- asNifti(data)
expect_equivalent(data, unclass(image))
expect_equivalent(data+1, unclass(image+1))
expect_equivalent(abs(data), unclass(abs(image)))
expect_equal(mean(image), NA_real_)
expect_equal(sum(image), NA_real_)
expect_equal(sum(is.na(image)), 1)
expect_equal(mean(image,na.rm=TRUE), mean(data,na.rm=TRUE))
expect_equal(range(image,na.rm=TRUE), range(data,na.rm=TRUE))

# Internal images are fundamentally character mode, so these tests rely on specialised methods
image <- asNifti(data, internal=TRUE)
expect_equivalent(data+1, unclass(image+1))
expect_equivalent(abs(data), unclass(abs(image)))
expect_equal(mean(image), NA_real_)
expect_equal(sum(image), NA_real_)
expect_equal(sum(is.na(image)), 1)
expect_equal(mean(image,na.rm=TRUE), mean(data,na.rm=TRUE))
expect_equal(range(image,na.rm=TRUE), range(data,na.rm=TRUE))
expect_equal(range(image+1,data,na.rm=TRUE), c(min(data,na.rm=TRUE), max(data,na.rm=TRUE)+1))

# Similar with integer-mode data (since NA_integer_ is not reflected naturally in C/C++)
data <- matrix(sample(-5:15,9), 3, 3)
data[sample(1:9,1)] <- NA_integer_
image <- asNifti(data, internal=TRUE)
expect_equal(sum(image), NA_integer_)
expect_equal(sum(image,na.rm=TRUE), sum(data,na.rm=TRUE))
expect_equal(range(image,na.rm=TRUE), range(data,na.rm=TRUE))
