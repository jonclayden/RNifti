# Check that NIfTI extensions can be written and read
image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
comments <- c("A first comment", "A second comment")
rawComments <- lapply(comments, function(x) structure(c(charToRaw(x), rep(as.raw(0L), 16 - (nchar(x)+8) %% 16)), code=6L))

expect_equal(extensions(image), list())
extension(image, "comment") <- comments[1]
expect_equal(extensions(image), rawComments[1])
expect_equal(extension(image,comment), rawComments[[1]])
expect_equal(extension(image,"comment",mode="character"), comments[1])
extension(image, comment) <- comments[2]
expect_equal(extensions(image), rawComments)
expect_equal(extension(image,"comment",mode="character"), as.list(comments))

# Check persistence through a write/read cycle
image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"), internal=TRUE)
tempPath <- paste(tempfile(), "nii.gz", sep=".")

comments <- c("A first comment", "A second comment")
rawComments <- lapply(comments, function(x) structure(c(charToRaw(x), rep(as.raw(0L), 16 - (nchar(x)+8) %% 16)), code=6L))
rawDouble <- writeBin(1000, raw())
extensions <- c(rawComments, list(structure(rawDouble, code=20L)))

extensions(image) <- rawComments
extension(image, "b_value") <- 1000
writeNifti(image, tempPath)
rereadImage <- readNifti(tempPath)
expect_equal(extensions(rereadImage), extensions)
