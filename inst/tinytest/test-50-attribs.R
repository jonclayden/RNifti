imagePath <- system.file("extdata", "example.nii.gz", package="RNifti")
jsonPath <- system.file("extdata", "example.json", package="RNifti")

image <- readNifti(imagePath)
expect_length(imageAttributes(image), 0L)

image <- readNifti(imagePath, json="read")
attribs <- imageAttributes(image)

# Should be one slice timing value per slice
expect_length(attribs$SliceTiming, dim(image)[3])

jsonString <- paste(readLines(jsonPath), collapse="\n")
jsonList <- jsonlite::fromJSON(jsonString)

expect_equivalent(jsonList, attribs[names(jsonList)])
index <- names(attribs)
# Round-trip conversion
expect_equivalent(attribs[index], RNifti:::renameToBids(RNifti:::renameFromBids(attribs))[index])
# Idempotence
expect_equivalent(attribs[index], RNifti:::renameToBids(attribs)[index])

# Setting `imageAttributes<-`(x,NULL) should scrub all but basic RNifti metadata
image <- readNifti(imagePath, json="convert")
copy <- image
imageAttributes(copy) <- NULL
expect_equal(attr(image,"pixdim"), attr(copy,"pixdim"))
expect_equal(attr(image,"echoTime"), 81)
expect_null(attr(copy, "echoTime"))

# Idempotence in the other direction
attribs <- imageAttributes(image)
index <- names(attribs)
expect_equivalent(attribs[index], RNifti:::renameFromBids(attribs)[index])
