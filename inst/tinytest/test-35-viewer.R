expect_null(view(system.file("extdata", "example.nii.gz", package="RNifti"), interactive=FALSE))

image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
expect_message(view(lyr(image,max=NA), interactive=FALSE), "(0, 748)")
expect_message(view(image[,,34], interactive=FALSE))

# Test fallback to a literal colour if the scale name isn't recognised
expect_null(view(image, lyr(image, mask=image>500, scale="red"), interactive=FALSE))

rgbImage <- readNifti(system.file("extdata", "example_rgb.nii.gz", package="RNifti"))
expect_equal(rgbImage$datatype, 128L)
expect_null(view(rgbImage, interactive=FALSE))
expect_null(view(lyr(rgbImage, mask=image<400), interactive=FALSE))

# NB: This conversion is flawed, since it does not restore sign, but it's fine for a test
vecImage <- structure(channels(rgbImage, c("red","green","blue"))/255, dim=c(dim(rgbImage),1L,3L))
vecImage <- updateNifti(vecImage, rgbImage)
vecImage$intent_code <- 1007
expect_null(view(image, vecImage, interactive=FALSE))

volImage <- readNifti(system.file("extdata", "example_4d.nii.gz", package="RNifti"))
expect_null(view(lyr(volImage,max=3), point=c(48,48,20), infoPanel=timeSeriesPanel, interactive=FALSE))

# The "label" argument overrides the default deparsed-expression label
expect_equal(lyr(image, label="Custom label")$label, "Custom label")
expect_equal(lyr(image)$label, "image")

# The "dict" argument must be named, and translates values to text in the info panel, falling
# back to the raw value for anything not present in the dictionary
expect_error(lyr(image, dict=c("Background","Region")), "named")
point <- c(48,48,30)
value <- image[point[1], point[2], point[3]]
dict <- structure("Test region", names=value)
captured <- NULL
capturePanel <- function(point, data, labels) captured <<- list(point=point, data=data, labels=labels)
expect_null(view(lyr(image,dict=dict), point=point, infoPanel=capturePanel, interactive=FALSE))
expect_equal(captured$point, point)
expect_equal(captured$data[[1]], paste0("Test region (", value, ")"))

missingDict <- structure("Should not be shown", names=value+1)
expect_null(view(lyr(image,dict=missingDict), point=point, infoPanel=capturePanel, interactive=FALSE))
expect_equal(captured$data[[1]], as.character(value))

# A named "scale" is a discrete colour lookup table rather than a continuous scale, and disables
# windowing; unnamed colours (the default) are unaffected and remain unnamed
scale <- structure(c("black","red","green"), names=c(0,value,value+1))
layer <- lyr(image, scale=scale)
expect_identical(layer$colours, scale)
expect_null(layer$window)
expect_true(!is.null(names(layer$colours)))
expect_null(names(lyr(image)$colours))

# The "alpha" argument requires the shades package, and adds an alpha channel to the colours,
# whether the scale is continuous or discrete
if (requireNamespace("shades", quietly=TRUE))
{
    expect_equal(nchar(lyr(image, alpha=0.5)$colours), rep(9L,100))
    expect_equal(nchar(lyr(image)$colours), rep(7L,100))
    expect_equal(unname(nchar(lyr(image, scale=scale, alpha=0.5)$colours)), rep(9L,3))
    expect_identical(names(lyr(image, scale=scale, alpha=0.5)$colours), names(scale))
    expect_null(view(image, lyr(image,scale=scale,dict=dict,alpha=0.5), interactive=FALSE))
}
