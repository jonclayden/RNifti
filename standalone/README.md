# RNifti in C++ projects

`RNifti` is structured as an R package, but its code can also be used in standalone C and C++ projects. The contents of this directory (including the symlinked source files) provides a minimal example of just such a standalone project, and can therefore be copied wholesale—following all symlinks—to create a skeleton project to work from.

`RNifti` builds on, and provides an object-oriented wrapper for, the [reference NIfTI library for C](https://github.com/NIFTI-Imaging/nifti_clib), which are in the `niftilib` subdirectory. Some minimal modifications have been made to avoid name clashes between the NIfTI-1 and NIfTI-2 specialised versions of the library. For use in pure-C projects, the `RNifti.h` header helps with selecting between the NIfTI-1 and NIfTI-2 interfaces (via the `RNIFTI_NIFTILIB_VERSION` compile-time constant), but otherwise adds little.

For C++ projects, however, `RNifti` provides a more natural, object-oriented interface to the library, and various associated niceties. Some of these are outlined below.

## Basic usage

Including the main `RNifti` header will in turn include the appropriate `niftilib` headers.

```c++
// This constant determines which version of the library your client code sees,
// and should be set (if needed) before including RNifti.h. The default is 1.
#define RNIFTI_NIFTILIB_VERSION 2
#include "RNifti.h"
```

This provides access to the `niftilib` C functions, which mostly handle C-style `nifti_image` `struct`s, but in C++ code will also provide access to the `NiftiImage` and `NiftiImageData` wrapper classes. Full Doxygen documentation for these classes is available at <https://doxygen.flakery.org/RNifti/>.

## Reading, writing and creating images

The core image type, `NiftiImage`, can be constructed from a `nifti_image` (which it wraps around and provides memory management for), from another `NiftiImage`, or from a file path:

```c++
RNifti::NiftiImage image("example1.nii");
RNifti::NiftiImage header("example2.nii.gz", false);

```

Note that directly reading from gzipped NIfTI files is supported. A Boolean second argument specifies whether or not to read the pixel data from the file; if `false`, only metadata in the image header will be read. Alternatively, for images of more than three dimensions, a vector of integers can be specified to specify a subset of volumes to read.

Alternatively, an image can be created from its dimensions and a datatype. Data values can then be inserted.

```c++
const std::vector<int> dim = { 128, 128, 3 };
const size_t size = 128 * 128 * 3;

RNifti::NiftiImage image(dim, DT_UINT8);
std::vector<uint8_t> data(size, 1);
image.replaceData(data);
```

Writing back to file is done using the `toFile` member function. An output datatype can be specified if required.

```c++
image.toFile("output.nii.gz")
image.toFile("output_i16.nii.gz", DT_INT16);
```

## More on pixel data

The NIfTI file format stores image pixel values in a weakly typed blob whose interpretation depends on the image's declared datatype and intercept and slope values that are added to and multiplied by the stored elements to convert them to their final intended values. The `NiftiImageData` class handles and abstracts away much of this detail, and provides read and write access to the values, which are scaled appropriately if an intercept and slope applies. The class provides indexing and iterator support.

```c++
RNifti::NiftiImageData idata = image.data();
int sum = std::accumulate(idata.begin(), idata.end(), 0);
```

If the unscaled data values of a scaled image are needed, a variant data object can be obtained:

```c++
RNifti::NiftiImageData unscaled = idata.unscaled();
```

The data in an image can also be set or replaced through this route.

```c++
std::vector<uint8_t> data(size, 1);
RNifti::NiftiImageData idata(data.begin(), data.end(), DT_UINT8);
image.data() = idata;       // Or image.replaceData(idata);
```

## Interoperability with niftilib

An object of class `NiftiImage` can be implicitly converted to a `nifti_image*` (and actually contains one internally), and therefore can be passed to any `niftilib` standard function that expects this type. The `RNifti` and `niftilib` APIs can be used in parallel as long as a valid and consistent `nifti_image` is always produced.
