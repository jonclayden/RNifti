# RNifti: Fast R and C++ Access to NIfTI Images

The [NIfTI-1 format](http://www.nitrc.org/docman/view.php/26/64/nifti1.h) is a popular file format for storing medical imaging data, widely used in medical research and related fields. Conceptually, a NIfTI-1 file incorporates multidimensional numeric data, like an R `array`, but with additional metadata describing the real-space resolution of the image, the physical orientation of the image, and how the image should be interpreted.

There are several packages available for reading and writing NIfTI-1 files in R, and these are summarised in the [Medical Imaging task view](https://cran.r-project.org/view=MedicalImaging). However, `RNifti` is distinguished by its

- [extremely strong performance](#performance), in terms of speed;
- [C/C++ API](#api), allowing access to NIfTI images in compiled code in other R packages, or even in [standalone C++ code](#use-in-pure-c-projects); and
- modest dependencies, consisting of only R itself and the very widely-used [Rcpp](https://cran.r-project.org/package=Rcpp) C++ wrapper library.

The latest development version of the package can always be installed from GitHub using the `devtools` package.

```r
## install.packages("devtools")
devtools::install_github("jonclayden/RNifti")
```

## Usage

The primary role of `RNifti` is in reading and writing NIfTI-1 files, either `gzip`-compressed or uncompressed, and providing access to image data and metadata. An image may be read into R using the `readNifti` function.

```r
library(RNifti)
image <- readNifti(system.file("extdata", "example.nii.gz", package="RNifti"))
```

This image is an R array with some additional attributes containing information such as its dimensions and the size of its pixels (or voxels, in this case, since it is a 3D image). There are auxiliary functions for extracting this information: the standard `dim()`, plus `pixdim()` and `pixunits()`.

```r
dim(image)
# [1] 96 96 60

pixdim(image)
# [1] 2.5 2.5 2.5

pixunits(image)
# [1] "mm" "s"
```

So this image is of size 96 x 96 x 60 voxels, with each voxel representing 2.5 x 2.5 x 2.5 mm in real space. (The temporal unit, seconds here, only applies to the fourth dimension, if it is present.) Replacement versions of the latter functions are also available, for modifying the metadata.

A fuller list of the raw metadata stored in the file can be obtained using the `dumpNifti` function.

```r
dumpNifti(image)
# NIfTI-1 header
#     sizeof_hdr: 348
#       dim_info: 0
#            dim: 3  96  96  60  1  1  1  1
#      intent_p1: 0
#      intent_p2: 0
#      intent_p3: 0
#    intent_code: 0
#       datatype: 64
#         bitpix: 64
#    slice_start: 0
#         pixdim: -1.0  2.5  2.5  2.5  0.0  0.0  0.0  0.0
#     vox_offset: 352
#      scl_slope: 1
#      scl_inter: 0
#      slice_end: 0
#     slice_code: 0
#     xyzt_units: 10
#        cal_max: 2503
#        cal_min: 0
# slice_duration: 0
#        toffset: 0
#        descrip: FSL5.0
#       aux_file: 
#     qform_code: 2
#     sform_code: 2
#      quatern_b: 0
#      quatern_c: 1
#      quatern_d: 0
#      qoffset_x: 122.0339
#      qoffset_y: -95.18523
#      qoffset_z: -55.03814
#         srow_x: -2.5000  0.0000  0.0000  122.0339
#         srow_y: 0.00000  2.50000  0.00000  -95.18523
#         srow_z: 0.00000  0.00000  2.50000  -55.03814
#    intent_name: 
#          magic: n+1
```

Advanced users who know the NIfTI format well may want to alter elements of this metadata directly, and the `updateNifti` function provides a mechanism for this, either by passing a second, "template" image, or by providing lists with elements named as above, as in

```r
image <- updateNifti(image, list(intent_code=1L))
```

An image can be written back to NIfTI-1 format using the `writeNifti` function.

```r
writeNifti(image, "file.nii.gz")
```

## Performance

The `RNifti` package uses the robust NIfTI-1 reference implementation, which is written in C, to read and write NIfTI files. It also uses the standard NIfTI-1 data structure as its canonical representation of a file in memory. Together, these make the package extremely fast, as the following benchmark against packages [`AnalyzeFMRI`](https://cran.r-project.org/package=AnalyzeFMRI), [`ANTsR`](https://github.com/stnava/ANTsR),  [`neuroim`](https://cran.r-project.org/package=neuroim), [`oro.nifti`](https://cran.r-project.org/package=oro.nifti) and [`tractor.base`](https://cran.r-project.org/package=tractor.base) shows.

```r
installed.packages()[c("AnalyzeFMRI","ANTsR","neuroim","oro.nifti","RNifti",
                       "tractor.base"), "Version"]
#  AnalyzeFMRI        ANTsR      neuroim    oro.nifti       RNifti tractor.base
#     "1.1-16"      "0.3.3"      "0.0.6"    "0.5.5.2"      "0.1.0"      "3.0.5"

library(microbenchmark)
microbenchmark(AnalyzeFMRI::f.read.volume("example.nii"),
               ANTsR::antsImageRead("example.nii"),
               neuroim::loadVolume("example.nii"),
               oro.nifti::readNIfTI("example.nii"),
               RNifti::readNifti("example.nii"),
               tractor.base::readImageFile("example.nii"))
# Output level is not set; defaulting to "Info"
# Unit: milliseconds
#                                        expr       min        lq       mean
#   AnalyzeFMRI::f.read.volume("example.nii") 25.532280 27.216641  46.646448
#         ANTsR::antsImageRead("example.nii")  1.805555  2.282496   3.430417
#          neuroim::loadVolume("example.nii") 33.639202 36.367228  87.101533
#         oro.nifti::readNIfTI("example.nii") 45.836130 49.692145 126.697245
#            RNifti::readNifti("example.nii")  1.312822  1.646371   7.239342
#  tractor.base::readImageFile("example.nii") 38.091666 39.628323  48.573645
#      median         uq        max neval
#   28.486199  31.441231  535.60771   100
#    2.434518   2.614578   96.20713   100
#   38.556286  56.650639 2071.82043   100
#  154.908146 162.407670  494.26035   100
#    1.970080   3.104831  122.78295   100
#   40.752468  42.635031  249.59664   100
```

With a median runtime of less than 2 ms, `RNifti` is typically at least ten times as fast as the alternatives to read this image into R. The exception is `ANTsR`, which uses a similar low-level pointer-based arrangement as `RNifti`, and is therefore comparable in speed. However, `ANTsR` has substantial dependencies, which may affect its suitability in some applications.

## Implementation details

The package does not fully duplicate the NIfTI-1 structure's contents in R-visible objects. Instead, it passes key metadata back to R, such as the image dimensions and pixel dimensions, and it also passes back the pixel values where they are needed. It also creates an [external pointer](http://r-manuals.flakery.org/R-exts.html#External-pointers-and-weak-references) to the native data structure, which is stored in an attribute. This pointer is dereferenced whenever the object is passed back to the C++ code, thereby avoiding unnecessary duplication and ensuring that all metadata remains intact. The full NIfTI-1 header can be obtained using the `dumpNifti` R function, if it is needed.

This arrangement is efficient and generally works well, but many R operations strip attributes—in which case the external pointer will be removed. The internal structure will be built again when necessary, but using default metadata. In these cases, if it is important to keep the original metadata, the `updateNifti` function should be called explicitly, with a template object. This reconstructs the NIfTI-1 data structure, using the template as a starting point.

## API

It is possible to use the package's NIfTI-handling code in other R packages' compiled code, thereby obviating the need to duplicate the reference implementation. Moreover, `RNifti` provides a C++ wrapper class, `NiftiImage`, which simplifies memory management, supports the package's internal image pointers and persistence, and provides syntactic sugar. Full doxygen documentation for this class is available at <http://doxygen.flakery.org/RNifti/>, and is also provided with package releases.

A third-party package can use the `NiftiImage` class by including

```
LinkingTo: Rcpp, RNifti
```

in its `DESCRIPTION` file, and then including the `RNifti.h` header file. For example,

```c++
#include "RNifti.h"

void myfunction ()
{
    RNifti::NiftiImage image("example.nii.gz");
    // Do something with the image
}
```

If you're using the `sourceCpp` function from `Rcpp`, you may also need to add the attribute line

```c++
// [[Rcpp::depends(RNifti)]]
```

to the top of your C++ source file.

In addition to the one taking a file path, there are also constructors taking a `SEXP` (i.e., an R object), another `NiftiImage`, or a `nifti_image` structure from the reference implementation. `NiftiImage` objects can be implicitly cast to pointers to `nifti_image` structs, meaning that they can be directly used in calls to the reference implementation's own API. The latter is accessed through the separate `RNiftiAPI.h` header file.

```c++
#include "RNifti.h"
#include "RNiftiAPI.h"

void myfunction (SEXP image_)
{
    RNifti::NiftiImage image(image_);
    const size_t volsize = nifti_get_volsize(image);
}
```

(`RNifti` will also have to be added to the `Imports` list in the package's `DESCRIPTION` file, as well as `LinkingTo`.) The `RNiftiAPI.h` header should only be included once per package, since it contains function implementations. Multiple includes will lead to duplicate symbol warnings from your linker. Therefore, if multiple source files require access to the NIfTI-1 reference implementation, it is recommended that the API header be included alone in a separate ".c" or ".cpp" file, while others only include the main `RNifti.h`.

`RNifti` is not specifically designed to be thread-safe, and R itself is expressly single-threaded. However, some effort has been made to try to minimise problems associated with parallelisation, such as putting R API calls within a critical region if OpenMP is being used. If you are using the API in a package that does use OpenMP or another form of threads, it is wise to preregister the functions exported by `RNifti` before use, by calling `niftilib_register_all()`. In single-threaded contexts this is optional, and will be performed when required.

## Use in pure C++ projects

Thanks to contributions from [@soolijoo](https://github.com/soolijoo), it is possible (as of package version 0.7.0) to use the `NiftiImage` C++ class in standalone C++ projects. You will need the following files:

| Path                            | Purpose                                                                                             |
| ------------------------------- | --------------------------------------------------------------------------------------------------- |
| `inst/include/lib/*.h`          | Headers defining the `NiftiImage` class itself, related functions and macros                        |
| `inst/include/niftilib/*.h`     | Headers for the NIfTI-1 reference implementation                                                    |
| `inst/include/znzlib/znzlib.h`  | Header for I/O functions from the NIfTI-1 reference implementation                                  |
| `inst/include/zlib/*.h`         | `zlib` headers for reading and writing gzipped files (optional; system `zlib` can be used instead)  |
| `src/niftilib/nifti1_io.c`      | Source file for the NIfTI-1 reference implementation                                                |
| `src/znzlib/znzlib.c`           | Source for I/O functions from the NIfTI-1 reference implementation                                  |
| `src/zlib/*`                    | `zlib` source files for reading and writing gzipped files (optional, as above)                      |

Note that the `NiftiImage` class is header-only, but C code from the NIfTI-1 reference implementation will need to be compiled and linked into the project. The `print.h` header should be included before including `NiftiImage.h`, so that the R API is not used for printing error messages. The [`standalone` directory](https://github.com/jonclayden/RNifti/tree/master/standalone) provides a minimal example.

## The NIfTI-2 format

The [NIfTI-2 format](https://nifti.nimh.nih.gov/nifti-2) is a evolution of the far more widely-used NIfTI-1. It primarily [uses wider types for various fields](https://www.nitrc.org/forum/forum.php?thread_id=2148&forum_id=1941), to support large datasets and improve precision.

Unfortunately, the NIfTI-1 version of the reference library is not forwards-compatible with NIfTI-2, and does not recognise NIfTI-2 files as valid, while the NIfTI-2 version changes the definition of the `nifti_image` data structure, and hence the return type of several core functions, rendering it potentially incompatible with software written for the original library. As a result, adding full NIfTI-2 support to `RNifti` without breaking existing code is not straightforward. Nevertheless, as of `RNifti` version 0.8.0, R function `niftiVersion()` and C++ static method `NiftiImage::fileVersion()` offer a forwards-compatible way to determine the version of the format used by a particular file, so that calling functions can take action accordingly.
