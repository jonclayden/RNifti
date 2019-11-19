#ifndef _RNIFTI_H_
#define _RNIFTI_H_

// RNiftyReg and divest have used HAVE_R, so accept this variant for compatibility
#if !defined(USING_R) && defined(HAVE_R)
#define USING_R
#endif

// Defined since RNifti v0.10.0
// Equal to 100 * (major version) + (minor version)
// Does not change with patch level, since the API should not change
#define RNIFTI_VERSION 100

#include "niftilib/nifti1_io.h"

#ifdef __cplusplus
#include "RNifti/NiftiImage.h"

// Defined since RNifti v0.3.0
#define HAVE_RNIFTI_NAMESPACE

extern "C" {
#endif

extern void niftilib_register_all ();

#ifdef __cplusplus
} // extern "C"
#endif

#endif
