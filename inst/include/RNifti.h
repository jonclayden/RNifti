#ifndef _RNIFTI_H_
#define _RNIFTI_H_

// RNiftyReg and divest have used HAVE_R, so accept this variant for compatibility
#if !defined(USING_R) && defined(HAVE_R)
#define USING_R
#endif

#include "niftilib/nifti1_io.h"

#ifdef __cplusplus
#include "lib/NiftiImage.h"

#define HAVE_RNIFTI_NAMESPACE

extern "C" {
#endif

extern void niftilib_register_all ();

#ifdef __cplusplus
} // extern "C"
#endif

#endif
