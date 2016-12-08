#ifndef _RNIFTI_H_
#define _RNIFTI_H_

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
