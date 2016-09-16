#ifndef _RNIFTI_H_
#define _RNIFTI_H_

#include "niftilib/nifti1_io.h"

#ifdef __cplusplus
#include "lib/NiftiImage.h"

extern "C" {
#endif

extern void niftilib_register_all ();

#ifdef __cplusplus
}
#endif

#endif
