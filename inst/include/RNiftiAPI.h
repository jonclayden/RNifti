#ifndef _RNIFTI_API_H_
#define _RNIFTI_API_H_

#include <R_ext/Rdynload.h>

#include "niftilib/nifti1_io.h"

#ifdef __cplusplus
extern "C" {
#endif

nifti_1_header * nifti_make_new_header (const int arg_dims[], int arg_dtype)
{
    static nifti_1_header*(*fun)(const int*, int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (nifti_1_header*(*)(const int*, int)) R_GetCCallable("RNifti","nii_make_new_header");
    return fun(arg_dims, arg_dtype);
}

nifti_image * nifti_make_new_nim (const int dims[], int datatype, int data_fill)
{
    static nifti_image*(*fun)(const int*, int, int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (nifti_image*(*)(const int*, int, int)) R_GetCCallable("RNifti","nii_make_new_nim");
    return fun(dims, datatype, data_fill);
}

nifti_image * nifti_convert_nhdr2nim (struct nifti_1_header nhdr, const char *fname)
{
    static nifti_image*(*fun)(struct nifti_1_header, const char *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (nifti_image*(*)(struct nifti_1_header, const char *)) R_GetCCallable("RNifti","nii_convert_nhdr2nim");
    return fun(nhdr, fname);
}

struct nifti_1_header nifti_convert_nim2nhdr (const nifti_image *nim)
{
    static struct nifti_1_header(*fun)(const nifti_image *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (struct nifti_1_header(*)(const nifti_image *)) R_GetCCallable("RNifti","nii_convert_nim2nhdr");
    return fun(nim);
}

nifti_image * nifti_copy_nim_info (const nifti_image *src)
{
    static nifti_image*(*fun)(const nifti_image*) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (nifti_image*(*)(const nifti_image*)) R_GetCCallable("RNifti","nii_copy_nim_info");
    return fun(src);
}

int nifti_copy_extensions (nifti_image *nim_dest, const nifti_image *nim_src)
{
    static int(*fun)(nifti_image*, const nifti_image*) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (int(*)(nifti_image*, const nifti_image*)) R_GetCCallable("RNifti","nii_copy_extensions");
    return fun(nim_dest, nim_src);
}

void nifti_image_free (nifti_image *nim)
{
    static void(*fun)(nifti_image*) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (void(*)(nifti_image*)) R_GetCCallable("RNifti","nii_image_free");
    fun(nim);
}

void nifti_datatype_sizes (int datatype, int *nbyper, int *swapsize)
{
    static void(*fun)(int, int*, int*) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (void(*)(int, int*, int*)) R_GetCCallable("RNifti","nii_datatype_sizes");
    return fun(datatype, nbyper, swapsize);
}

char const * nifti_datatype_string (int dt)
{
    static char const *(*fun)(int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (char const *(*)(int)) R_GetCCallable("RNifti","nii_datatype_string");
    return fun(dt);
}

char const * nifti_units_string (int uu)
{
    static char const *(*fun)(int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (char const *(*)(int)) R_GetCCallable("RNifti","nii_units_string");
    return fun(uu);
}

size_t nifti_get_volsize (const nifti_image *nim)
{
    static size_t(*fun)(const nifti_image *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (size_t(*)(const nifti_image *)) R_GetCCallable("RNifti","nii_get_volsize");
    return fun(nim);
}

int nifti_update_dims_from_array (nifti_image *nim)
{
    static int(*fun)(nifti_image *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (int(*)(nifti_image *)) R_GetCCallable("RNifti","nii_update_dims_from_array");
    return fun(nim);
}

int nifti_set_filenames (nifti_image *nim, const char *prefix, int check, int set_byte_order)
{
    static int(*fun)(nifti_image*, const char*, int, int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (int(*)(nifti_image*, const char*, int, int)) R_GetCCallable("RNifti","nii_set_filenames");
    return fun(nim, prefix, check, set_byte_order);
}

nifti_image * nifti_image_read (const char *hname, int read_data)
{
    static nifti_image *(*fun)(const char*, int) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (nifti_image *(*)(const char*, int)) R_GetCCallable("RNifti","nii_image_read");
    return fun(hname, read_data);
}

void nifti_image_write (nifti_image *nim)
{
    static void(*fun)(nifti_image *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (void(*)(nifti_image *)) R_GetCCallable("RNifti","nii_image_write");
    fun(nim);
}

float nifti_mat33_determ (mat33 R)
{
    static float(*fun)(mat33) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (float(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_determ");
    return fun(R);
}

mat33 nifti_mat33_inverse (mat33 R)
{
    static mat33(*fun)(mat33) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (mat33(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_inverse");
    return fun(R);
}

mat33 nifti_mat33_mul (mat33 A, mat33 B)
{
    static mat33(*fun)(mat33, mat33) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (mat33(*)(mat33, mat33)) R_GetCCallable("RNifti","nii_mat33_mul");
    return fun(A, B);
}

mat33 nifti_mat33_polar (mat33 A)
{
    static mat33(*fun)(mat33) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (mat33(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_polar");
    return fun(A);
}

mat44 nifti_mat44_inverse (mat44 R)
{
    static mat44(*fun)(mat44) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (mat44(*)(mat44)) R_GetCCallable("RNifti","nii_mat44_inverse");
    return fun(R);
}

void nifti_mat44_to_quatern (mat44 R, float *qb, float *qc, float *qd, float *qx, float *qy, float *qz, float *dx, float *dy, float *dz, float *qfac)
{
    static void(*fun)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (void(*)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *)) R_GetCCallable("RNifti","nii_mat44_to_quatern");
    fun(R, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac);
}

mat44 nifti_quatern_to_mat44 (float qb, float qc, float qd, float qx, float qy, float qz, float dx, float dy, float dz, float qfac)
{
    static mat44(*fun)(float, float, float, float, float, float, float, float, float, float) = NULL;
#ifdef _OPENMP
#pragma omp critical
#endif
    if (fun == NULL)
        fun = (mat44(*)(float, float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti","nii_quatern_to_mat44");
    return fun(qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif
