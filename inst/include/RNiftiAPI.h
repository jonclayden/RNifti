#ifndef _RNIFTI_API_H_
#define _RNIFTI_API_H_

#include <R_ext/Rdynload.h>

#include "niftilib/nifti1_io.h"

#ifdef __cplusplus
extern "C" {
#endif

static nifti_1_header*(*_nifti_make_new_header)(const int*, int) = NULL;
static nifti_image*(*_nifti_make_new_nim)(const int*, int, int) = NULL;
static nifti_image*(*_nifti_convert_nhdr2nim)(struct nifti_1_header, const char *) = NULL;
static struct nifti_1_header(*_nifti_convert_nim2nhdr)(const nifti_image *) = NULL;
static nifti_image*(*_nifti_copy_nim_info)(const nifti_image*) = NULL;
static int(*_nifti_copy_extensions)(nifti_image*, const nifti_image*) = NULL;
static void(*_nifti_image_unload)(nifti_image*) = NULL;
static void(*_nifti_image_free)(nifti_image*) = NULL;

static void(*_nifti_datatype_sizes)(int, int*, int*) = NULL;
static char const *(*_nifti_datatype_string)(int) = NULL;
static char const *(*_nifti_units_string)(int) = NULL;
static size_t(*_nifti_get_volsize)(const nifti_image *) = NULL;
static int(*_nifti_update_dims_from_array)(nifti_image *) = NULL;

static int(*_nifti_set_filenames)(nifti_image*, const char*, int, int) = NULL;
static nifti_image *(*_nifti_image_read)(const char*, int) = NULL;
static void(*_nifti_image_write)(nifti_image *) = NULL;

static float(*_nifti_mat33_rownorm)(mat33) = NULL;
static float(*_nifti_mat33_colnorm)(mat33) = NULL;
static float(*_nifti_mat33_determ)(mat33) = NULL;
static mat33(*_nifti_mat33_inverse)(mat33) = NULL;
static mat33(*_nifti_mat33_mul)(mat33, mat33) = NULL;
static mat33(*_nifti_mat33_polar)(mat33) = NULL;
static mat44(*_nifti_mat44_inverse)(mat44) = NULL;
static void(*_nifti_mat44_to_quatern)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *) = NULL;
static mat44(*_nifti_quatern_to_mat44)(float, float, float, float, float, float, float, float, float, float) = NULL;
static void(*_nifti_mat44_to_orientation)(mat44, int *, int *, int *) = NULL;

static void(*_nifti_swap_2bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_4bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_8bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_16bytes)(size_t, void *) = NULL;

void niftilib_register_all ()
{
#ifdef _OPENMP
#pragma omp critical
#endif
    if (_nifti_make_new_header == NULL)
    {
        _nifti_make_new_header = (nifti_1_header*(*)(const int*, int)) R_GetCCallable("RNifti","nii_make_new_header");
        _nifti_make_new_nim = (nifti_image*(*)(const int*, int, int)) R_GetCCallable("RNifti","nii_make_new_nim");
        _nifti_convert_nhdr2nim = (nifti_image*(*)(struct nifti_1_header, const char *)) R_GetCCallable("RNifti","nii_convert_nhdr2nim");
        _nifti_convert_nim2nhdr = (struct nifti_1_header(*)(const nifti_image *)) R_GetCCallable("RNifti","nii_convert_nim2nhdr");
        _nifti_copy_nim_info = (nifti_image*(*)(const nifti_image*)) R_GetCCallable("RNifti","nii_copy_nim_info");
        _nifti_copy_extensions = (int(*)(nifti_image*, const nifti_image*)) R_GetCCallable("RNifti","nii_copy_extensions");
        _nifti_image_unload = (void(*)(nifti_image*)) R_GetCCallable("RNifti","nii_image_unload");
        _nifti_image_free = (void(*)(nifti_image*)) R_GetCCallable("RNifti","nii_image_free");
        
        _nifti_datatype_sizes = (void(*)(int, int*, int*)) R_GetCCallable("RNifti","nii_datatype_sizes");
        _nifti_datatype_string = (char const *(*)(int)) R_GetCCallable("RNifti","nii_datatype_string");
        _nifti_units_string = (char const *(*)(int)) R_GetCCallable("RNifti","nii_units_string");
        _nifti_get_volsize = (size_t(*)(const nifti_image *)) R_GetCCallable("RNifti","nii_get_volsize");
        _nifti_update_dims_from_array = (int(*)(nifti_image *)) R_GetCCallable("RNifti","nii_update_dims_from_array");
        
        _nifti_set_filenames = (int(*)(nifti_image*, const char*, int, int)) R_GetCCallable("RNifti","nii_set_filenames");
        _nifti_image_read = (nifti_image *(*)(const char*, int)) R_GetCCallable("RNifti","nii_image_read");
        _nifti_image_write = (void(*)(nifti_image *)) R_GetCCallable("RNifti","nii_image_write");
        
        _nifti_mat33_rownorm = (float(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_rownorm");
        _nifti_mat33_colnorm = (float(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_colnorm");
        _nifti_mat33_determ = (float(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_determ");
        _nifti_mat33_inverse = (mat33(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_inverse");
        _nifti_mat33_mul = (mat33(*)(mat33, mat33)) R_GetCCallable("RNifti","nii_mat33_mul");
        _nifti_mat33_polar = (mat33(*)(mat33)) R_GetCCallable("RNifti","nii_mat33_polar");
        _nifti_mat44_inverse = (mat44(*)(mat44)) R_GetCCallable("RNifti","nii_mat44_inverse");
        _nifti_mat44_to_quatern = (void(*)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *)) R_GetCCallable("RNifti","nii_mat44_to_quatern");
        _nifti_quatern_to_mat44 = (mat44(*)(float, float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti","nii_quatern_to_mat44");
        _nifti_mat44_to_orientation = (void(*)(mat44, int *, int *, int *)) R_GetCCallable("RNifti","nii_mat44_to_orientation");
        
        _nifti_swap_2bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti","nii_swap_2bytes");
        _nifti_swap_4bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti","nii_swap_4bytes");
        _nifti_swap_8bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti","nii_swap_8bytes");
        _nifti_swap_16bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti","nii_swap_16bytes");
    }
}

nifti_1_header * nifti_make_new_header (const int arg_dims[], int arg_dtype)
{
    if (_nifti_make_new_header == NULL)
        niftilib_register_all();
    return _nifti_make_new_header(arg_dims, arg_dtype);
}

nifti_image * nifti_make_new_nim (const int dims[], int datatype, int data_fill)
{
    if (_nifti_make_new_nim == NULL)
        niftilib_register_all();
    return _nifti_make_new_nim(dims, datatype, data_fill);
}

nifti_image * nifti_convert_nhdr2nim (struct nifti_1_header nhdr, const char *fname)
{
    if (_nifti_convert_nhdr2nim == NULL)
        niftilib_register_all();
    return _nifti_convert_nhdr2nim(nhdr, fname);
}

struct nifti_1_header nifti_convert_nim2nhdr (const nifti_image *nim)
{
    if (_nifti_convert_nim2nhdr == NULL)
        niftilib_register_all();
    return _nifti_convert_nim2nhdr(nim);
}

nifti_image * nifti_copy_nim_info (const nifti_image *src)
{
    if (_nifti_copy_nim_info == NULL)
        niftilib_register_all();
    return _nifti_copy_nim_info(src);
}

int nifti_copy_extensions (nifti_image *nim_dest, const nifti_image *nim_src)
{
    if (_nifti_copy_extensions == NULL)
        niftilib_register_all();
    return _nifti_copy_extensions(nim_dest, nim_src);
}

void nifti_image_unload (nifti_image *nim)
{
    if (_nifti_image_unload == NULL)
        niftilib_register_all();
    _nifti_image_unload(nim);
}

void nifti_image_free (nifti_image *nim)
{
    if (_nifti_image_free == NULL)
        niftilib_register_all();
    _nifti_image_free(nim);
}

void nifti_datatype_sizes (int datatype, int *nbyper, int *swapsize)
{
    if (_nifti_datatype_sizes == NULL)
        niftilib_register_all();
    return _nifti_datatype_sizes(datatype, nbyper, swapsize);
}

char const * nifti_datatype_string (int dt)
{
    if (_nifti_datatype_string == NULL)
        niftilib_register_all();
    return _nifti_datatype_string(dt);
}

char const * nifti_units_string (int uu)
{
    if (_nifti_units_string == NULL)
        niftilib_register_all();
    return _nifti_units_string(uu);
}

size_t nifti_get_volsize (const nifti_image *nim)
{
    if (_nifti_get_volsize == NULL)
        niftilib_register_all();
    return _nifti_get_volsize(nim);
}

int nifti_update_dims_from_array (nifti_image *nim)
{
    if (_nifti_update_dims_from_array == NULL)
        niftilib_register_all();
    return _nifti_update_dims_from_array(nim);
}

int nifti_set_filenames (nifti_image *nim, const char *prefix, int check, int set_byte_order)
{
    if (_nifti_set_filenames == NULL)
        niftilib_register_all();
    return _nifti_set_filenames(nim, prefix, check, set_byte_order);
}

nifti_image * nifti_image_read (const char *hname, int read_data)
{
    if (_nifti_image_read == NULL)
        niftilib_register_all();
    return _nifti_image_read(hname, read_data);
}

void nifti_image_write (nifti_image *nim)
{
    if (_nifti_image_write == NULL)
        niftilib_register_all();
    _nifti_image_write(nim);
}

float nifti_mat33_rownorm (mat33 A)
{
    if (_nifti_mat33_rownorm == NULL)
        niftilib_register_all();
    return _nifti_mat33_rownorm(A);
}

float nifti_mat33_colnorm (mat33 A)
{
    if (_nifti_mat33_colnorm == NULL)
        niftilib_register_all();
    return _nifti_mat33_colnorm(A);
}

float nifti_mat33_determ (mat33 R)
{
    if (_nifti_mat33_determ == NULL)
        niftilib_register_all();
    return _nifti_mat33_determ(R);
}

mat33 nifti_mat33_inverse (mat33 R)
{
    if (_nifti_mat33_inverse == NULL)
        niftilib_register_all();
    return _nifti_mat33_inverse(R);
}

mat33 nifti_mat33_mul (mat33 A, mat33 B)
{
    if (_nifti_mat33_mul == NULL)
        niftilib_register_all();
    return _nifti_mat33_mul(A, B);
}

mat33 nifti_mat33_polar (mat33 A)
{
    if (_nifti_mat33_polar == NULL)
        niftilib_register_all();
    return _nifti_mat33_polar(A);
}

mat44 nifti_mat44_inverse (mat44 R)
{
    if (_nifti_mat44_inverse == NULL)
        niftilib_register_all();
    return _nifti_mat44_inverse(R);
}

void nifti_mat44_to_quatern (mat44 R, float *qb, float *qc, float *qd, float *qx, float *qy, float *qz, float *dx, float *dy, float *dz, float *qfac)
{
    if (_nifti_mat44_to_quatern == NULL)
        niftilib_register_all();
    _nifti_mat44_to_quatern(R, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac);
}

mat44 nifti_quatern_to_mat44 (float qb, float qc, float qd, float qx, float qy, float qz, float dx, float dy, float dz, float qfac)
{
    if (_nifti_quatern_to_mat44 == NULL)
        niftilib_register_all();
    return _nifti_quatern_to_mat44(qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac);
}

void nifti_mat44_to_orientation(mat44 R, int *icod, int *jcod, int *kcod)
{
    if (_nifti_mat44_to_orientation == NULL)
        niftilib_register_all();
    return _nifti_mat44_to_orientation(R, icod, jcod, kcod);
}

void nifti_swap_2bytes (size_t n, void *ar)
{
    if (_nifti_swap_2bytes == NULL)
        niftilib_register_all();
    _nifti_swap_2bytes(n, ar);
}

void nifti_swap_4bytes (size_t n, void *ar)
{
    if (_nifti_swap_4bytes == NULL)
        niftilib_register_all();
    _nifti_swap_4bytes(n, ar);
}

void nifti_swap_8bytes (size_t n, void *ar)
{
    if (_nifti_swap_8bytes == NULL)
        niftilib_register_all();
    _nifti_swap_8bytes(n, ar);
}

void nifti_swap_16bytes (size_t n, void *ar)
{
    if (_nifti_swap_16bytes == NULL)
        niftilib_register_all();
    _nifti_swap_16bytes(n, ar);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif
