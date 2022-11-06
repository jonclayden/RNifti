#ifndef _RNIFTI_API_H_
#define _RNIFTI_API_H_

#include <R_ext/Rdynload.h>
#include "RNifti.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NIFTILIB_WRAPPER_BODY(symbol, ...)  \
    if (symbol == NULL)                     \
        niftilib_register_all();            \
    return symbol(__VA_ARGS__);

#define NIFTILIB_WRAPPER_BODY_VOID(symbol)  \
    if (symbol == NULL)                     \
        niftilib_register_all();            \
    return symbol();

#define NIFTILIB_WRAPPER_BODY_NORETURN(symbol, ...) \
    if (symbol == NULL)                             \
        niftilib_register_all();                    \
    symbol(__VA_ARGS__);

#define NIFTILIB_WRAPPER_BODY_VOID_NORETURN(symbol) \
    if (symbol == NULL)                             \
        niftilib_register_all();                    \
    symbol();

static int registered = 0;

#if RNIFTI_NIFTILIB_VERSION == 1
static char const *(*_nifti_datatype_string)(int) = NULL;
static char const *(*_nifti_units_string)(int) = NULL;
static char const *(*_nifti_intent_string)(int) = NULL;
static char const *(*_nifti_xform_string)(int) = NULL;
static char const *(*_nifti_slice_string)(int) = NULL;
static char const *(*_nifti_orientation_string)(int) = NULL;
static int(*_nifti_is_inttype)(int) = NULL;
static mat44(*_nifti_mat44_inverse)(mat44) = NULL;
static mat33(*_nifti_mat33_inverse)(mat33) = NULL;
static mat33(*_nifti_mat33_polar)(mat33) = NULL;
static float(*_nifti_mat33_rownorm)(mat33) = NULL;
static float(*_nifti_mat33_colnorm)(mat33) = NULL;
static float(*_nifti_mat33_determ)(mat33) = NULL;
static mat33(*_nifti_mat33_mul)(mat33, mat33) = NULL;
static void(*_nifti_swap_2bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_4bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_8bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_16bytes)(size_t, void *) = NULL;
static void(*_nifti_swap_Nbytes)(size_t, int, void *) = NULL;
static int(*_nifti_datatype_is_valid)(int, int) = NULL;
static int(*_nifti_datatype_from_string)(const char *) = NULL;
static const char *(*_nifti_datatype_to_string)(int) = NULL;
static int(*_nifti_get_filesize)(const char *) = NULL;
static void(*_swap_nifti_header)(struct nifti_1_header *, int) = NULL;
static void(*_old_swap_nifti_header)(struct nifti_1_header *, int) = NULL;
static int(*_nifti_swap_as_analyze)(nifti_analyze75 *) = NULL;
static nifti_image *(*_nifti_image_read_bricks)(const char *, int, const int *, nifti_brick_list *) = NULL;
static int(*_nifti_image_load_bricks)(nifti_image *, int, const int *, nifti_brick_list *) = NULL;
static void(*_nifti_free_NBL)(nifti_brick_list *) = NULL;
static nifti_image *(*_nifti_image_read)(const char *, int) = NULL;
static int(*_nifti_image_load)(nifti_image *) = NULL;
static void(*_nifti_image_unload)(nifti_image *) = NULL;
static void(*_nifti_image_free)(nifti_image *) = NULL;
static int(*_nifti_read_collapsed_image)(nifti_image *, const int[8], void **) = NULL;
static int(*_nifti_read_subregion_image)(nifti_image *, const int *, const int *, void **) = NULL;
static void(*_nifti_image_write)(nifti_image *) = NULL;
static void(*_nifti_image_write_bricks)(nifti_image *, const nifti_brick_list *) = NULL;
static void(*_nifti_image_infodump)(const nifti_image *) = NULL;
static void(*_nifti_disp_lib_hist)(void) = NULL;
static void(*_nifti_disp_lib_version)(void) = NULL;
static int(*_nifti_disp_matrix_orient)(const char *, mat44) = NULL;
static int(*_nifti_disp_type_list)(int) = NULL;
static char *(*_nifti_image_to_ascii)(const nifti_image *) = NULL;
static nifti_image *(*_nifti_image_from_ascii)(const char *, int *) = NULL;
static size_t(*_nifti_get_volsize)(const nifti_image *) = NULL;
static int(*_nifti_set_filenames)(nifti_image *, const char *, int, int) = NULL;
static char *(*_nifti_makehdrname)(const char *, int, int, int) = NULL;
static char *(*_nifti_makeimgname)(const char *, int, int, int) = NULL;
static int(*_is_nifti_file)(const char *) = NULL;
static char *(*_nifti_find_file_extension)(const char *) = NULL;
static int(*_nifti_is_complete_filename)(const char*) = NULL;
static int(*_nifti_validfilename)(const char*) = NULL;
static int(*_disp_nifti_1_header)(const char *, const nifti_1_header *) = NULL;
static void(*_nifti_set_debug_level)(int) = NULL;
static void(*_nifti_set_skip_blank_ext)(int) = NULL;
static void(*_nifti_set_allow_upper_fext)(int) = NULL;
static int(*_valid_nifti_brick_list)(nifti_image *, int, const int *, int) = NULL;
static znzFile(*_nifti_image_open)(const char *, const char *, nifti_image **) = NULL;
static znzFile(*_nifti_image_write_hdr_img)(nifti_image *, int, const char*) = NULL;
static znzFile(*_nifti_image_write_hdr_img2)(nifti_image *, int, const char*, znzFile, const nifti_brick_list *) = NULL;
static size_t(*_nifti_read_buffer)(znzFile, void*, size_t, nifti_image *) = NULL;
static int(*_nifti_write_all_data)(znzFile, nifti_image *, const nifti_brick_list *) = NULL;
static size_t(*_nifti_write_buffer)(znzFile, const void *, size_t) = NULL;
static nifti_image *(*_nifti_read_ascii_image)(znzFile, char *, int, int) = NULL;
static znzFile(*_nifti_write_ascii_image)(nifti_image *, const nifti_brick_list *, const char *, int, int) = NULL;
static void(*_nifti_datatype_sizes)(int, int *, int *) = NULL;
static void(*_nifti_mat44_to_quatern)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *) = NULL;
static mat44(*_nifti_quatern_to_mat44)(float, float, float, float, float, float, float, float, float, float) = NULL;
static mat44(*_nifti_make_orthog_mat44)(float, float, float, float, float, float, float, float, float) = NULL;
static int(*_nifti_short_order)(void) = NULL;
static void(*_nifti_mat44_to_orientation)(mat44, int *, int *, int *) = NULL;
static char *(*_nifti_findhdrname)(const char*) = NULL;
static char *(*_nifti_findimgname)(const char*, int) = NULL;
static int(*_nifti_is_gzfile)(const char*) = NULL;
static char *(*_nifti_makebasename)(const char*) = NULL;
static struct nifti_1_header(*_nifti_convert_nim2nhdr)(const nifti_image*) = NULL;
static nifti_1_header *(*_nifti_make_new_header)(const int[], int) = NULL;
static nifti_1_header *(*_nifti_read_header)(const char *, int *, int) = NULL;
static nifti_image *(*_nifti_copy_nim_info)(const nifti_image *) = NULL;
static nifti_image *(*_nifti_make_new_nim)(const int[], int, int) = NULL;
static nifti_image *(*_nifti_simple_init_nim)(void) = NULL;
static nifti_image *(*_nifti_convert_nhdr2nim)(struct nifti_1_header, const char *) = NULL;
static int(*_nifti_hdr_looks_good)(const nifti_1_header *) = NULL;
static int(*_nifti_is_valid_datatype)(int) = NULL;
static int(*_nifti_is_valid_ecode)(int) = NULL;
static int(*_nifti_nim_is_valid)(nifti_image *, int) = NULL;
static int(*_nifti_nim_has_valid_dims)(nifti_image *, int) = NULL;
static int(*_is_valid_nifti_type)(int) = NULL;
static int(*_nifti_test_datatype_sizes)(int) = NULL;
static int(*_nifti_type_and_names_match)(nifti_image *, int) = NULL;
static int(*_nifti_update_dims_from_array)(nifti_image *) = NULL;
static void(*_nifti_set_iname_offset)(nifti_image *) = NULL;
static int(*_nifti_set_type_from_names)(nifti_image *) = NULL;
static int(*_nifti_add_extension)(nifti_image *, const char *, int, int) = NULL;
static int(*_nifti_compiled_with_zlib)(void) = NULL;
static int(*_nifti_copy_extensions)(nifti_image *, const nifti_image *) = NULL;
static int(*_nifti_free_extensions)(nifti_image *) = NULL;
static int *(*_nifti_get_intlist)(int, const char *) = NULL;
static char *(*_nifti_strdup)(const char *) = NULL;
static int(*_valid_nifti_extensions)(const nifti_image *) = NULL;
#elif RNIFTI_NIFTILIB_VERSION == 2
static char const *(*_nifti_datatype_string)(int) = NULL;
static char const *(*_nifti_units_string)(int) = NULL;
static char const *(*_nifti_intent_string)(int) = NULL;
static char const *(*_nifti_xform_string)(int) = NULL;
static char const *(*_nifti_slice_string)(int) = NULL;
static char const *(*_nifti_orientation_string)(int) = NULL;
static int(*_nifti_is_inttype)(int) = NULL;
static mat44(*_nifti_mat44_inverse)(mat44) = NULL;
static mat44(*_nifti_mat44_mul)(mat44, mat44) = NULL;
static nifti_dmat44(*_nifti_dmat44_inverse)(nifti_dmat44) = NULL;
static int(*_nifti_mat44_to_dmat44)(mat44 *, nifti_dmat44 *) = NULL;
static int(*_nifti_dmat44_to_mat44)(nifti_dmat44 *, mat44 *) = NULL;
static nifti_dmat44(*_nifti_dmat44_mul)(nifti_dmat44, nifti_dmat44) = NULL;
static nifti_dmat33(*_nifti_dmat33_inverse)(nifti_dmat33) = NULL;
static nifti_dmat33(*_nifti_dmat33_polar)(nifti_dmat33) = NULL;
static double(*_nifti_dmat33_rownorm)(nifti_dmat33) = NULL;
static double(*_nifti_dmat33_colnorm)(nifti_dmat33) = NULL;
static double(*_nifti_dmat33_determ)(nifti_dmat33) = NULL;
static nifti_dmat33(*_nifti_dmat33_mul)(nifti_dmat33, nifti_dmat33) = NULL;
static mat33(*_nifti_mat33_inverse)(mat33) = NULL;
static mat33(*_nifti_mat33_polar)(mat33) = NULL;
static float(*_nifti_mat33_rownorm)(mat33) = NULL;
static float(*_nifti_mat33_colnorm)(mat33) = NULL;
static float(*_nifti_mat33_determ)(mat33) = NULL;
static mat33(*_nifti_mat33_mul)(mat33, mat33) = NULL;
static void(*_nifti_swap_2bytes)(int64_t, void *) = NULL;
static void(*_nifti_swap_4bytes)(int64_t, void *) = NULL;
static void(*_nifti_swap_8bytes)(int64_t, void *) = NULL;
static void(*_nifti_swap_16bytes)(int64_t, void *) = NULL;
static void(*_nifti_swap_Nbytes)(int64_t, int, void *) = NULL;
static int(*_nifti_datatype_is_valid)(int, int) = NULL;
static int(*_nifti_datatype_from_string)(const char *) = NULL;
static const char *(*_nifti_datatype_to_string)(int) = NULL;
static int(*_nifti_header_version)(const char *, size_t) = NULL;
static int64_t(*_nifti2_get_filesize)(const char *) = NULL;
static void(*_swap_nifti_header)(void *, int) = NULL;
static void(*_old_swap_nifti_header)(struct nifti_1_header *, int) = NULL;
static void(*_nifti_swap_as_analyze)(nifti_analyze75 *) = NULL;
static void(*_nifti_swap_as_nifti1)(nifti_1_header *) = NULL;
static void(*_nifti_swap_as_nifti2)(nifti_2_header *) = NULL;
static nifti_image *(*_nifti2_image_read_bricks)(const char *, int64_t, const int64_t *, nifti_brick_list *) = NULL;
static int(*_nifti2_image_load_bricks)(nifti_image *, int64_t, const int64_t *, nifti_brick_list *) = NULL;
static void(*_nifti2_free_NBL)(nifti_brick_list *) = NULL;
static nifti_image *(*_nifti2_image_read)(const char *, int) = NULL;
static int(*_nifti2_image_load)(nifti_image *) = NULL;
static void(*_nifti2_image_unload)(nifti_image *) = NULL;
static void(*_nifti2_image_free)(nifti_image *) = NULL;
static int64_t(*_nifti2_read_collapsed_image)(nifti_image *, const int64_t[8], void **) = NULL;
static int64_t(*_nifti2_read_subregion_image)(nifti_image *, const int64_t *, const int64_t *, void **) = NULL;
static void(*_nifti2_image_write)(nifti_image *) = NULL;
static void(*_nifti2_image_write_bricks)(nifti_image *, const nifti_brick_list *) = NULL;
static void(*_nifti2_image_infodump)(const nifti_image *) = NULL;
static void(*_nifti2_disp_lib_hist)(int) = NULL;
static void(*_nifti_disp_lib_version)(void) = NULL;
static int(*_nifti2_disp_matrix_orient)(const char *, nifti_dmat44) = NULL;
static int(*_nifti_disp_type_list)(int) = NULL;
static char *(*_nifti2_image_to_ascii)(const nifti_image *) = NULL;
static nifti_image *(*_nifti2_image_from_ascii)(const char *, int *) = NULL;
static int64_t(*_nifti2_get_volsize)(const nifti_image *) = NULL;
static int(*_nifti2_set_filenames)(nifti_image *, const char *, int, int) = NULL;
static char *(*_nifti_makehdrname)(const char *, int, int, int) = NULL;
static char *(*_nifti_makeimgname)(const char *, int, int, int) = NULL;
static int(*_is_nifti_file)(const char *) = NULL;
static char *(*_nifti_find_file_extension)(const char *) = NULL;
static int(*_nifti_is_complete_filename)(const char*) = NULL;
static int(*_nifti_validfilename)(const char*) = NULL;
static int(*_disp_nifti_1_header)(const char *, const nifti_1_header *) = NULL;
static int(*_disp_nifti_2_header)(const char *, const nifti_2_header *) = NULL;
static void(*_nifti_set_debug_level)(int) = NULL;
static void(*_nifti_set_skip_blank_ext)(int) = NULL;
static void(*_nifti_set_allow_upper_fext)(int) = NULL;
static int(*_nifti_get_alter_cifti)(void) = NULL;
static void(*_nifti_set_alter_cifti)(int) = NULL;
static int(*_nifti_alter_cifti_dims)(nifti_image *) = NULL;
static int(*_valid_nifti2_brick_list)(nifti_image *, int64_t, const int64_t *, int) = NULL;
static znzFile(*_nifti2_image_open)(const char *, char *, nifti_image **) = NULL;
static znzFile(*_nifti2_image_write_hdr_img)(nifti_image *, int, const char*) = NULL;
static znzFile(*_nifti2_image_write_hdr_img2)(nifti_image *, int, const char*, znzFile, const nifti_brick_list *) = NULL;
static int64_t(*_nifti2_read_buffer)(znzFile, void*, int64_t, nifti_image *) = NULL;
static int(*_nifti2_write_all_data)(znzFile, nifti_image *, const nifti_brick_list *) = NULL;
static int64_t(*_nifti2_write_buffer)(znzFile, const void *, int64_t) = NULL;
static nifti_image *(*_nifti2_read_ascii_image)(znzFile, const char *, int, int) = NULL;
static znzFile(*_nifti2_write_ascii_image)(nifti_image *, const nifti_brick_list *, const char *, int, int) = NULL;
static void(*_nifti_datatype_sizes)(int, int *, int *) = NULL;
static void(*_nifti_dmat44_to_quatern)(nifti_dmat44, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *) = NULL;
static nifti_dmat44(*_nifti_quatern_to_dmat44)(double, double, double, double, double, double, double, double, double, double) = NULL;
static nifti_dmat44(*_nifti_make_orthog_dmat44)(double, double, double, double, double, double, double, double, double) = NULL;
static void(*_nifti_mat44_to_quatern)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *) = NULL;
static mat44(*_nifti_quatern_to_mat44)(float, float, float, float, float, float, float, float, float, float) = NULL;
static mat44(*_nifti_make_orthog_mat44)(float, float, float, float, float, float, float, float, float) = NULL;
static int(*_nifti_short_order)(void) = NULL;
static void(*_nifti_mat44_to_orientation)(mat44, int *, int *, int *) = NULL;
static void(*_nifti_dmat44_to_orientation)(nifti_dmat44, int *, int *, int *) = NULL;
static char *(*_nifti_findhdrname)(const char*) = NULL;
static char *(*_nifti_findimgname)(const char*, int) = NULL;
static int(*_nifti_is_gzfile)(const char*) = NULL;
static char *(*_nifti_makebasename)(const char*) = NULL;
static int(*_nifti_convert_nim2n1hdr)(const nifti_image*, nifti_1_header *) = NULL;
static int(*_nifti_convert_nim2n2hdr)(const nifti_image*, nifti_2_header *) = NULL;
static nifti_1_header *(*_nifti_make_new_n1_header)(const int64_t[], int) = NULL;
static nifti_2_header *(*_nifti_make_new_n2_header)(const int64_t[], int) = NULL;
static void *(*_nifti2_read_header)(const char *, int *, int) = NULL;
static nifti_1_header *(*_nifti_read_n1_hdr)(const char *, int *, int) = NULL;
static nifti_2_header *(*_nifti_read_n2_hdr)(const char *, int *, int) = NULL;
static nifti_image *(*_nifti2_copy_nim_info)(const nifti_image *) = NULL;
static nifti_image *(*_nifti2_make_new_nim)(const int64_t[], int, int) = NULL;
static nifti_image *(*_nifti2_simple_init_nim)(void) = NULL;
static nifti_image *(*_nifti_convert_n1hdr2nim)(nifti_1_header, const char *) = NULL;
static nifti_image *(*_nifti_convert_n2hdr2nim)(nifti_2_header, const char *) = NULL;
static int(*_nifti_looks_like_cifti)(nifti_image *) = NULL;
static int(*_nifti_hdr1_looks_good)(const nifti_1_header *) = NULL;
static int(*_nifti_hdr2_looks_good)(const nifti_2_header *) = NULL;
static int(*_nifti_is_valid_datatype)(int) = NULL;
static int(*_nifti_is_valid_ecode)(int) = NULL;
static int(*_nifti2_nim_is_valid)(nifti_image *, int) = NULL;
static int(*_nifti2_nim_has_valid_dims)(nifti_image *, int) = NULL;
static int(*_is_valid_nifti2_type)(int) = NULL;
static int(*_nifti_test_datatype_sizes)(int) = NULL;
static int(*_nifti2_type_and_names_match)(nifti_image *, int) = NULL;
static int(*_nifti2_update_dims_from_array)(nifti_image *) = NULL;
static void(*_nifti2_set_iname_offset)(nifti_image *, int) = NULL;
static int(*_nifti2_set_type_from_names)(nifti_image *) = NULL;
static int(*_nifti2_add_extension)(nifti_image *, const char *, int, int) = NULL;
static int(*_nifti_compiled_with_zlib)(void) = NULL;
static int(*_nifti2_copy_extensions)(nifti_image *, const nifti_image *) = NULL;
static int(*_nifti2_free_extensions)(nifti_image *) = NULL;
static int64_t *(*_nifti_get_int64list)(int64_t, const char *) = NULL;
static int *(*_nifti_get_intlist)(int, const char *) = NULL;
static char *(*_nifti_strdup)(const char *) = NULL;
static int(*_valid_nifti2_extensions)(const nifti_image *) = NULL;
static int(*_nifti_valid_header_size)(int, int) = NULL;
#endif

void niftilib_register_all (void)
{
#ifdef _OPENMP
#pragma omp critical
#endif
    if (!registered)
    {
#if RNIFTI_NIFTILIB_VERSION == 1
        _nifti_datatype_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_datatype_string");
        _nifti_units_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_units_string");
        _nifti_intent_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_intent_string");
        _nifti_xform_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_xform_string");
        _nifti_slice_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_slice_string");
        _nifti_orientation_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_orientation_string");
        _nifti_is_inttype = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_inttype");
        _nifti_mat44_inverse = (mat44(*)(mat44)) R_GetCCallable("RNifti", "nii_mat44_inverse");
        _nifti_mat33_inverse = (mat33(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_inverse");
        _nifti_mat33_polar = (mat33(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_polar");
        _nifti_mat33_rownorm = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_rownorm");
        _nifti_mat33_colnorm = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_colnorm");
        _nifti_mat33_determ = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_determ");
        _nifti_mat33_mul = (mat33(*)(mat33, mat33)) R_GetCCallable("RNifti", "nii_mat33_mul");
        _nifti_swap_2bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti", "nii_swap_2bytes");
        _nifti_swap_4bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti", "nii_swap_4bytes");
        _nifti_swap_8bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti", "nii_swap_8bytes");
        _nifti_swap_16bytes = (void(*)(size_t, void *)) R_GetCCallable("RNifti", "nii_swap_16bytes");
        _nifti_swap_Nbytes = (void(*)(size_t, int, void *)) R_GetCCallable("RNifti", "nii_swap_Nbytes");
        _nifti_datatype_is_valid = (int(*)(int, int)) R_GetCCallable("RNifti", "nii_datatype_is_valid");
        _nifti_datatype_from_string = (int(*)(const char *)) R_GetCCallable("RNifti", "nii_datatype_from_string");
        _nifti_datatype_to_string = (const char *(*)(int)) R_GetCCallable("RNifti", "nii_datatype_to_string");
        _nifti_get_filesize = (int(*)(const char *)) R_GetCCallable("RNifti", "nii_get_filesize");
        _swap_nifti_header = (void(*)(struct nifti_1_header *, int)) R_GetCCallable("RNifti", "swap_nii_header");
        _old_swap_nifti_header = (void(*)(struct nifti_1_header *, int)) R_GetCCallable("RNifti", "old_swap_nii_header");
        _nifti_swap_as_analyze = (int(*)(nifti_analyze75 *)) R_GetCCallable("RNifti", "nii_swap_as_analyze");
        _nifti_image_read_bricks = (nifti_image *(*)(const char *, int, const int *, nifti_brick_list *)) R_GetCCallable("RNifti", "nii_image_read_bricks");
        _nifti_image_load_bricks = (int(*)(nifti_image *, int, const int *, nifti_brick_list *)) R_GetCCallable("RNifti", "nii_image_load_bricks");
        _nifti_free_NBL = (void(*)(nifti_brick_list *)) R_GetCCallable("RNifti", "nii_free_NBL");
        _nifti_image_read = (nifti_image *(*)(const char *, int)) R_GetCCallable("RNifti", "nii_image_read");
        _nifti_image_load = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_image_load");
        _nifti_image_unload = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_image_unload");
        _nifti_image_free = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_image_free");
        _nifti_read_collapsed_image = (int(*)(nifti_image *, const int[8], void **)) R_GetCCallable("RNifti", "nii_read_collapsed_image");
        _nifti_read_subregion_image = (int(*)(nifti_image *, const int *, const int *, void **)) R_GetCCallable("RNifti", "nii_read_subregion_image");
        _nifti_image_write = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_image_write");
        _nifti_image_write_bricks = (void(*)(nifti_image *, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii_image_write_bricks");
        _nifti_image_infodump = (void(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii_image_infodump");
        _nifti_disp_lib_hist = (void(*)(void)) R_GetCCallable("RNifti", "nii_disp_lib_hist");
        _nifti_disp_lib_version = (void(*)(void)) R_GetCCallable("RNifti", "nii_disp_lib_version");
        _nifti_disp_matrix_orient = (int(*)(const char *, mat44)) R_GetCCallable("RNifti", "nii_disp_matrix_orient");
        _nifti_disp_type_list = (int(*)(int)) R_GetCCallable("RNifti", "nii_disp_type_list");
        _nifti_image_to_ascii = (char *(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii_image_to_ascii");
        _nifti_image_from_ascii = (nifti_image *(*)(const char *, int *)) R_GetCCallable("RNifti", "nii_image_from_ascii");
        _nifti_get_volsize = (size_t(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii_get_volsize");
        _nifti_set_filenames = (int(*)(nifti_image *, const char *, int, int)) R_GetCCallable("RNifti", "nii_set_filenames");
        _nifti_makehdrname = (char *(*)(const char *, int, int, int)) R_GetCCallable("RNifti", "nii_makehdrname");
        _nifti_makeimgname = (char *(*)(const char *, int, int, int)) R_GetCCallable("RNifti", "nii_makeimgname");
        _is_nifti_file = (int(*)(const char *)) R_GetCCallable("RNifti", "is_nii_file");
        _nifti_find_file_extension = (char *(*)(const char *)) R_GetCCallable("RNifti", "nii_find_file_extension");
        _nifti_is_complete_filename = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_is_complete_filename");
        _nifti_validfilename = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_validfilename");
        _disp_nifti_1_header = (int(*)(const char *, const nifti_1_header *)) R_GetCCallable("RNifti", "disp_nii_1_header");
        _nifti_set_debug_level = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_debug_level");
        _nifti_set_skip_blank_ext = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_skip_blank_ext");
        _nifti_set_allow_upper_fext = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_allow_upper_fext");
        _valid_nifti_brick_list = (int(*)(nifti_image *, int, const int *, int)) R_GetCCallable("RNifti", "valid_nii_brick_list");
        _nifti_image_open = (znzFile(*)(const char *, const char *, nifti_image **)) R_GetCCallable("RNifti", "nii_image_open");
        _nifti_image_write_hdr_img = (znzFile(*)(nifti_image *, int, const char*)) R_GetCCallable("RNifti", "nii_image_write_hdr_img");
        _nifti_image_write_hdr_img2 = (znzFile(*)(nifti_image *, int, const char*, znzFile, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii_image_write_hdr_img2");
        _nifti_read_buffer = (size_t(*)(znzFile, void*, size_t, nifti_image *)) R_GetCCallable("RNifti", "nii_read_buffer");
        _nifti_write_all_data = (int(*)(znzFile, nifti_image *, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii_write_all_data");
        _nifti_write_buffer = (size_t(*)(znzFile, const void *, size_t)) R_GetCCallable("RNifti", "nii_write_buffer");
        _nifti_read_ascii_image = (nifti_image *(*)(znzFile, char *, int, int)) R_GetCCallable("RNifti", "nii_read_ascii_image");
        _nifti_write_ascii_image = (znzFile(*)(nifti_image *, const nifti_brick_list *, const char *, int, int)) R_GetCCallable("RNifti", "nii_write_ascii_image");
        _nifti_datatype_sizes = (void(*)(int, int *, int *)) R_GetCCallable("RNifti", "nii_datatype_sizes");
        _nifti_mat44_to_quatern = (void(*)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *)) R_GetCCallable("RNifti", "nii_mat44_to_quatern");
        _nifti_quatern_to_mat44 = (mat44(*)(float, float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti", "nii_quatern_to_mat44");
        _nifti_make_orthog_mat44 = (mat44(*)(float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti", "nii_make_orthog_mat44");
        _nifti_short_order = (int(*)(void)) R_GetCCallable("RNifti", "nii_short_order");
        _nifti_mat44_to_orientation = (void(*)(mat44, int *, int *, int *)) R_GetCCallable("RNifti", "nii_mat44_to_orientation");
        _nifti_findhdrname = (char *(*)(const char*)) R_GetCCallable("RNifti", "nii_findhdrname");
        _nifti_findimgname = (char *(*)(const char*, int)) R_GetCCallable("RNifti", "nii_findimgname");
        _nifti_is_gzfile = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_is_gzfile");
        _nifti_makebasename = (char *(*)(const char*)) R_GetCCallable("RNifti", "nii_makebasename");
        _nifti_convert_nim2nhdr = (struct nifti_1_header(*)(const nifti_image*)) R_GetCCallable("RNifti", "nii_convert_nim2nhdr");
        _nifti_make_new_header = (nifti_1_header *(*)(const int[], int)) R_GetCCallable("RNifti", "nii_make_new_header");
        _nifti_read_header = (nifti_1_header *(*)(const char *, int *, int)) R_GetCCallable("RNifti", "nii_read_header");
        _nifti_copy_nim_info = (nifti_image *(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii_copy_nim_info");
        _nifti_make_new_nim = (nifti_image *(*)(const int[], int, int)) R_GetCCallable("RNifti", "nii_make_new_nim");
        _nifti_simple_init_nim = (nifti_image *(*)(void)) R_GetCCallable("RNifti", "nii_simple_init_nim");
        _nifti_convert_nhdr2nim = (nifti_image *(*)(struct nifti_1_header, const char *)) R_GetCCallable("RNifti", "nii_convert_nhdr2nim");
        _nifti_hdr_looks_good = (int(*)(const nifti_1_header *)) R_GetCCallable("RNifti", "nii_hdr_looks_good");
        _nifti_is_valid_datatype = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_valid_datatype");
        _nifti_is_valid_ecode = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_valid_ecode");
        _nifti_nim_is_valid = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii_nim_is_valid");
        _nifti_nim_has_valid_dims = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii_nim_has_valid_dims");
        _is_valid_nifti_type = (int(*)(int)) R_GetCCallable("RNifti", "is_valid_nii_type");
        _nifti_test_datatype_sizes = (int(*)(int)) R_GetCCallable("RNifti", "nii_test_datatype_sizes");
        _nifti_type_and_names_match = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii_type_and_names_match");
        _nifti_update_dims_from_array = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_update_dims_from_array");
        _nifti_set_iname_offset = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_set_iname_offset");
        _nifti_set_type_from_names = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_set_type_from_names");
        _nifti_add_extension = (int(*)(nifti_image *, const char *, int, int)) R_GetCCallable("RNifti", "nii_add_extension");
        _nifti_compiled_with_zlib = (int(*)(void)) R_GetCCallable("RNifti", "nii_compiled_with_zlib");
        _nifti_copy_extensions = (int(*)(nifti_image *, const nifti_image *)) R_GetCCallable("RNifti", "nii_copy_extensions");
        _nifti_free_extensions = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_free_extensions");
        _nifti_get_intlist = (int *(*)(int, const char *)) R_GetCCallable("RNifti", "nii_get_intlist");
        _nifti_strdup = (char *(*)(const char *)) R_GetCCallable("RNifti", "nii_strdup");
        _valid_nifti_extensions = (int(*)(const nifti_image *)) R_GetCCallable("RNifti", "valid_nii_extensions");
#elif RNIFTI_NIFTILIB_VERSION == 2
        _nifti_datatype_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_datatype_string");
        _nifti_units_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_units_string");
        _nifti_intent_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_intent_string");
        _nifti_xform_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_xform_string");
        _nifti_slice_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_slice_string");
        _nifti_orientation_string = (char const *(*)(int)) R_GetCCallable("RNifti", "nii_orientation_string");
        _nifti_is_inttype = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_inttype");
        _nifti_mat44_inverse = (mat44(*)(mat44)) R_GetCCallable("RNifti", "nii_mat44_inverse");
        _nifti_mat44_mul = (mat44(*)(mat44, mat44)) R_GetCCallable("RNifti", "nii_mat44_mul");
        _nifti_dmat44_inverse = (nifti_dmat44(*)(nifti_dmat44)) R_GetCCallable("RNifti", "nii_dmat44_inverse");
        _nifti_mat44_to_dmat44 = (int(*)(mat44 *, nifti_dmat44 *)) R_GetCCallable("RNifti", "nii_mat44_to_dmat44");
        _nifti_dmat44_to_mat44 = (int(*)(nifti_dmat44 *, mat44 *)) R_GetCCallable("RNifti", "nii_dmat44_to_mat44");
        _nifti_dmat44_mul = (nifti_dmat44(*)(nifti_dmat44, nifti_dmat44)) R_GetCCallable("RNifti", "nii_dmat44_mul");
        _nifti_dmat33_inverse = (nifti_dmat33(*)(nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_inverse");
        _nifti_dmat33_polar = (nifti_dmat33(*)(nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_polar");
        _nifti_dmat33_rownorm = (double(*)(nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_rownorm");
        _nifti_dmat33_colnorm = (double(*)(nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_colnorm");
        _nifti_dmat33_determ = (double(*)(nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_determ");
        _nifti_dmat33_mul = (nifti_dmat33(*)(nifti_dmat33, nifti_dmat33)) R_GetCCallable("RNifti", "nii_dmat33_mul");
        _nifti_mat33_inverse = (mat33(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_inverse");
        _nifti_mat33_polar = (mat33(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_polar");
        _nifti_mat33_rownorm = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_rownorm");
        _nifti_mat33_colnorm = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_colnorm");
        _nifti_mat33_determ = (float(*)(mat33)) R_GetCCallable("RNifti", "nii_mat33_determ");
        _nifti_mat33_mul = (mat33(*)(mat33, mat33)) R_GetCCallable("RNifti", "nii_mat33_mul");
        _nifti_swap_2bytes = (void(*)(int64_t, void *)) R_GetCCallable("RNifti", "nii_swap_2bytes");
        _nifti_swap_4bytes = (void(*)(int64_t, void *)) R_GetCCallable("RNifti", "nii_swap_4bytes");
        _nifti_swap_8bytes = (void(*)(int64_t, void *)) R_GetCCallable("RNifti", "nii_swap_8bytes");
        _nifti_swap_16bytes = (void(*)(int64_t, void *)) R_GetCCallable("RNifti", "nii_swap_16bytes");
        _nifti_swap_Nbytes = (void(*)(int64_t, int, void *)) R_GetCCallable("RNifti", "nii_swap_Nbytes");
        _nifti_datatype_is_valid = (int(*)(int, int)) R_GetCCallable("RNifti", "nii_datatype_is_valid");
        _nifti_datatype_from_string = (int(*)(const char *)) R_GetCCallable("RNifti", "nii_datatype_from_string");
        _nifti_datatype_to_string = (const char *(*)(int)) R_GetCCallable("RNifti", "nii_datatype_to_string");
        _nifti_header_version = (int(*)(const char *, size_t)) R_GetCCallable("RNifti", "nii_header_version");
        _nifti2_get_filesize = (int64_t(*)(const char *)) R_GetCCallable("RNifti", "nii2_get_filesize");
        _swap_nifti_header = (void(*)(void *, int)) R_GetCCallable("RNifti", "swap_nii_header");
        _old_swap_nifti_header = (void(*)(struct nifti_1_header *, int)) R_GetCCallable("RNifti", "old_swap_nii_header");
        _nifti_swap_as_analyze = (void(*)(nifti_analyze75 *)) R_GetCCallable("RNifti", "nii_swap_as_analyze");
        _nifti_swap_as_nifti1 = (void(*)(nifti_1_header *)) R_GetCCallable("RNifti", "nii_swap_as_nifti1");
        _nifti_swap_as_nifti2 = (void(*)(nifti_2_header *)) R_GetCCallable("RNifti", "nii_swap_as_nifti2");
        _nifti2_image_read_bricks = (nifti_image *(*)(const char *, int64_t, const int64_t *, nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_image_read_bricks");
        _nifti2_image_load_bricks = (int(*)(nifti_image *, int64_t, const int64_t *, nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_image_load_bricks");
        _nifti2_free_NBL = (void(*)(nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_free_NBL");
        _nifti2_image_read = (nifti_image *(*)(const char *, int)) R_GetCCallable("RNifti", "nii2_image_read");
        _nifti2_image_load = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_image_load");
        _nifti2_image_unload = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_image_unload");
        _nifti2_image_free = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_image_free");
        _nifti2_read_collapsed_image = (int64_t(*)(nifti_image *, const int64_t[8], void **)) R_GetCCallable("RNifti", "nii2_read_collapsed_image");
        _nifti2_read_subregion_image = (int64_t(*)(nifti_image *, const int64_t *, const int64_t *, void **)) R_GetCCallable("RNifti", "nii2_read_subregion_image");
        _nifti2_image_write = (void(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_image_write");
        _nifti2_image_write_bricks = (void(*)(nifti_image *, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_image_write_bricks");
        _nifti2_image_infodump = (void(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii2_image_infodump");
        _nifti2_disp_lib_hist = (void(*)(int)) R_GetCCallable("RNifti", "nii2_disp_lib_hist");
        _nifti_disp_lib_version = (void(*)(void)) R_GetCCallable("RNifti", "nii_disp_lib_version");
        _nifti2_disp_matrix_orient = (int(*)(const char *, nifti_dmat44)) R_GetCCallable("RNifti", "nii2_disp_matrix_orient");
        _nifti_disp_type_list = (int(*)(int)) R_GetCCallable("RNifti", "nii_disp_type_list");
        _nifti2_image_to_ascii = (char *(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii2_image_to_ascii");
        _nifti2_image_from_ascii = (nifti_image *(*)(const char *, int *)) R_GetCCallable("RNifti", "nii2_image_from_ascii");
        _nifti2_get_volsize = (int64_t(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii2_get_volsize");
        _nifti2_set_filenames = (int(*)(nifti_image *, const char *, int, int)) R_GetCCallable("RNifti", "nii2_set_filenames");
        _nifti_makehdrname = (char *(*)(const char *, int, int, int)) R_GetCCallable("RNifti", "nii_makehdrname");
        _nifti_makeimgname = (char *(*)(const char *, int, int, int)) R_GetCCallable("RNifti", "nii_makeimgname");
        _is_nifti_file = (int(*)(const char *)) R_GetCCallable("RNifti", "is_nii_file");
        _nifti_find_file_extension = (char *(*)(const char *)) R_GetCCallable("RNifti", "nii_find_file_extension");
        _nifti_is_complete_filename = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_is_complete_filename");
        _nifti_validfilename = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_validfilename");
        _disp_nifti_1_header = (int(*)(const char *, const nifti_1_header *)) R_GetCCallable("RNifti", "disp_nii_1_header");
        _disp_nifti_2_header = (int(*)(const char *, const nifti_2_header *)) R_GetCCallable("RNifti", "disp_nii_2_header");
        _nifti_set_debug_level = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_debug_level");
        _nifti_set_skip_blank_ext = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_skip_blank_ext");
        _nifti_set_allow_upper_fext = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_allow_upper_fext");
        _nifti_get_alter_cifti = (int(*)(void)) R_GetCCallable("RNifti", "nii_get_alter_cifti");
        _nifti_set_alter_cifti = (void(*)(int)) R_GetCCallable("RNifti", "nii_set_alter_cifti");
        _nifti_alter_cifti_dims = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_alter_cifti_dims");
        _valid_nifti2_brick_list = (int(*)(nifti_image *, int64_t, const int64_t *, int)) R_GetCCallable("RNifti", "valid_nii2_brick_list");
        _nifti2_image_open = (znzFile(*)(const char *, char *, nifti_image **)) R_GetCCallable("RNifti", "nii2_image_open");
        _nifti2_image_write_hdr_img = (znzFile(*)(nifti_image *, int, const char*)) R_GetCCallable("RNifti", "nii2_image_write_hdr_img");
        _nifti2_image_write_hdr_img2 = (znzFile(*)(nifti_image *, int, const char*, znzFile, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_image_write_hdr_img2");
        _nifti2_read_buffer = (int64_t(*)(znzFile, void*, int64_t, nifti_image *)) R_GetCCallable("RNifti", "nii2_read_buffer");
        _nifti2_write_all_data = (int(*)(znzFile, nifti_image *, const nifti_brick_list *)) R_GetCCallable("RNifti", "nii2_write_all_data");
        _nifti2_write_buffer = (int64_t(*)(znzFile, const void *, int64_t)) R_GetCCallable("RNifti", "nii2_write_buffer");
        _nifti2_read_ascii_image = (nifti_image *(*)(znzFile, const char *, int, int)) R_GetCCallable("RNifti", "nii2_read_ascii_image");
        _nifti2_write_ascii_image = (znzFile(*)(nifti_image *, const nifti_brick_list *, const char *, int, int)) R_GetCCallable("RNifti", "nii2_write_ascii_image");
        _nifti_datatype_sizes = (void(*)(int, int *, int *)) R_GetCCallable("RNifti", "nii_datatype_sizes");
        _nifti_dmat44_to_quatern = (void(*)(nifti_dmat44, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *)) R_GetCCallable("RNifti", "nii_dmat44_to_quatern");
        _nifti_quatern_to_dmat44 = (nifti_dmat44(*)(double, double, double, double, double, double, double, double, double, double)) R_GetCCallable("RNifti", "nii_quatern_to_dmat44");
        _nifti_make_orthog_dmat44 = (nifti_dmat44(*)(double, double, double, double, double, double, double, double, double)) R_GetCCallable("RNifti", "nii_make_orthog_dmat44");
        _nifti_mat44_to_quatern = (void(*)(mat44, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *)) R_GetCCallable("RNifti", "nii_mat44_to_quatern");
        _nifti_quatern_to_mat44 = (mat44(*)(float, float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti", "nii_quatern_to_mat44");
        _nifti_make_orthog_mat44 = (mat44(*)(float, float, float, float, float, float, float, float, float)) R_GetCCallable("RNifti", "nii_make_orthog_mat44");
        _nifti_short_order = (int(*)(void)) R_GetCCallable("RNifti", "nii_short_order");
        _nifti_mat44_to_orientation = (void(*)(mat44, int *, int *, int *)) R_GetCCallable("RNifti", "nii_mat44_to_orientation");
        _nifti_dmat44_to_orientation = (void(*)(nifti_dmat44, int *, int *, int *)) R_GetCCallable("RNifti", "nii_dmat44_to_orientation");
        _nifti_findhdrname = (char *(*)(const char*)) R_GetCCallable("RNifti", "nii_findhdrname");
        _nifti_findimgname = (char *(*)(const char*, int)) R_GetCCallable("RNifti", "nii_findimgname");
        _nifti_is_gzfile = (int(*)(const char*)) R_GetCCallable("RNifti", "nii_is_gzfile");
        _nifti_makebasename = (char *(*)(const char*)) R_GetCCallable("RNifti", "nii_makebasename");
        _nifti_convert_nim2n1hdr = (int(*)(const nifti_image*, nifti_1_header *)) R_GetCCallable("RNifti", "nii_convert_nim2n1hdr");
        _nifti_convert_nim2n2hdr = (int(*)(const nifti_image*, nifti_2_header *)) R_GetCCallable("RNifti", "nii_convert_nim2n2hdr");
        _nifti_make_new_n1_header = (nifti_1_header *(*)(const int64_t[], int)) R_GetCCallable("RNifti", "nii_make_new_n1_header");
        _nifti_make_new_n2_header = (nifti_2_header *(*)(const int64_t[], int)) R_GetCCallable("RNifti", "nii_make_new_n2_header");
        _nifti2_read_header = (void *(*)(const char *, int *, int)) R_GetCCallable("RNifti", "nii2_read_header");
        _nifti_read_n1_hdr = (nifti_1_header *(*)(const char *, int *, int)) R_GetCCallable("RNifti", "nii_read_n1_hdr");
        _nifti_read_n2_hdr = (nifti_2_header *(*)(const char *, int *, int)) R_GetCCallable("RNifti", "nii_read_n2_hdr");
        _nifti2_copy_nim_info = (nifti_image *(*)(const nifti_image *)) R_GetCCallable("RNifti", "nii2_copy_nim_info");
        _nifti2_make_new_nim = (nifti_image *(*)(const int64_t[], int, int)) R_GetCCallable("RNifti", "nii2_make_new_nim");
        _nifti2_simple_init_nim = (nifti_image *(*)(void)) R_GetCCallable("RNifti", "nii2_simple_init_nim");
        _nifti_convert_n1hdr2nim = (nifti_image *(*)(nifti_1_header, const char *)) R_GetCCallable("RNifti", "nii_convert_n1hdr2nim");
        _nifti_convert_n2hdr2nim = (nifti_image *(*)(nifti_2_header, const char *)) R_GetCCallable("RNifti", "nii_convert_n2hdr2nim");
        _nifti_looks_like_cifti = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii_looks_like_cifti");
        _nifti_hdr1_looks_good = (int(*)(const nifti_1_header *)) R_GetCCallable("RNifti", "nii_hdr1_looks_good");
        _nifti_hdr2_looks_good = (int(*)(const nifti_2_header *)) R_GetCCallable("RNifti", "nii_hdr2_looks_good");
        _nifti_is_valid_datatype = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_valid_datatype");
        _nifti_is_valid_ecode = (int(*)(int)) R_GetCCallable("RNifti", "nii_is_valid_ecode");
        _nifti2_nim_is_valid = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii2_nim_is_valid");
        _nifti2_nim_has_valid_dims = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii2_nim_has_valid_dims");
        _is_valid_nifti2_type = (int(*)(int)) R_GetCCallable("RNifti", "is_valid_nii2_type");
        _nifti_test_datatype_sizes = (int(*)(int)) R_GetCCallable("RNifti", "nii_test_datatype_sizes");
        _nifti2_type_and_names_match = (int(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii2_type_and_names_match");
        _nifti2_update_dims_from_array = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_update_dims_from_array");
        _nifti2_set_iname_offset = (void(*)(nifti_image *, int)) R_GetCCallable("RNifti", "nii2_set_iname_offset");
        _nifti2_set_type_from_names = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_set_type_from_names");
        _nifti2_add_extension = (int(*)(nifti_image *, const char *, int, int)) R_GetCCallable("RNifti", "nii2_add_extension");
        _nifti_compiled_with_zlib = (int(*)(void)) R_GetCCallable("RNifti", "nii_compiled_with_zlib");
        _nifti2_copy_extensions = (int(*)(nifti_image *, const nifti_image *)) R_GetCCallable("RNifti", "nii2_copy_extensions");
        _nifti2_free_extensions = (int(*)(nifti_image *)) R_GetCCallable("RNifti", "nii2_free_extensions");
        _nifti_get_int64list = (int64_t *(*)(int64_t, const char *)) R_GetCCallable("RNifti", "nii_get_int64list");
        _nifti_get_intlist = (int *(*)(int, const char *)) R_GetCCallable("RNifti", "nii_get_intlist");
        _nifti_strdup = (char *(*)(const char *)) R_GetCCallable("RNifti", "nii_strdup");
        _valid_nifti2_extensions = (int(*)(const nifti_image *)) R_GetCCallable("RNifti", "valid_nii2_extensions");
        _nifti_valid_header_size = (int(*)(int, int)) R_GetCCallable("RNifti", "nii_valid_header_size");
#endif
        
        registered = 1;
    }
}

#if RNIFTI_NIFTILIB_VERSION == 1
char const * nifti_datatype_string (int dt) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_string, dt) }
char const * nifti_units_string (int uu) { NIFTILIB_WRAPPER_BODY(_nifti_units_string, uu) }
char const * nifti_intent_string (int ii) { NIFTILIB_WRAPPER_BODY(_nifti_intent_string, ii) }
char const * nifti_xform_string (int xx) { NIFTILIB_WRAPPER_BODY(_nifti_xform_string, xx) }
char const * nifti_slice_string (int ss) { NIFTILIB_WRAPPER_BODY(_nifti_slice_string, ss) }
char const * nifti_orientation_string (int ii) { NIFTILIB_WRAPPER_BODY(_nifti_orientation_string, ii) }
int nifti_is_inttype (int dt) { NIFTILIB_WRAPPER_BODY(_nifti_is_inttype, dt) }
mat44 nifti_mat44_inverse (mat44 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat44_inverse, R) }
mat33 nifti_mat33_inverse (mat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_inverse, R) }
mat33 nifti_mat33_polar (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_polar, A) }
float nifti_mat33_rownorm (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_rownorm, A) }
float nifti_mat33_colnorm (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_colnorm, A) }
float nifti_mat33_determ (mat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_determ, R) }
mat33 nifti_mat33_mul (mat33 A, mat33 B) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_mul, A, B) }
void nifti_swap_2bytes (size_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_2bytes, n, ar) }
void nifti_swap_4bytes (size_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_4bytes, n, ar) }
void nifti_swap_8bytes (size_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_8bytes, n, ar) }
void nifti_swap_16bytes (size_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_16bytes, n, ar) }
void nifti_swap_Nbytes (size_t n, int siz, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_Nbytes, n, siz, ar) }
int nifti_datatype_is_valid (int dtype, int for_nifti) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_is_valid, dtype, for_nifti) }
int nifti_datatype_from_string (const char * name) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_from_string, name) }
const char * nifti_datatype_to_string (int dtype) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_to_string, dtype) }
int nifti_get_filesize (const char * pathname) { NIFTILIB_WRAPPER_BODY(_nifti_get_filesize, pathname) }
void swap_nifti_header (struct nifti_1_header * h, int is_nifti) { NIFTILIB_WRAPPER_BODY_NORETURN(_swap_nifti_header, h, is_nifti) }
void old_swap_nifti_header (struct nifti_1_header * h, int is_nifti) { NIFTILIB_WRAPPER_BODY_NORETURN(_old_swap_nifti_header, h, is_nifti) }
int nifti_swap_as_analyze (nifti_analyze75 * h) { NIFTILIB_WRAPPER_BODY(_nifti_swap_as_analyze, h) }
nifti_image * nifti_image_read_bricks (const char * hname, int nbricks, const int * blist, nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti_image_read_bricks, hname, nbricks, blist, NBL) }
int nifti_image_load_bricks (nifti_image * nim, int nbricks, const int * blist, nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti_image_load_bricks, nim, nbricks, blist, NBL) }
void nifti_free_NBL (nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_free_NBL, NBL) }
nifti_image * nifti_image_read (const char * hname, int read_data) { NIFTILIB_WRAPPER_BODY(_nifti_image_read, hname, read_data) }
int nifti_image_load (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_image_load, nim) }
void nifti_image_unload (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_image_unload, nim) }
void nifti_image_free (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_image_free, nim) }
int nifti_read_collapsed_image (nifti_image * nim, const int dims[8], void ** data) { NIFTILIB_WRAPPER_BODY(_nifti_read_collapsed_image, nim, dims, data) }
int nifti_read_subregion_image (nifti_image * nim, const int * start_index, const int * region_size, void ** data) { NIFTILIB_WRAPPER_BODY(_nifti_read_subregion_image, nim, start_index, region_size, data) }
void nifti_image_write (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_image_write, nim) }
void nifti_image_write_bricks (nifti_image * nim, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_image_write_bricks, nim, NBL) }
void nifti_image_infodump (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_image_infodump, nim) }
void nifti_disp_lib_hist (void) { NIFTILIB_WRAPPER_BODY_VOID_NORETURN(_nifti_disp_lib_hist) }
void nifti_disp_lib_version (void) { NIFTILIB_WRAPPER_BODY_VOID_NORETURN(_nifti_disp_lib_version) }
int nifti_disp_matrix_orient (const char * mesg, mat44 mat) { NIFTILIB_WRAPPER_BODY(_nifti_disp_matrix_orient, mesg, mat) }
int nifti_disp_type_list (int which) { NIFTILIB_WRAPPER_BODY(_nifti_disp_type_list, which) }
char * nifti_image_to_ascii (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_image_to_ascii, nim) }
nifti_image * nifti_image_from_ascii (const char * str, int * bytes_read) { NIFTILIB_WRAPPER_BODY(_nifti_image_from_ascii, str, bytes_read) }
size_t nifti_get_volsize (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_get_volsize, nim) }
int nifti_set_filenames (nifti_image * nim, const char * prefix, int check, int set_byte_order) { NIFTILIB_WRAPPER_BODY(_nifti_set_filenames, nim, prefix, check, set_byte_order) }
char * nifti_makehdrname (const char * prefix, int nifti_type, int check, int comp) { NIFTILIB_WRAPPER_BODY(_nifti_makehdrname, prefix, nifti_type, check, comp) }
char * nifti_makeimgname (const char * prefix, int nifti_type, int check, int comp) { NIFTILIB_WRAPPER_BODY(_nifti_makeimgname, prefix, nifti_type, check, comp) }
int is_nifti_file (const char * hname) { NIFTILIB_WRAPPER_BODY(_is_nifti_file, hname) }
char * nifti_find_file_extension (const char * name) { NIFTILIB_WRAPPER_BODY(_nifti_find_file_extension, name) }
int nifti_is_complete_filename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_is_complete_filename, fname) }
int nifti_validfilename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_validfilename, fname) }
int disp_nifti_1_header (const char * info, const nifti_1_header * hp) { NIFTILIB_WRAPPER_BODY(_disp_nifti_1_header, info, hp) }
void nifti_set_debug_level (int level) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_debug_level, level) }
void nifti_set_skip_blank_ext (int skip) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_skip_blank_ext, skip) }
void nifti_set_allow_upper_fext (int allow) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_allow_upper_fext, allow) }
int valid_nifti_brick_list (nifti_image * nim, int nbricks, const int * blist, int disp_error) { NIFTILIB_WRAPPER_BODY(_valid_nifti_brick_list, nim, nbricks, blist, disp_error) }
znzFile nifti_image_open (const char * hname, const char * opts, nifti_image ** nim) { NIFTILIB_WRAPPER_BODY(_nifti_image_open, hname, opts, nim) }
znzFile nifti_image_write_hdr_img (nifti_image * nim, int write_data, const char* opts) { NIFTILIB_WRAPPER_BODY(_nifti_image_write_hdr_img, nim, write_data, opts) }
znzFile nifti_image_write_hdr_img2 (nifti_image * nim, int write_opts, const char* opts, znzFile imgfile, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti_image_write_hdr_img2, nim, write_opts, opts, imgfile, NBL) }
size_t nifti_read_buffer (znzFile fp, void* dataptr, size_t ntot, nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_read_buffer, fp, dataptr, ntot, nim) }
int nifti_write_all_data (znzFile fp, nifti_image * nim, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti_write_all_data, fp, nim, NBL) }
size_t nifti_write_buffer (znzFile fp, const void * buffer, size_t numbytes) { NIFTILIB_WRAPPER_BODY(_nifti_write_buffer, fp, buffer, numbytes) }
nifti_image * nifti_read_ascii_image (znzFile fp, char * fname, int flen, int read_data) { NIFTILIB_WRAPPER_BODY(_nifti_read_ascii_image, fp, fname, flen, read_data) }
znzFile nifti_write_ascii_image (nifti_image * nim, const nifti_brick_list * NBL, const char * opts, int write_data, int leave_open) { NIFTILIB_WRAPPER_BODY(_nifti_write_ascii_image, nim, NBL, opts, write_data, leave_open) }
void nifti_datatype_sizes (int datatype, int * nbyper, int * swapsize) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_datatype_sizes, datatype, nbyper, swapsize) }
void nifti_mat44_to_quatern (mat44 R, float * qb, float * qc, float * qd, float * qx, float * qy, float * qz, float * dx, float * dy, float * dz, float * qfac) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_mat44_to_quatern, R, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
mat44 nifti_quatern_to_mat44 (float qb, float qc, float qd, float qx, float qy, float qz, float dx, float dy, float dz, float qfac) { NIFTILIB_WRAPPER_BODY(_nifti_quatern_to_mat44, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
mat44 nifti_make_orthog_mat44 (float r11, float r12, float r13, float r21, float r22, float r23, float r31, float r32, float r33) { NIFTILIB_WRAPPER_BODY(_nifti_make_orthog_mat44, r11, r12, r13, r21, r22, r23, r31, r32, r33) }
int nifti_short_order (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_short_order) }
void nifti_mat44_to_orientation (mat44 R, int * icod, int * jcod, int * kcod) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_mat44_to_orientation, R, icod, jcod, kcod) }
char * nifti_findhdrname (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_findhdrname, fname) }
char * nifti_findimgname (const char* fname, int nifti_type) { NIFTILIB_WRAPPER_BODY(_nifti_findimgname, fname, nifti_type) }
int nifti_is_gzfile (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_is_gzfile, fname) }
char * nifti_makebasename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_makebasename, fname) }
struct nifti_1_header nifti_convert_nim2nhdr (const nifti_image* nim) { NIFTILIB_WRAPPER_BODY(_nifti_convert_nim2nhdr, nim) }
nifti_1_header * nifti_make_new_header (const int arg_dims[], int arg_dtype) { NIFTILIB_WRAPPER_BODY(_nifti_make_new_header, arg_dims, arg_dtype) }
nifti_1_header * nifti_read_header (const char * hname, int * swapped, int check) { NIFTILIB_WRAPPER_BODY(_nifti_read_header, hname, swapped, check) }
nifti_image * nifti_copy_nim_info (const nifti_image * src) { NIFTILIB_WRAPPER_BODY(_nifti_copy_nim_info, src) }
nifti_image * nifti_make_new_nim (const int dims[], int datatype, int data_fill) { NIFTILIB_WRAPPER_BODY(_nifti_make_new_nim, dims, datatype, data_fill) }
nifti_image * nifti_simple_init_nim (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_simple_init_nim) }
nifti_image * nifti_convert_nhdr2nim (struct nifti_1_header nhdr, const char * fname) { NIFTILIB_WRAPPER_BODY(_nifti_convert_nhdr2nim, nhdr, fname) }
int nifti_hdr_looks_good (const nifti_1_header * hdr) { NIFTILIB_WRAPPER_BODY(_nifti_hdr_looks_good, hdr) }
int nifti_is_valid_datatype (int dtype) { NIFTILIB_WRAPPER_BODY(_nifti_is_valid_datatype, dtype) }
int nifti_is_valid_ecode (int ecode) { NIFTILIB_WRAPPER_BODY(_nifti_is_valid_ecode, ecode) }
int nifti_nim_is_valid (nifti_image * nim, int complain) { NIFTILIB_WRAPPER_BODY(_nifti_nim_is_valid, nim, complain) }
int nifti_nim_has_valid_dims (nifti_image * nim, int complain) { NIFTILIB_WRAPPER_BODY(_nifti_nim_has_valid_dims, nim, complain) }
int is_valid_nifti_type (int nifti_type) { NIFTILIB_WRAPPER_BODY(_is_valid_nifti_type, nifti_type) }
int nifti_test_datatype_sizes (int verb) { NIFTILIB_WRAPPER_BODY(_nifti_test_datatype_sizes, verb) }
int nifti_type_and_names_match (nifti_image * nim, int show_warn) { NIFTILIB_WRAPPER_BODY(_nifti_type_and_names_match, nim, show_warn) }
int nifti_update_dims_from_array (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_update_dims_from_array, nim) }
void nifti_set_iname_offset (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_iname_offset, nim) }
int nifti_set_type_from_names (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_set_type_from_names, nim) }
int nifti_add_extension (nifti_image * nim, const char * data, int len, int ecode) { NIFTILIB_WRAPPER_BODY(_nifti_add_extension, nim, data, len, ecode) }
int nifti_compiled_with_zlib (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_compiled_with_zlib) }
int nifti_copy_extensions (nifti_image * nim_dest, const nifti_image * nim_src) { NIFTILIB_WRAPPER_BODY(_nifti_copy_extensions, nim_dest, nim_src) }
int nifti_free_extensions (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_free_extensions, nim) }
int * nifti_get_intlist (int nvals, const char * str) { NIFTILIB_WRAPPER_BODY(_nifti_get_intlist, nvals, str) }
char * nifti_strdup (const char * str) { NIFTILIB_WRAPPER_BODY(_nifti_strdup, str) }
int valid_nifti_extensions (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_valid_nifti_extensions, nim) }
#elif RNIFTI_NIFTILIB_VERSION == 2
char const * nifti_datatype_string (int dt) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_string, dt) }
char const * nifti_units_string (int uu) { NIFTILIB_WRAPPER_BODY(_nifti_units_string, uu) }
char const * nifti_intent_string (int ii) { NIFTILIB_WRAPPER_BODY(_nifti_intent_string, ii) }
char const * nifti_xform_string (int xx) { NIFTILIB_WRAPPER_BODY(_nifti_xform_string, xx) }
char const * nifti_slice_string (int ss) { NIFTILIB_WRAPPER_BODY(_nifti_slice_string, ss) }
char const * nifti_orientation_string (int ii) { NIFTILIB_WRAPPER_BODY(_nifti_orientation_string, ii) }
int nifti_is_inttype (int dt) { NIFTILIB_WRAPPER_BODY(_nifti_is_inttype, dt) }
mat44 nifti_mat44_inverse (mat44 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat44_inverse, R) }
mat44 nifti_mat44_mul (mat44 A, mat44 B) { NIFTILIB_WRAPPER_BODY(_nifti_mat44_mul, A, B) }
nifti_dmat44 nifti_dmat44_inverse (nifti_dmat44 R) { NIFTILIB_WRAPPER_BODY(_nifti_dmat44_inverse, R) }
int nifti_mat44_to_dmat44 (mat44 * fm, nifti_dmat44 * dm) { NIFTILIB_WRAPPER_BODY(_nifti_mat44_to_dmat44, fm, dm) }
int nifti_dmat44_to_mat44 (nifti_dmat44 * dm, mat44 * fm) { NIFTILIB_WRAPPER_BODY(_nifti_dmat44_to_mat44, dm, fm) }
nifti_dmat44 nifti_dmat44_mul (nifti_dmat44 A, nifti_dmat44 B) { NIFTILIB_WRAPPER_BODY(_nifti_dmat44_mul, A, B) }
nifti_dmat33 nifti_dmat33_inverse (nifti_dmat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_inverse, R) }
nifti_dmat33 nifti_dmat33_polar (nifti_dmat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_polar, A) }
double nifti_dmat33_rownorm (nifti_dmat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_rownorm, A) }
double nifti_dmat33_colnorm (nifti_dmat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_colnorm, A) }
double nifti_dmat33_determ (nifti_dmat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_determ, R) }
nifti_dmat33 nifti_dmat33_mul (nifti_dmat33 A, nifti_dmat33 B) { NIFTILIB_WRAPPER_BODY(_nifti_dmat33_mul, A, B) }
mat33 nifti_mat33_inverse (mat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_inverse, R) }
mat33 nifti_mat33_polar (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_polar, A) }
float nifti_mat33_rownorm (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_rownorm, A) }
float nifti_mat33_colnorm (mat33 A) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_colnorm, A) }
float nifti_mat33_determ (mat33 R) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_determ, R) }
mat33 nifti_mat33_mul (mat33 A, mat33 B) { NIFTILIB_WRAPPER_BODY(_nifti_mat33_mul, A, B) }
void nifti_swap_2bytes (int64_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_2bytes, n, ar) }
void nifti_swap_4bytes (int64_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_4bytes, n, ar) }
void nifti_swap_8bytes (int64_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_8bytes, n, ar) }
void nifti_swap_16bytes (int64_t n, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_16bytes, n, ar) }
void nifti_swap_Nbytes (int64_t n, int siz, void * ar) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_Nbytes, n, siz, ar) }
int nifti_datatype_is_valid (int dtype, int for_nifti) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_is_valid, dtype, for_nifti) }
int nifti_datatype_from_string (const char * name) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_from_string, name) }
const char * nifti_datatype_to_string (int dtype) { NIFTILIB_WRAPPER_BODY(_nifti_datatype_to_string, dtype) }
int nifti_header_version (const char * buf, size_t nbytes) { NIFTILIB_WRAPPER_BODY(_nifti_header_version, buf, nbytes) }
int64_t nifti2_get_filesize (const char * pathname) { NIFTILIB_WRAPPER_BODY(_nifti2_get_filesize, pathname) }
void swap_nifti_header (void * hdr, int ni_ver) { NIFTILIB_WRAPPER_BODY_NORETURN(_swap_nifti_header, hdr, ni_ver) }
void old_swap_nifti_header (struct nifti_1_header * h, int is_nifti) { NIFTILIB_WRAPPER_BODY_NORETURN(_old_swap_nifti_header, h, is_nifti) }
void nifti_swap_as_analyze (nifti_analyze75 * h) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_as_analyze, h) }
void nifti_swap_as_nifti1 (nifti_1_header * h) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_as_nifti1, h) }
void nifti_swap_as_nifti2 (nifti_2_header * h) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_swap_as_nifti2, h) }
nifti_image * nifti2_image_read_bricks (const char * hname, int64_t nbricks, const int64_t * blist, nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti2_image_read_bricks, hname, nbricks, blist, NBL) }
int nifti2_image_load_bricks (nifti_image * nim, int64_t nbricks, const int64_t * blist, nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti2_image_load_bricks, nim, nbricks, blist, NBL) }
void nifti2_free_NBL (nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_free_NBL, NBL) }
nifti_image * nifti2_image_read (const char * hname, int read_data) { NIFTILIB_WRAPPER_BODY(_nifti2_image_read, hname, read_data) }
int nifti2_image_load (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_image_load, nim) }
void nifti2_image_unload (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_image_unload, nim) }
void nifti2_image_free (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_image_free, nim) }
int64_t nifti2_read_collapsed_image (nifti_image * nim, const int64_t dims[8], void ** data) { NIFTILIB_WRAPPER_BODY(_nifti2_read_collapsed_image, nim, dims, data) }
int64_t nifti2_read_subregion_image (nifti_image * nim, const int64_t * start_index, const int64_t * region_size, void ** data) { NIFTILIB_WRAPPER_BODY(_nifti2_read_subregion_image, nim, start_index, region_size, data) }
void nifti2_image_write (nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_image_write, nim) }
void nifti2_image_write_bricks (nifti_image * nim, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_image_write_bricks, nim, NBL) }
void nifti2_image_infodump (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_image_infodump, nim) }
void nifti2_disp_lib_hist (int ver) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_disp_lib_hist, ver) }
void nifti_disp_lib_version (void) { NIFTILIB_WRAPPER_BODY_VOID_NORETURN(_nifti_disp_lib_version) }
int nifti2_disp_matrix_orient (const char * mesg, nifti_dmat44 mat) { NIFTILIB_WRAPPER_BODY(_nifti2_disp_matrix_orient, mesg, mat) }
int nifti_disp_type_list (int which) { NIFTILIB_WRAPPER_BODY(_nifti_disp_type_list, which) }
char * nifti2_image_to_ascii (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_image_to_ascii, nim) }
nifti_image * nifti2_image_from_ascii (const char * str, int * bytes_read) { NIFTILIB_WRAPPER_BODY(_nifti2_image_from_ascii, str, bytes_read) }
int64_t nifti2_get_volsize (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_get_volsize, nim) }
int nifti2_set_filenames (nifti_image * nim, const char * prefix, int check, int set_byte_order) { NIFTILIB_WRAPPER_BODY(_nifti2_set_filenames, nim, prefix, check, set_byte_order) }
char * nifti_makehdrname (const char * prefix, int nifti_type, int check, int comp) { NIFTILIB_WRAPPER_BODY(_nifti_makehdrname, prefix, nifti_type, check, comp) }
char * nifti_makeimgname (const char * prefix, int nifti_type, int check, int comp) { NIFTILIB_WRAPPER_BODY(_nifti_makeimgname, prefix, nifti_type, check, comp) }
int is_nifti_file (const char * hname) { NIFTILIB_WRAPPER_BODY(_is_nifti_file, hname) }
char * nifti_find_file_extension (const char * name) { NIFTILIB_WRAPPER_BODY(_nifti_find_file_extension, name) }
int nifti_is_complete_filename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_is_complete_filename, fname) }
int nifti_validfilename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_validfilename, fname) }
int disp_nifti_1_header (const char * info, const nifti_1_header * hp) { NIFTILIB_WRAPPER_BODY(_disp_nifti_1_header, info, hp) }
int disp_nifti_2_header (const char * info, const nifti_2_header * hp) { NIFTILIB_WRAPPER_BODY(_disp_nifti_2_header, info, hp) }
void nifti_set_debug_level (int level) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_debug_level, level) }
void nifti_set_skip_blank_ext (int skip) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_skip_blank_ext, skip) }
void nifti_set_allow_upper_fext (int allow) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_allow_upper_fext, allow) }
int nifti_get_alter_cifti (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_get_alter_cifti) }
void nifti_set_alter_cifti (int alter_cifti) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_set_alter_cifti, alter_cifti) }
int nifti_alter_cifti_dims (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_alter_cifti_dims, nim) }
int valid_nifti2_brick_list (nifti_image * nim, int64_t nbricks, const int64_t * blist, int disp_error) { NIFTILIB_WRAPPER_BODY(_valid_nifti2_brick_list, nim, nbricks, blist, disp_error) }
znzFile nifti2_image_open (const char * hname, char * opts, nifti_image ** nim) { NIFTILIB_WRAPPER_BODY(_nifti2_image_open, hname, opts, nim) }
znzFile nifti2_image_write_hdr_img (nifti_image * nim, int write_data, const char* opts) { NIFTILIB_WRAPPER_BODY(_nifti2_image_write_hdr_img, nim, write_data, opts) }
znzFile nifti2_image_write_hdr_img2 (nifti_image * nim, int write_opts, const char* opts, znzFile imgfile, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti2_image_write_hdr_img2, nim, write_opts, opts, imgfile, NBL) }
int64_t nifti2_read_buffer (znzFile fp, void* dataptr, int64_t ntot, nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_read_buffer, fp, dataptr, ntot, nim) }
int nifti2_write_all_data (znzFile fp, nifti_image * nim, const nifti_brick_list * NBL) { NIFTILIB_WRAPPER_BODY(_nifti2_write_all_data, fp, nim, NBL) }
int64_t nifti2_write_buffer (znzFile fp, const void * buffer, int64_t numbytes) { NIFTILIB_WRAPPER_BODY(_nifti2_write_buffer, fp, buffer, numbytes) }
nifti_image * nifti2_read_ascii_image (znzFile fp, const char * fname, int flen, int read_data) { NIFTILIB_WRAPPER_BODY(_nifti2_read_ascii_image, fp, fname, flen, read_data) }
znzFile nifti2_write_ascii_image (nifti_image * nim, const nifti_brick_list * NBL, const char * opts, int write_data, int leave_open) { NIFTILIB_WRAPPER_BODY(_nifti2_write_ascii_image, nim, NBL, opts, write_data, leave_open) }
void nifti_datatype_sizes (int datatype, int * nbyper, int * swapsize) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_datatype_sizes, datatype, nbyper, swapsize) }
void nifti_dmat44_to_quatern (nifti_dmat44 R, double * qb, double * qc, double * qd, double * qx, double * qy, double * qz, double * dx, double * dy, double * dz, double * qfac) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_dmat44_to_quatern, R, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
nifti_dmat44 nifti_quatern_to_dmat44 (double qb, double qc, double qd, double qx, double qy, double qz, double dx, double dy, double dz, double qfac) { NIFTILIB_WRAPPER_BODY(_nifti_quatern_to_dmat44, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
nifti_dmat44 nifti_make_orthog_dmat44 (double r11, double r12, double r13, double r21, double r22, double r23, double r31, double r32, double r33) { NIFTILIB_WRAPPER_BODY(_nifti_make_orthog_dmat44, r11, r12, r13, r21, r22, r23, r31, r32, r33) }
void nifti_mat44_to_quatern (mat44 R, float * qb, float * qc, float * qd, float * qx, float * qy, float * qz, float * dx, float * dy, float * dz, float * qfac) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_mat44_to_quatern, R, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
mat44 nifti_quatern_to_mat44 (float qb, float qc, float qd, float qx, float qy, float qz, float dx, float dy, float dz, float qfac) { NIFTILIB_WRAPPER_BODY(_nifti_quatern_to_mat44, qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac) }
mat44 nifti_make_orthog_mat44 (float r11, float r12, float r13, float r21, float r22, float r23, float r31, float r32, float r33) { NIFTILIB_WRAPPER_BODY(_nifti_make_orthog_mat44, r11, r12, r13, r21, r22, r23, r31, r32, r33) }
int nifti_short_order (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_short_order) }
void nifti_mat44_to_orientation (mat44 R, int * icod, int * jcod, int * kcod) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_mat44_to_orientation, R, icod, jcod, kcod) }
void nifti_dmat44_to_orientation (nifti_dmat44 R, int * icod, int * jcod, int * kcod) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti_dmat44_to_orientation, R, icod, jcod, kcod) }
char * nifti_findhdrname (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_findhdrname, fname) }
char * nifti_findimgname (const char* fname, int nifti_type) { NIFTILIB_WRAPPER_BODY(_nifti_findimgname, fname, nifti_type) }
int nifti_is_gzfile (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_is_gzfile, fname) }
char * nifti_makebasename (const char* fname) { NIFTILIB_WRAPPER_BODY(_nifti_makebasename, fname) }
int nifti_convert_nim2n1hdr (const nifti_image* nim, nifti_1_header * hdr) { NIFTILIB_WRAPPER_BODY(_nifti_convert_nim2n1hdr, nim, hdr) }
int nifti_convert_nim2n2hdr (const nifti_image* nim, nifti_2_header * hdr) { NIFTILIB_WRAPPER_BODY(_nifti_convert_nim2n2hdr, nim, hdr) }
nifti_1_header * nifti_make_new_n1_header (const int64_t arg_dims[], int arg_dtype) { NIFTILIB_WRAPPER_BODY(_nifti_make_new_n1_header, arg_dims, arg_dtype) }
nifti_2_header * nifti_make_new_n2_header (const int64_t arg_dims[], int arg_dtype) { NIFTILIB_WRAPPER_BODY(_nifti_make_new_n2_header, arg_dims, arg_dtype) }
void * nifti2_read_header (const char * hname, int * nver, int check) { NIFTILIB_WRAPPER_BODY(_nifti2_read_header, hname, nver, check) }
nifti_1_header * nifti_read_n1_hdr (const char * hname, int * swapped, int check) { NIFTILIB_WRAPPER_BODY(_nifti_read_n1_hdr, hname, swapped, check) }
nifti_2_header * nifti_read_n2_hdr (const char * hname, int * swapped, int check) { NIFTILIB_WRAPPER_BODY(_nifti_read_n2_hdr, hname, swapped, check) }
nifti_image * nifti2_copy_nim_info (const nifti_image * src) { NIFTILIB_WRAPPER_BODY(_nifti2_copy_nim_info, src) }
nifti_image * nifti2_make_new_nim (const int64_t dims[], int datatype, int data_fill) { NIFTILIB_WRAPPER_BODY(_nifti2_make_new_nim, dims, datatype, data_fill) }
nifti_image * nifti2_simple_init_nim (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti2_simple_init_nim) }
nifti_image * nifti_convert_n1hdr2nim (nifti_1_header nhdr, const char * fname) { NIFTILIB_WRAPPER_BODY(_nifti_convert_n1hdr2nim, nhdr, fname) }
nifti_image * nifti_convert_n2hdr2nim (nifti_2_header nhdr, const char * fname) { NIFTILIB_WRAPPER_BODY(_nifti_convert_n2hdr2nim, nhdr, fname) }
int nifti_looks_like_cifti (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti_looks_like_cifti, nim) }
int nifti_hdr1_looks_good (const nifti_1_header * hdr) { NIFTILIB_WRAPPER_BODY(_nifti_hdr1_looks_good, hdr) }
int nifti_hdr2_looks_good (const nifti_2_header * hdr) { NIFTILIB_WRAPPER_BODY(_nifti_hdr2_looks_good, hdr) }
int nifti_is_valid_datatype (int dtype) { NIFTILIB_WRAPPER_BODY(_nifti_is_valid_datatype, dtype) }
int nifti_is_valid_ecode (int ecode) { NIFTILIB_WRAPPER_BODY(_nifti_is_valid_ecode, ecode) }
int nifti2_nim_is_valid (nifti_image * nim, int complain) { NIFTILIB_WRAPPER_BODY(_nifti2_nim_is_valid, nim, complain) }
int nifti2_nim_has_valid_dims (nifti_image * nim, int complain) { NIFTILIB_WRAPPER_BODY(_nifti2_nim_has_valid_dims, nim, complain) }
int is_valid_nifti2_type (int nifti_type) { NIFTILIB_WRAPPER_BODY(_is_valid_nifti2_type, nifti_type) }
int nifti_test_datatype_sizes (int verb) { NIFTILIB_WRAPPER_BODY(_nifti_test_datatype_sizes, verb) }
int nifti2_type_and_names_match (nifti_image * nim, int show_warn) { NIFTILIB_WRAPPER_BODY(_nifti2_type_and_names_match, nim, show_warn) }
int nifti2_update_dims_from_array (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_update_dims_from_array, nim) }
void nifti2_set_iname_offset (nifti_image * nim, int nifti_ver) { NIFTILIB_WRAPPER_BODY_NORETURN(_nifti2_set_iname_offset, nim, nifti_ver) }
int nifti2_set_type_from_names (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_set_type_from_names, nim) }
int nifti2_add_extension (nifti_image * nim, const char * data, int len, int ecode) { NIFTILIB_WRAPPER_BODY(_nifti2_add_extension, nim, data, len, ecode) }
int nifti_compiled_with_zlib (void) { NIFTILIB_WRAPPER_BODY_VOID(_nifti_compiled_with_zlib) }
int nifti2_copy_extensions (nifti_image * nim_dest, const nifti_image * nim_src) { NIFTILIB_WRAPPER_BODY(_nifti2_copy_extensions, nim_dest, nim_src) }
int nifti2_free_extensions (nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_nifti2_free_extensions, nim) }
int64_t * nifti_get_int64list (int64_t nvals, const char * str) { NIFTILIB_WRAPPER_BODY(_nifti_get_int64list, nvals, str) }
int * nifti_get_intlist (int nvals, const char * str) { NIFTILIB_WRAPPER_BODY(_nifti_get_intlist, nvals, str) }
char * nifti_strdup (const char * str) { NIFTILIB_WRAPPER_BODY(_nifti_strdup, str) }
int valid_nifti2_extensions (const nifti_image * nim) { NIFTILIB_WRAPPER_BODY(_valid_nifti2_extensions, nim) }
int nifti_valid_header_size (int ni_ver, int whine) { NIFTILIB_WRAPPER_BODY(_nifti_valid_header_size, ni_ver, whine) }
#endif

#ifdef __cplusplus
} // extern "C"
#endif

#endif
