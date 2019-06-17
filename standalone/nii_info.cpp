// The NiftiImage_print.h header should be included first. The niftilib header
// nifti1_io.h contains prototypes for niftilib functions like
// nifti_convert_nim2nhdr and disp_nifti_1_header, but it is included by
// lib/NiftiImage.h, so is only here for explicitness
#include "lib/NiftiImage_print.h"
#include "niftilib/nifti1_io.h"
#include "lib/NiftiImage.h"

int main (int argc, char* argv[])
{
    // The program requires an image as an argument
    if (argc < 2)
    {
        // This macro is defined in NiftiImage_print.h. There is no obligation to use it
        Rc_fputs_stderr("Error: A NIfTI-1 image file must be specified\n");
        return 1;
    }
    
    int status = 0;
    
    try
    {
        // Create a string from the first argument
        const std::string path(argv[1]);
        // Read the image metadata (but not pixel data) from file
        RNifti::NiftiImage image(path, false);
        // Convert between struct types using niftilib
        ::nifti_1_header header = nifti_convert_nim2nhdr(image);
        // Print the header information
        status = disp_nifti_1_header(NULL, &header);
    }
    catch (...)
    {
        status = 2;
    }
    
    return status;
}
