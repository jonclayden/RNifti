#include "lib/print.h"
#include "niftilib/nifti1_io.h"
#include "lib/NiftiImage.h"

int main (int argc, char* argv[])
{
    if (argc < 2)
        return 1;
    
    int status = 0;
    
    try
    {
        const std::string path(argv[1]);
        RNifti::NiftiImage image(path, false);
        ::nifti_1_header header = nifti_convert_nim2nhdr(image);
        status = disp_nifti_1_header(NULL, &header);
    }
    catch (...)
    {
        status = 2;
    }
    
    return status;
}
