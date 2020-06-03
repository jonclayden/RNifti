CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I.
CXXFLAGS += -I.

NIFTILIB2_CPPFLAGS = -DRNIFTI_NIFTILIB_VERSION=2 -DNO_REMAP_NIFTI2_FUNCTIONS

RM = rm -f

NIFTILIB1_OBJECTS = niftilib/nifti1_io.o znzlib/znzlib.o
NIFTILIB2_OBJECTS = niftilib/nifti2_io.o znzlib/znzlib.o
ZLIB_OBJECTS = zlib/adler32.o zlib/compress.o zlib/crc32.o zlib/deflate.o zlib/gzclose.o zlib/gzlib.o zlib/gzread.o zlib/gzwrite.o zlib/infback.o zlib/inffast.o zlib/inflate.o zlib/inftrees.o zlib/trees.o zlib/uncompr.o zlib/zutil.o

all: nii_info nii2_info

nii_info: nii_info.cpp $(NIFTILIB1_OBJECTS) $(ZLIB_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

nii2_info: nii_info.cpp $(NIFTILIB2_OBJECTS) $(ZLIB_OBJECTS)
	$(CXX) $(CPPFLAGS) $(NIFTILIB2_CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(NIFTILIB1_OBJECTS) $(NIFTILIB2_OBJECTS) $(ZLIB_OBJECTS) nii_info nii2_info
