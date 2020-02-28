CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I.
CXXFLAGS += -I.

RM = rm -f

NIFTILIB_OBJECTS = niftilib/nifti1_io.o niftilib/nifti2_io.o znzlib/znzlib.o
ZLIB_OBJECTS = zlib/adler32.o zlib/compress.o zlib/crc32.o zlib/deflate.o zlib/gzclose.o zlib/gzlib.o zlib/gzread.o zlib/gzwrite.o zlib/infback.o zlib/inffast.o zlib/inflate.o zlib/inftrees.o zlib/trees.o zlib/uncompr.o zlib/zutil.o

all: nii_info

nii_info: nii_info.cpp $(NIFTILIB_OBJECTS) $(ZLIB_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(NIFTILIB_OBJECTS) $(ZLIB_OBJECTS) nii_info
