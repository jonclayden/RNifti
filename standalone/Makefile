CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I../inst/include
CXXFLAGS += -I../inst/include

RM = rm -f

NIFTILIB_OBJECTS = ../src/niftilib/nifti1_io.o ../src/znzlib/znzlib.o
ZLIB_OBJECTS = ../src/zlib/adler32.o ../src/zlib/compress.o ../src/zlib/crc32.o ../src/zlib/deflate.o ../src/zlib/gzclose.o ../src/zlib/gzlib.o ../src/zlib/gzread.o ../src/zlib/gzwrite.o ../src/zlib/infback.o ../src/zlib/inffast.o ../src/zlib/inflate.o ../src/zlib/inftrees.o ../src/zlib/trees.o ../src/zlib/uncompr.o ../src/zlib/zutil.o

all: nii_info

nii_info: nii_info.cpp $(NIFTILIB_OBJECTS) $(ZLIB_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(NIFTILIB_OBJECTS) $(ZLIB_OBJECTS) nii_info
