CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I. -Izlib
CXXFLAGS += -I.
LIBS += -lz -Lzlib

NIFTILIB2_CPPFLAGS = -DRNIFTI_NIFTILIB_VERSION=2 -DNO_REMAP_NIFTI2_FUNCTIONS

NIFTILIB1_OBJECTS = niftilib/nifti1_io.o znzlib/znzlib.o
NIFTILIB2_OBJECTS = niftilib/nifti2_io.o znzlib/znzlib.o

all: zlib/libz.a nii_info nii2_info

zlib/libz.a:
	rm -rf zlib
	mkdir zlib && cd zlib && ln -s ../../src/zlib/* .
	cd zlib && ./configure --static && $(MAKE) libz.a CC="$(CC)"

nii_info: nii_info.cpp $(NIFTILIB1_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

nii2_info: nii_info.cpp $(NIFTILIB2_OBJECTS)
	$(CXX) $(CPPFLAGS) $(NIFTILIB2_CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -rf $(NIFTILIB1_OBJECTS) $(NIFTILIB2_OBJECTS) zlib nii_info nii2_info
