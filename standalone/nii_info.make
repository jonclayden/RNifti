CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I. -Izlib
CXXFLAGS += -I.
LIBS += -lz -Lzlib

NIFTILIB2_CPPFLAGS = -DRNIFTI_NIFTILIB_VERSION=2

NIFTILIB1_OBJECTS = znzlib/znzlib.o niftilib/nifti1_io.o
NIFTILIB2_OBJECTS = znzlib/znzlib.o niftilib/nifti2_io.o

all: zlib/libz.a nii_info nii2_info

zlib/libz.a:
	rm -rf zlib
	mkdir zlib && cd zlib && ln -s ../../src/zlib/* .
	cd zlib && ./configure --static && $(MAKE) libz.a CC="$(CC)"

niftilib/nifti1_io.c:
	mkdir -p niftilib
	cd niftilib && ln -s ../../src/niftilib/* ../../inst/include/niftilib/* .
	rm -f niftilib/*.o

niftilib/nifti2_io.c:
	mkdir -p niftilib
	cd niftilib && ln -s ../../src/niftilib/* ../../inst/include/niftilib/* .
	rm -f niftilib/*.o

znzlib/znzlib.c:
	mkdir -p znzlib
	cd znzlib && ln -s ../../src/znzlib/* ../../inst/include/znzlib/* .
	rm -f znzlib/*.o

nii_info: nii_info.cpp $(NIFTILIB1_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

nii2_info: nii_info.cpp $(NIFTILIB2_OBJECTS)
	$(CXX) $(CPPFLAGS) $(NIFTILIB2_CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -rf $(NIFTILIB1_OBJECTS) $(NIFTILIB2_OBJECTS) znzlib niftilib zlib nii_info nii2_info
