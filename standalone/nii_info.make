XFR = ln -sf

CPPFLAGS += -DNDEBUG -DHAVE_ZLIB
CFLAGS += -I. -Izlib
CXXFLAGS += -I.
LIBS += -lz -Lzlib

NIFTILIB2_CPPFLAGS = -DRNIFTI_NIFTILIB_VERSION=2

NIFTILIB1_OBJECTS = znzlib/znzlib.o niftilib/nifti1_io.o
NIFTILIB2_OBJECTS = znzlib/znzlib.o niftilib/nifti2_io.o

all: fetch zlib/libz.a nii_info nii2_info

fetch:
	mkdir -p zlib && cd zlib && $(XFR) ../../src/zlib/* .
	rm -f zlib/*.o
	mkdir -p znzlib
	cd znzlib && $(XFR) ../../src/znzlib/* ../../inst/include/znzlib/* .
	rm -f znzlib/*.o
	mkdir -p niftilib
	cd niftilib && $(XFR) ../../src/niftilib/* ../../inst/include/niftilib/* .
	rm -f niftilib/*.o
	mkdir -p RNifti
	cd RNifti && $(XFR) ../../inst/include/RNifti/* .
	$(XFR) ../inst/include/RNifti.h .

zlib/libz.a:
	cd zlib && ./configure --static && $(MAKE) libz.a CC="$(CC)"

nii_info: nii_info.cpp $(NIFTILIB1_OBJECTS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

nii2_info: nii_info.cpp $(NIFTILIB2_OBJECTS)
	$(CXX) $(CPPFLAGS) $(NIFTILIB2_CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -rf $(NIFTILIB1_OBJECTS) $(NIFTILIB2_OBJECTS) RNifti.h RNifti znzlib niftilib zlib nii_info nii2_info
