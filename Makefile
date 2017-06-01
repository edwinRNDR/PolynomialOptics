
all: bin/ex1-postprocess bin/ex0-basicarithmetic

HEADERS=\
./TruncPoly/TruncPolySystem.hh \
./OpticalElements/OpticalMaterial.hh \
./OpticalElements/Cylindrical5.hh \
./OpticalElements/Propagation5.hh \
./OpticalElements/Spherical5.hh \
./OpticalElements/FindFocus.hh \
./OpticalElements/PointToPupil5.hh \
./OpticalElements/TwoPlane5.hh

LDFLAGS=-lm
LDFLAGS+=${shell pkg-config OpenEXR --libs}
CXXFLAGS=-fPIC -D_REENTRANT -D_THREAD_SAFE -D_GNU_SOURCE -Dcimg_use_openexr
CXXFLAGS+=${shell pkg-config OpenEXR --cflags}
CXXFLAGS+=-I. -ITruncPoly -IOpticalElements -Iinclude -g -Wall -fno-strict-aliasing

# Define this if not debugging:
OPTFLAGS=-O3 #-ffast-math -mfpmath=sse -march=native -fexpensive-optimizations -DNDEBUG

bin/ex1-postprocess: ${HEADERS} Example_PostprocessImage.cpp
	mkdir -p bin
	mkdir -p OutputPFM
	g++ ${LDFLAGS} ${CXXFLAGS} ${OPTFLAGS} Example_PostprocessImage.cpp -o bin/ex1-postprocess

bin/ex0-basicarithmetic: ${HEADERS} Example_BasicArithmetic.cpp
	mkdir -p bin
	g++ ${LDFLAGS} ${CXXFLAGS} ${OPTFLAGS} Example_BasicArithmetic.cpp -o bin/ex0-basicarithmetic

clean:
	rm -rf bin
