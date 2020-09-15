
OUTDIR = $(shell pwd)/out
LIBDIR = $(shell pwd)/lib


# define include path
INCLUDES := -I./include

# GNU scientific library stuff
GSLCFLAGS := $(shell gsl-config --cflags)

GSLLIBS := $(shell gsl-config --libs)

# define additional switches to be passed to the compiler
CXXFLAGS := $(GSLCFLAGS) -m64 -std=c++11

# define additional switches to be passed to the linker
LDFLAGS := $(GSLLIBS)


help:
	@echo ''
	@echo 'GAMERA available make targets:'
	@echo ''
	@echo '  gamera           Create the GAMERA shared object (lib/libgamera.so)'
	@echo '  gappa            Wrap the python module (lib/gappa.py and lib/_gappa.so)'
	@echo '  clean            remove temporary files'


all: gamera gappa

gamera: Radiation Particles Utils Astro bicubic bilinear interp2d interp2d_spline libgamera




# create the individual .o files
Radiation : src/Radiation.C
	$(CXX) -g -O2 -fpic -Wall -c src/Radiation.C -o $(OUTDIR)/Radiation.o $(CXXFLAGS) $(INCLUDES)
Particles : src/Particles.C
	$(CXX) -g -O2 -fpic -Wall -c src/Particles.C -o $(OUTDIR)/Particles.o $(CXXFLAGS) $(INCLUDES)
Utils : src/Utils.C
	$(CXX) -g -O2 -fpic -Wall -c src/Utils.C -o $(OUTDIR)/Utils.o $(CXXFLAGS) $(INCLUDES)
Astro : src/Astro.C
	$(CXX) -g -O2 -fpic -Wall -c src/Astro.C -o $(OUTDIR)/Astro.o $(CXXFLAGS) $(INCLUDES)



bicubic : src/2D_interp/bicubic.C
	$(CXX) -g -O2 -fpic -Wall -c src/2D_interp/bicubic.C -o $(OUTDIR)/bicubic.o $(INCLUDES)
bilinear : src/2D_interp/bilinear.C
	$(CXX) -g -O2 -fpic -Wall -c src/2D_interp/bilinear.C -o $(OUTDIR)/bilinear.o $(INCLUDES)
interp2d : src/2D_interp/interp2d.C
	$(CXX) -g -O2 -fpic -Wall -c src/2D_interp/interp2d.C -o $(OUTDIR)/interp2d.o $(INCLUDES)
interp2d_spline : src/2D_interp/interp2d_spline.C
	$(CXX) -g -O2 -fpic -Wall -c src/2D_interp/interp2d_spline.C -o $(OUTDIR)/interp2d_spline.o $(INCLUDES)


# create the shared library
objectsSO = $(OUTDIR)/Radiation.o $(OUTDIR)/Particles.o $(OUTDIR)/Utils.o $(OUTDIR)/Astro.o $(OUTDIR)/bicubic.o $(OUTDIR)/bilinear.o $(OUTDIR)/interp2d.o $(OUTDIR)/interp2d_spline.o
libgamera : $(objectsSO)
	@echo "Generating library $@..."
	$(CXX) -shared -o $(LIBDIR)/libgamera.so $(objectsSO) $(CXXFLAGS) $(LDFLAGS)
	@echo "-> done!"


# make gappa package
gappa:
	cd python && \
	swig -python -c++ -outdir ../lib -o _gappa.cc gappa.i && \
	python setup.py build_ext --build-lib ../lib


clean-out:
	-rm -f out/*
clean-bin:
	-rm -f bin/*
clean-lib:
	-rm -f lib/*
clean-python:
	-rm -f python/_gappa.cc python/_gappa.so python/gappa.py
	-rm -rf python/build/

clean: clean-out clean-bin clean-lib clean-python
