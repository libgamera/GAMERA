
OUTDIR = $(shell pwd)/out
LIBDIR = $(shell pwd)/lib


# define include path
INCLUDES := -I./include

# GNU scientific library stuff
GSLCFLAGS := $(shell gsl-config --cflags)

GSLLIBS := $(shell gsl-config --libs)

# define additional switches to be passed to the compiler
CXXFLAGS := $(GSLCFLAGS) -std=c++11

# define additional switches to be passed to the linker
LDFLAGS := $(GSLLIBS)


help:
	@echo ''
	@echo 'GAMERA available make targets:'
	@echo ''
	@echo '  gamera           Create the GAMERA shared object (lib/libgamera.so)'
	@echo '  gamerapy         Wrap the python module (lib/gamerapy.py and lib/_gamerapy.so)'
	@echo '  tutorial1        Compile a simple tutorial program'
	@echo '                   (-> bin/TutorialLevel1, see http://joachimhahn.github.io/GAMERA/)'
	@echo '  documentation    Generate documentation'
	@echo '                   (-> docu/doxygen/html/ and docu/sphinx/build/html/)'
	@echo '  clean            remove temporary files'


all: gamera TutorialLevel1 gamerapy documentation

gamera: Radiation Particles Utils libgamera

tutorial1: gamera TutorialLevel1

# create the individual .o files
Radiation : src/Radiation.C
	$(CXX) -g -O2 -fpic -Wall -c src/Radiation.C -o $(OUTDIR)/Radiation.o $(CXXFLAGS) $(INCLUDES)
Particles : src/Particles.C
	$(CXX) -g -O2 -fpic -Wall -c src/Particles.C -o $(OUTDIR)/Particles.o $(CXXFLAGS) $(INCLUDES)
Utils : src/Utils.C
	$(CXX) -g -O2 -fpic -Wall -c src/Utils.C -o $(OUTDIR)/Utils.o $(CXXFLAGS) $(INCLUDES)

# create the shared library
objectsSO = $(OUTDIR)/Radiation.o $(OUTDIR)/Particles.o $(OUTDIR)/Utils.o
libgamera : $(objectsSO)
	@echo "Generating library $@..."
	$(CXX) -shared -o $(LIBDIR)/libgamera.so $(objectsSO) $(CXXFLAGS) $(LDFLAGS)
	@echo "-> done!"

# This is an example on how to use the shared library libgamera.so
TutorialLevel1 : docu/tutorial/TutorialLevel1.C
	$(CXX) -g -O2 -Wall -c docu/tutorial/TutorialLevel1.C -o $(OUTDIR)/TutorialLevel1.o $(CXXFLAGS) $(INCLUDES)
	$(CXX) -g -Wall -o bin/TutorialLevel1 $(OUTDIR)/TutorialLevel1.o -lgamera -L$(LIBDIR)

# make gamerapy package
gamerapy:
	cd python;\
	swig -python -c++ -nosafecstrings $(GSLCFLAGS) -outdir ../lib -o _gamerapy.cc gamerapy.i;\
	python setup.py build_ext --build-lib ../lib;\
	cd ..;\

# make documentation
documentation:
	echo pwd;\
	cd docu/doxygen;\
	doxygen doxygenconfig.cfg;\
	cd ../sphinx/;\
	make html;\

clean-out: 
	-rm -f out/*
clean-bin: 
	-rm -f bin/*
clean-lib: 
	-rm -f lib/*
clean-python: 
	-rm -f python/_gamerapy.cc python/_gamerapy.so python/gamerapy.py
	-rm -rf python/build/
clean-docu:
	-rm -rf docu/doxygen/xml/* docu/doxygen/html/* docu/sphinx/build/*

clean: clean-out clean-bin clean-lib clean-python clean-docu
