
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

INCS	:=	$(wildcard include/*.h)
SRCS	:=	$(patsubst include%,src%,$(patsubst %.h,%.C,$(INCS)))
OBJS    :=	$(patsubst include%,$(OUTDIR)%,$(patsubst %.h,%.o,$(INCS)))
LIB		:=	$(LIBDIR)/libgamera.so
help:
	@echo ''
	@echo 'GAMERA available make targets:'
	@echo ''
	@echo '  gamera           Create the GAMERA shared object (lib/libgamera.so)'
	@echo '  gappa            Wrap the python module (lib/gappa.py and lib/_gappa.so)'
	@echo '  tutorial1        Compile a simple tutorial program'
	@echo '                   (-> bin/TutorialLevel1, see http://joachimhahn.github.io/GAMERA/)'
	@echo '  documentation    Generate documentation'
	@echo '                   (-> docu/doxygen/html/ and docu/sphinx/build/html/)'
	@echo '  clean            remove temporary files'


all: gamera TutorialLevel1 gappa

gamera: $(LIB)

tutorial1: gamera TutorialLevel1

$(OUTDIR)/%.o:	src/%.C include/%.h
				@echo "Compiling" $<
				$(CXX) -g -O2 -fpic -Wall -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

# create the shared library
$(LIB) : $(OBJS)
	@echo "Generating library $@..."
	$(CXX) -shared -o $@ $(OBJS) $(CXXFLAGS) $(LDFLAGS)
	@echo "-> done!"

# This is an example on how to use the shared library libgamera.so
TutorialLevel1 : docu/tutorial/TutorialLevel1.C
	$(CXX) -g -O2 -Wall -c docu/tutorial/TutorialLevel1.C -o $(OUTDIR)/TutorialLevel1.o $(CXXFLAGS) $(INCLUDES)
	$(CXX) -g -Wall -o bin/TutorialLevel1 $(OUTDIR)/TutorialLevel1.o -lgamera -L$(LIBDIR)

# make gappa package
gappa:
	cd python;\
	swig -python -c++ -nosafecstrings -outdir ../lib -o _gappa.cc gappa.i;\
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
	-rm -f python/_gappa.cc python/_gappa.so python/gappa.py
	-rm -rf python/build/
clean-docu:
	-rm -rf docu/doxygen/xml/* docu/doxygen/html/* docu/sphinx/build/*

clean: clean-out clean-bin clean-lib clean-python clean-docu
