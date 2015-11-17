from distutils.core import setup, Extension
import pkgconfig

GSL_CFLAGS = pkgconfig.cflags('gsl')
GSL_LIBS = pkgconfig.libs('gsl')

if(GSL_CFLAGS):
  COMPILEARGS=['-std=c++11', GSL_CFLAGS]
else:
  COMPILEARGS=['-std=c++11']

extension_mod = Extension("_gappa",
                          ["_gappa.cc", "../src/Radiation.C","../src/Particles.C","../src/Utils.C","../src/Astro.C"],
                          extra_compile_args=COMPILEARGS,
                          libraries=['gsl','gslcblas'],
                          include_dirs=['../include'],
)

setup(name = "gamerapy", ext_modules=[extension_mod])
