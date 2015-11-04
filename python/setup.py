from distutils.core import setup, Extension
import pkgconfig

GSL_CFLAGS = pkgconfig.cflags('gsl')
GSL_LIBS = pkgconfig.libs('gsl')

extension_mod = Extension("_gamerapy", 
                          ["_gamerapy.cc", "../src/Radiation.C","../src/Particles.C","../src/Utils.C"],
                          extra_compile_args=['-std=c++11', GSL_CFLAGS],
                        #   extra_link_args=[GSL_LIBS],
                          libraries=['gsl','gslcblas'],
                          include_dirs=['../include'],
)

setup(name = "gamerapy", ext_modules=[extension_mod])
