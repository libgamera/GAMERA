from distutils.core import setup, Extension
import string
import subprocess
GSL_LIBS = string.split(subprocess.check_output(['gsl-config','--libs']))
GSL_CFLAGS = string.split(subprocess.check_output(['gsl-config','--cflags']))
#GSL_LIBS = pkgconfig.libs('gsl')
print "->",GSL_LIBS
if(GSL_CFLAGS):
  COMPILEARGS=GSL_CFLAGS
  COMPILEARGS.append('-std=c++11')
else:
  COMPILEARGS=['-std=c++11']

extension_mod = Extension("_gappa",
                          ["_gappa.cc", "../src/Radiation.C","../src/Particles.C","../src/Utils.C","../src/Astro.C"],
                          extra_compile_args=COMPILEARGS,
                          #libraries=['gsl','gslcblas','m'],
                          #libraries=GSL_LIBS,
                          extra_link_args=GSL_LIBS,
                          include_dirs=['../include'],
)

setup(name = "gappa", ext_modules=[extension_mod])
