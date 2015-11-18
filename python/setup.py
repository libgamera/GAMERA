from distutils.core import setup, Extension
import string
import sys
import subprocess
GSL_LIBS = subprocess.check_output(['gsl-config','--libs']).split()
GSL_CFLAGS = subprocess.check_output(['gsl-config','--cflags']).split()

if(GSL_CFLAGS):
  COMPILEARGS=GSL_CFLAGS
  COMPILEARGS.append('-std=c++11')
else:
  COMPILEARGS=['-std=c++11']

if sys.platform == 'darwin':
  COMPILEARGS.append('-Wno-error=shorten-64-to-32')

extension_mod = Extension("_gappa",
                          ["_gappa.cc", "../src/Radiation.C","../src/Particles.C","../src/Utils.C","../src/Astro.C"],
                          extra_compile_args=COMPILEARGS,
                          extra_link_args=GSL_LIBS,
                          include_dirs=['../include'],
)

setup(name = "gappa", ext_modules=[extension_mod])
