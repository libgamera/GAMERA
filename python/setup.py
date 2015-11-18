from distutils.core import setup, Extension
import os
import sys

# Use the bundeled pkgconfig
here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, here)
import pkgconfig

extra_link_args = pkgconfig.libs('gsl').split()

extra_compile_args = pkgconfig.cflags('gsl').split()
extra_compile_args.append('-std=c++11')

extension_mod = Extension(
    "_gappa",
    ["_gappa.cc", "../src/Radiation.C", "../src/Particles.C", "../src/Utils.C", "../src/Astro.C"],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    include_dirs=['../include'],
)

setup(
    name="gappa",
    ext_modules=[extension_mod],
)
