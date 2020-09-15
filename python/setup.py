from distutils.core import setup, Extension
import distutils.sysconfig
import os
import sys


# Use the bundeled pkgconfig
here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, here)
import pkgconfig

cfg_vars = distutils.sysconfig.get_config_vars()

for key,value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes","")

extra_link_args = pkgconfig.libs('gsl').split()

extra_compile_args = pkgconfig.cflags('gsl').split()
extra_compile_args.append('-m64')
extra_compile_args.append('-std=c++11')

extension_mod = Extension(
    "_gappa",
    ["_gappa.cc",
     "../src/Radiation.C",
     "../src/Particles.C",
     "../src/Utils.C",
     "../src/Astro.C",
     "../src/2D_interp/interp2d.C",
     "../src/2D_interp/interp2d_spline.C",
     "../src/2D_interp/bilinear.C",
     "../src/2D_interp/bicubic.C"],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    include_dirs=['../include'],
)

setup(
    name="gappa",
    ext_modules=[extension_mod],
)
