# GAMERA

is a project launched by the [Max Planck Insitute for Nuclear Physics in Heidelberg (MPIK)](https://www.mpi-hd.mpg.de/mpi/en/),
an open-source C++/python package which handles the spectral modelling of non-thermally emitting astrophysical sources in a simple and modular way. It allows the user to devise time-dependent models of leptonic and hadronic particle populations in a general astrophysical context (including SNRs, PWNs and AGNs) and to compute their subsequent photon emission. 

For more info and a turorial, see the [GAMERA docu!](http://libgamera.github.io/GAMERA/docs/main_page.html)
(currently being updated)

GAMERA is written in C++ and can be wrapped to python.

The software is listed in the Astrophysics Source Code Library <a href="https://ascl.net/2203.007"><img src="https://img.shields.io/badge/ascl-2203.007-blue.svg?colorB=262255" alt="ascl:2203.007" /></a>

## Quick Start

### Dependencies

You need to have [GSL](http://www.gnu.org/software/gsl/) installed on your system.

### Building the CPP library

This project is using `cmake` as build system.
Do the usual configure / build / install steps:

```
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~/.local/gamera
$ cmake --build build
$ cmake --install build
```

### Building the python bindings

This project is using [`scikit-build-core`](https://scikit-build-core.readthedocs.io/en/latest/).
To install the python bindings to GAMERA (gappa), just run:

```
$ pip install git+https://github.com/libgamera/GAMERA
```

For development, you can use an editable mode installation that will
automatically rebuild the project on import using.
First create a virtual environment, install the build dependencies into it
and then install gappa in editable mode:

```
$ python -m venv venv
$ source venv/bin/activate
$ pip install cmake swig scikit-build-core
$ pip install -e . --no-build-isolation --config-settings=editable.rebuild=true --config-settings=build-dir=build
```

## Licence


The GAMERA library and the GAMERA programs are free software;
you can redistribute them and/or modify it under the terms of
the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License,
or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

A copy of the GNU Lesser General Public License version 2.1 can be found
[here](https://github.com/libgamera/GAMERA/blob/master/licenses/lgpl-2.1.txt).




