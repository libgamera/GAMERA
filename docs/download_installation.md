[(back to main page)](main_page.md)

Download
========

GAMERA is available on [github](https://github.com/libgamera/GAMERA).
There, you can download a [zip archive](https://github.com/libgamera/GAMERA/archive/master.zip) with the source code or clone it with git via 
```
$ git clone https://github.com/libgamera/GAMERA.git
```

Installation
============

Dependencies
------------
You will have to compile `GAMERA` from source. Apart from a `C++` compiler (e.g. using 
`gcc` on Linux or `clang` on Mac), you need to install the following dependencies: 
- the GNU scientific library `gsl`
- the `C/C++` python wrapper `swig`, if you want the python-wrapped version (i.e. `GAPPA`)
- python (including the headers) for `GAPPA`

On `Linux`:

- On Ubuntu 16.04 LTS
  
  These packages are readily available on all major Linux platforms. On Ubuntu 16.04 LTS,  
  you can install them via 
  ```
  sudo apt install libgsl2 libgsl-dev swig python-dev
  ```
- On Ubuntu 18.04 LTS
  ```
  sudo apt install libgsl23 libgsl-dev swig python-dev pkg-config
  ```

On `MacOS`: 

many repositories provide `gsl` and `swig`, for example [homebrew](https://brew.sh/): 
```
$ brew install gsl swig
```
on some MacOS systems, issues are encountered in the compilation phase. These can be fixed by adding 2 lines to the file \<GAMERA-path\>/python/setup.py

after line 20 add:
```
extra_compile_args.append(“-stdlib=libc++”)
extra_link_args.append(“-stdlib=libc++”)
```
on some other MacOS versions (e.g. Mojave 10.14.5) a further issue can stop the compilation of the gappa library: the swig version that can be installed does not support anymore the option `-nosafecstrings` option.
The workaround is to remove it from the Makefile (at line 70).

Building: `python`
------------------
__Compilation__

If you are interested in the `python` package (`GAPPA`), run 
```
$ make gappa
```
and it will generate `lib/_gappa.so` and `lib/gappa.py`. 

__Usage__

In your python script, You need to add the directory holding the module (default: `lib/` in the `GAMERA` directory) to the search  path for modules, e.g.:
```
import sys
sys.path.append('/home/user/Documents/GAMERA/lib')
```

Building: `C++`
---------------
__Compilation__

To build the `C++` library, go to the base GAMERA directory and run
```
$ make gamera
```
This will create the shared object `lib/libgamera.so` which can be used in `C++` programs. 

__Usage__

To use the library in you program, compile like e.g.:
```
$ gcc -o yourProgram yourProgramSource.C -lgamera -L$(LIBDIR)
```
where `LIBDIR` is the directory where `libgamera.so` is located (default: `lib/` in the `GAMERA` directory). 
Also, it might be necessary to export the lib directory to `LD_LIBRARY_PATH`, e.g.: 
```
$ export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH
```


Updating `GAMERA`
-----------------
After any new `git pull`, you should clean up the directory via 
```
$ make clean
```
and re-compile with the commands above.

[(back to main page)](main_page.md)
