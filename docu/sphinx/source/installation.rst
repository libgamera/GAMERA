Installation and Usage
======================


Dependencies
------------

``GAMERA`` has the following dependencies:

* `GSL <http://www.gnu.org/software/gsl/>`_
* `SWIG <http://www.swig.org/>`_ (to wrap the python package)
* `Sphinx <http://sphinx-doc.org/>`_, `Breathe <https://breathe.readthedocs.org/>`_ and `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ (optional, to build the docu)

**on Ubuntu**

using apt::

    $ sudo apt-get install libgsl0ldbl libgsl0-dev python-dev swig

(optional)::

    $ sudo apt-get install doxygen python-pip
    $ sudo pip install sphinx breathe

**on MacOS**
using `homebrew <http://brew.sh/>`_ ::

    $ brew install gsl swig

(optional)::

    $ brew install doxygen python-pip
    $ pip install sphinx breathe

Downloading ``GAMERA``
----------------------

via git::

    $ git clone https://github.com/JoachimHahn/GAMERA.git

Manual download from `github <https://github.com/JoachimHahn/GAMERA>`_

 .. note::

   Installation of ``gappa``, the python-wrapped version of ``GAMERA``, via
   ``pypi`` will be available soon.

Building ``GAMERA``
-------------------
To build the ``C++`` library, issue ::

    $ make gamera

which will create ``lib/libgamera.so``.

If you want to build ``gappa``, run::

    $ make gappa

and you will get ``lib/_gappa.so`` and ``lib/gappa.py``.

 .. warning::

   On MacOS, building ``GAMERA`` is only supported with ``clang``


Usage in ``C++``
----------------

To use the library, compile like e.g::

    $ gcc -o yourProgram yourProgramSource.C -lgamera -L$(LIBDIR)

where ``LIBDIR`` is the directory where ``libgamera.so`` is located (default: ``./lib``).
Also, it might be necessary to export the ``lib`` directory to ``LD_LIBRARY_PATH``, e.g.::

    $ export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH

Usage in ``python``
-------------------

You need to add the directory holding the module (default: ``./lib``) to the search path for modules in your script, e.g.:

.. sourcecode:: python

   import sys
   import os
   sys.path.append(os.path.abspath('./lib'))

.. warning::

  Currently, ``GAMERA`` doesn't work with ``python3``, use ``python2`` instead!
