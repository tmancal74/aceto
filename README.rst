
Accelerated Charge and Energy Transfer Objects (ACETO) 
======================================================
 
Aceto is a Fortran (2003) library containing routines for some common tasks in quantum theory of molecular charge
and energy transfer. 

Preferred pronounciation of the library's abbreviation is the Italian one. ACETO means vinegar
in Italian. The library is designed to add some extra flavour to another project, the Python
open quantum system theory package Quantarhei (see http://github.com/tmancal74/quantarhei).


HOW TO INSTALL ACETO
--------------------

Aceto can be installed on Linux and Mac

Installation from source:
-------------------------

The safest way of installing Aceto is by downloading source code from GitHub. You need 'git'
and a Fortran compiler installed on your system. Currently we support 'gnufortran', but
some more compiler will follow soon. Get the source code by typing:

> git clone https://github.com/tmancal74/aceto

This will create a directory 'aceto' in your working directory. To compile,
enter the directory:

> cd aceto

Aceto uses autotools to manage compilation of the Fortran source code. Compilation proceeds in the usual autotools
way, i.e.

> ./configure
> make

This results in compilation of the Fortran/C part of the aceto library. To create an installable python package
type the following:

> python setup.py bdist_wheel 

which compiles Python part of the package and links it with the Fortran part. It also creates a *.whl file in the 'dist/' subdirectory
To install the wheel, enter the 'dist' directory and use 'pip'

> cd dist
> pip install 'aceto_wheel_file'

Use the file you find in the 'dist' directory instead of 'aceto_wheel_file'.
I recommend to enter the 'dist' directory before you install - in this way
the installation is the same as if you just downloaded the wheel file from
internet. Some problems with 'pip uninstall' can occur if you install from 
the source directory.


Binary installation:
--------------------

The installation procedure of Aceto is still in development. Currently, without
any warranty, we provide binary distribution for macOS. Windows will follow
soon (it is hoped). Use 'pip' to install the binary version:

> pip install aceto

You might need to have gfortran installed on your system even if you install
binaries due to dependency on the fortran shared libraries.


