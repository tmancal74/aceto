
            Accelerated Charge and Energy Transfer Objects (ACETO) library
 
Contains Fortran (2003) routines for some common tasks in quantum theory of molecular charge
and energy transfer. 

Preferred pronounciation of the library's abbreviation is the Italian one. ACETO means vinegar
in Italian. The library is designed to add some extra flavour to another project, the Python
open quantum system theory package Quantarhei (see http://github.com/tmancal74/quantarhei).


HOW TO INSTALL ACETO
--------------------

Aceto can be installed on Linux and Mac

Installation from source:

The safest way of installing Aceto is by downloading source code from github.com. You need
to configur the Makefile by changing the content of the 'conf/conf.in' file to point to
a file containing gcc flags (gcc_linux.in and gcc_mac.in files are tested). Of course
your system has to have gcc and gfortran compilers installed. Then you
need to create a 'lib' directory in your home directory. This is a temporal 
fix, but right now, shared library 'libaceto-version-platform.so' is "installed" 
locally to this directory. Then you need to issue

> make
> make install

series of commnads. You will be asked to confirm that a shared library that was
created in the 'lib' subdirectory can be copied to the 'lib' directory in your 
home directory.

Binary installation:

The installation procedure of Aceto is still in development. Currently we
provide binary distribution for macOS and Linux compiled with gcc complilers
through a Python egg awailable from PyPa via the 'easy_install' command. Typing

> easy_install aceto

will instal aceto, but in order for it to run correctly, you have to type

> aceto_conf

after the installation is over.

You will be asked to confirm that the 'libaceto-version-platform.so' file can be
moved to the directory 'lib' in your home directory. If this directory is not present, you
will be asked to confirm its creation.


Linux specific:

On Linux it seems that LD_LIBRARY_PATH variable set in qrhei script which is
used to run the quantarhei input files (as a subprocess) does not influence the setting for
the subprocess. The solution is to export the LD_LIBRARY_PATH in something like
.bashrc to point to the ${HOME}/lib directory.