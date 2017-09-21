# -*- coding: utf-8 -*-
#import os, sys
import setuptools  # needed for interoperation with numpy.distutils
from setuptools import find_packages
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from sys import platform as _platform

# To use a consistent encoding
from codecs import open
from os import path

# FIXME: delete in the next version
# --- NO MORE NEEDED --- now the post installation task is handled by a
# configuration script
#from distutils.command.install import install as _install
#
#
#def _post_install(dir):
#    from subprocess import call
#    call([sys.executable, os.path.join('scripts','postinstall.py')],
#         cwd=os.path.join(dir, 'aceto'))
#
#
#class install(_install):
#    def run(self):
#        _install.run(self)
#        self.execute(_post_install, (self.install_lib,),
#                     msg="Running post install task")


# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
 
# library name is comprised of the version number and platform 
version = "0.0.4"

libname = "aceto-"+version+"-"+_platform
libfile = "lib"+libname+".so"
liblink = "-l"+libname  
liblocn = path.join(".","lib",libfile)  
Ilib = "-I"+path.join(".","lib")
Llib = "-L"+path.join(".","lib")
    

ext1 = Extension(name="aceto.nr3td_fic",
                 sources=[path.join("lib","nr3td_fic.f90")],
                 define_macros = [('F2PY_REPORT_ON_ARRAY_COPY','1')],
                 extra_f90_compile_args=[Ilib],
                 extra_f77_compile_args=[Ilib],                                  
                 extra_link_args=[Llib+" "+liblink],
                )



setup(name = "aceto",
      version=version,
      
      description = "Accelerated Charge and Energy Transfer Objects",
      long_description=long_description, 

      # The project's main homepage.
      url='https://github.com/tmancal74/aceto',

      # Author details
      author='Tomas Mancal',
      author_email='mancal@karlov.mff.cuni.cz',

      # Choose your license
      license='MIT',
      

      # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure        
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
      ],

      keywords='physics, chemistry, quantum mechanics, open quantum systems',
      
      packages = find_packages(exclude=['lib','conf','src','tests','docs']),
      
      data_files = [("lib/",[liblocn])],
      include_package_data = True, 
      
      ext_modules = [ext1],
       
      # FIXME: delete in next version             
      # NO MORE NEEDED
      # This will be use when I figure out how to copy a file from within an egg
      #
      #cmdclass={'install': install},

        
      entry_points={
        'console_scripts': ['aceto_conf=aceto.scripts.aceto_conf:main'],
      }


)

