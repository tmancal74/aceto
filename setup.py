# -*- coding: utf-8 -*-
from setuptools import find_packages
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
ext1 = Extension(name="aceto.nr3td_fi",
                 sources=["lib/nr3td_fi.f95"],
                 define_macros = [('F2PY_REPORT_ON_ARRAY_COPY','1')],
                 extra_f90_compile_args=["-I./lib"],
                 extra_f77_compile_args=["-I./lib"],
                 extra_link_args=["-L./lib -laceto"],
                )

setup(name = "aceto",
      version="0.0.1",
      
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

        # Specify the Python versions you support here. In particular, ensure        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
      ],

      keywords='physics, chemistry, quantum mechanics, open quantum systems',
      
      packages = find_packages(exclude=['lib','conf','src','tests','docs']),
      data_files = [("/usr/local/lib/",["./lib/libaceto.so"])],
      include_package_data = True, 
      ext_modules = [ext1],
                    

      
      )
