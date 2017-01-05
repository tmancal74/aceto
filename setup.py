# -*- coding: utf-8 -*-
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

ext1 = Extension(name="aceto",
                 sources=["lib/trp2.f95"],
                 define_macros = [('F2PY_REPORT_ON_ARRAY_COPY','1')],
                 # this is the flag gfortran needs to process OpenMP directives
                 extra_compile_args = ['-fopenmp '],
                 extra_link_args = ['-fopenmp '])

setup(name = "aceto",
      description = "Accelerated Charge and Energy Transfer Objects",
      ext_modules = [ext1]
      )

