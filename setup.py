# -*- coding: utf-8 -*-
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

ext1 = Extension(name="aceto",
                 sources=["lib/trp2.f95"])

setup(name = "aceto",
      description = "Accelerated Charge and Energy Transfer Objects",
      ext_modules = [ext1]
      )



