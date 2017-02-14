# -*- coding: utf-8 -*-
import numpy

from lab_settings import lab_settings
from band_system import band_system

import nr3td

#
# Define system
#
Ns = [1,2,1]
sys = band_system(3, Ns)

en = numpy.array([0.0, 0.9, 1.1, 2.0], dtype=numpy.float64)
sys.set_energies(en)

d_ge = numpy.zeros((3,1,2), dtype=numpy.float64, order='F')
d_ef = numpy.zeros((3,2,1), dtype=numpy.float64, order='F')
d_ge[:,0,1] = numpy.array([0.0, 2.0, 0.0],dtype=numpy.float64)
d_ef[:,1,0] = numpy.array([0.0, 3.0, 0.0],dtype=numpy.float64)
sys.set_dipoles(0,1,d_ge)
sys.set_dipoles(1,2,d_ef)

Kr = numpy.zeros((2,2), dtype=numpy.float64)
Kr[0,0] = -1.0/500.0
Kr[1,0] = 1.0/500.0
Kr[0,1] = 1.0/50.0
Kr[1,1] = -1.0/50.0
sys.init_dephasing_rates()
sys.set_relaxation_rates(1,Kr)

#
# define lab settings
#
lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(X, X, X, X)

#
# Initialize response storage
#
resp = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')

#
# Other parameters
#
t2 = 0.0
dt = 1.0
Nr = 1000
t1s = numpy.array([k*dt for k in range(Nr)],dtype=numpy.float64)
t3s = numpy.array([k*dt for k in range(Nr)],dtype=numpy.float64)

#
# calcute response
#
nr3td.nr3_r2g(lab, sys, t2, t1s, t3s, resp)

print(resp)


