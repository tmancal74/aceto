# -*- coding: utf-8 -*-
import numpy
import scipy
import time

from lab_settings import lab_settings
from band_system import band_system

import nr3td

#
# Define system
#
Ns = [1,2,1]
sys = band_system(3, Ns)

en = numpy.array([0.0, 0.92, 1.08, 2.0], dtype=numpy.float64)
rwa = 1.0
sys.set_energies(en)

d_ge = numpy.zeros((3,1,2), dtype=numpy.float64, order='F')
d_ef = numpy.zeros((3,2,1), dtype=numpy.float64, order='F')
d_ge[:,0,1] = numpy.array([0.0, 1.0, 0.0],dtype=numpy.float64)
d_ge[:,0,0] = numpy.array([1.0/numpy.sqrt(2.0), 1.0/numpy.sqrt(2.0), 0.0],
                            dtype=numpy.float64)
d_ef[:,1,0] = numpy.array([0.0, 3.0, 0.0],dtype=numpy.float64)
d_ef[:,0,0] = numpy.array([0.0, 3.0, 0.0],dtype=numpy.float64)
sys.set_dipoles(0,1,d_ge)
sys.set_dipoles(1,2,d_ef)

Kr = numpy.zeros((2,2), dtype=numpy.float64)
Kr[0,0] = -1.0/50.0
Kr[1,0] = 1.0/50.0
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
Nr = 1000
resp = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')

#
# Other parameters
#
t2 = 300.0
dt = 1.0
t1s = numpy.array([k*dt for k in range(Nr)],dtype=numpy.float64)
t3s = numpy.array([k*dt for k in range(Nr)],dtype=numpy.float64)

#
# calcute response
#
print("calculating response: ")
t1 = time.time()
nr3td.nr3_r2g(lab, sys, t2, t1s, t3s, rwa, resp)
t2 = time.time()
print("... calculated in ",t2-t1, " sec")
#
# Show corresponding 2D
#
import matplotlib.pyplot as plt
ftresp = numpy.fft.fft(resp,axis=0)
ftresp = numpy.fft.ifft(ftresp,axis=1)
ftresp = numpy.fft.fftshift(ftresp)
om1 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t1s), d=dt)) + rwa
om3 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t3s), d=dt)) + rwa

plt.contourf(om1,om3,numpy.real(ftresp),100)
plt.show()




