# -*- coding: utf-8 -*-
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction
from quantarhei import eigenbasis_of

from quantarhei import PopulationPropagator

from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
from quantarhei.spectroscopy.twod2 import TwoDSpectrumCalculator

import numpy
#import scipy
import time

from aceto.lab_settings import lab_settings
from aceto.band_system import band_system

import aceto.nr3td as nr3td

import matplotlib.pyplot as plt

t_start = time.time()

# Time axis for t1 and t3 times
Nr = 500
ta = TimeAxis(0.0, Nr, 2.0)

###############################################################################
#
# Define problem
#
###############################################################################

#
# define molecules
#
with energy_units("1/cm"):
    mol1 = Molecule(elenergies=[0.0, 12000.0])
    mol2 = Molecule(elenergies=[0.0, 12200.0])
    mol3 = Molecule(elenergies=[0.0, 12400.0])
    mol4 = Molecule(elenergies=[0.0, 12500.0])
mol1.position = [0.0, 0.0, 0.0]
mol2.position = [0.0, 10.0, 0.0]
mol3.position = [10.0, 0.0, 0.0]
mol4.position = [0.0, 0.0, 10.0]
mol1.set_dipole(0,1,[10.0, 0.0, 0.0])
mol2.set_dipole(0,1,[-10.0, 0.0, 0.0])
mol3.set_dipole(0,1,[0.0, 0.0, 10.0])
v1 = numpy.array([0.3, 0.3, 0.3])
v1 = 10.0*v1/numpy.sqrt(numpy.dot(v1,v1))
mol4.set_dipole(0,1,v1)

# rwa frequency as an average transition frequency
rwa = (mol1.elenergies[1]+mol2.elenergies[1])/2.0

lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
a_0 = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
a_90 = numpy.array([0.0, 1.0, 0.0], dtype=numpy.float64)
a_m45 = numpy.array([-1.0/numpy.sqrt(2.0), 
                     1.0/numpy.sqrt(2.0), 0.0], dtype=numpy.float64)
a_p45 = numpy.array([1.0/numpy.sqrt(2.0), 
                     1.0/numpy.sqrt(2.0), 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(a_0,a_90,a_m45,a_p45)


#
# System-bath interaction
#      

tsbi = TimeAxis(0.0, 3*Nr, 2.0)

params = dict(ftype="OverdampedBrownian", T=300, reorg=50.0, cortime=100.0)
with energy_units('1/cm'):
    cf = CorrelationFunction(tsbi, params)
    
mol1.set_transition_environment((0,1),cf)
mol2.set_transition_environment((0,1),cf)
mol3.set_transition_environment((0,1),cf)
mol4.set_transition_environment((0,1),cf)

#
# Creating aggregate
#      
agg = Aggregate("Dimer", molecules=[mol1, mol2, mol3, mol4])
#agg = Aggregate("Dimer", molecules=[mol1, mol2, mol3])

with energy_units("1/cm"):
    agg.set_resonance_coupling(0,1, 60.0)
    agg.set_resonance_coupling(1,2, 60.0)
    agg.set_resonance_coupling(0,2, 30.0)
    agg.set_resonance_coupling(1,3, 30.0)
    pass

#agg.set_coupling_by_dipole_dipole()    

agg.build(mult=2)

with energy_units("1/cm"):
    rwa_cm = agg.get_RWA_suggestion()
rwa = agg.get_RWA_suggestion()


#
# Calculate 2D spectra
#

# TimeAxis for t2 waiting time
t2s = TimeAxis(0.0, 5, 200.0)
tcalc = TwoDSpectrumCalculator(t1axis=ta, t2axis=t2s, t3axis=ta,
                               system=agg)
twods = tcalc.calculate(rwa, verbose=True, lab=lab)

#
# Show 2D spectra (save them)
#    

w1_min = rwa_cm - 700.0
w1_max = rwa_cm + 700.0
w3_min = rwa_cm - 700.0
w3_max = rwa_cm + 700.0

window_2D = [w1_min, w1_max, w3_min, w3_max]
print(window_2D)

with energy_units("1/cm"):

    
    # find info for normalization to the largest value across the whole set
#    k = 0
#    mx = numpy.zeros(len(t2s.data))
#    for tt2 in t2s.data:
#        spectrum = twods[k]
#        spect2D = numpy.real(spectrum.reph2D) + numpy.real(spectrum.nonr2D)
#        mx[k] = numpy.amax(spect2D)
#        k += 1
#        
#    mloc = numpy.argmax(mx)
#    mval = mx[mloc]
#    
#    print(mloc, mval)
#    mval = 1.0
#    
#    normalize = True
#
#    # normalize
#    k = 0
#    for tt2 in t2s.data:
#        
#        if normalize:
#            twods[k].devide_by(mval)
#        k += 1
#        
#    cbmax = twods[mloc].plot(axis=[w1_min, w1_max, w3_min, w3_max],vmax=1.0)
    
    k = 0
    for tt2 in t2s.data:
            
        #
        # Plotting with given units on axes
        #
        twods[k].plot(axis=[w1_min, w1_max, w3_min, w3_max]) #,vmax=1.0, cbmax=cbmax)

        figname = "fig"+str(round(tt2))+".png"
        print("saving file: ", figname, " with 2D spectrum at ", tt2, "fs")
        plt.savefig(figname)
        #plt.show()
        
        k += 1

t_end = time.time()
print("Finished at ", t_end-t_start, " secs") 




