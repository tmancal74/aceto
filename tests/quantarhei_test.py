# -*- coding: utf-8 -*-
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction
from quantarhei import eigenbasis_of

from quantarhei import PopulationPropagator

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
    mol4 = Molecule(elenergies=[0.0, 12800.0])
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

#
# System-bath interaction
#      

tsbi = TimeAxis(0.0, 3*Nr, 2.0)

params = dict(ftype="OverdampedBrownian", T=300, reorg=200.0, cortime=100.0)
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


###############################################################################
#
# Create band_system from quantarhei classes
#
###############################################################################

#
# hamiltonian and transition dipole moment operators
#
H = agg.get_Hamiltonian()
D = agg.get_TransitionDipoleMoment()

#
# Construct band_system object
#
Nb = 3
Ns = numpy.zeros(Nb, dtype=numpy.int)
Ns[0] = 1
Ns[1] = agg.nmono
Ns[2] = Ns[1]*(Ns[1]-1)/2
sys = band_system(Nb, Ns)


#
# Set energies
#
en = numpy.zeros(sys.Ne, dtype=numpy.float64)
#if True:
with eigenbasis_of(H):
    for i in range(sys.Ne):
        en[i] = H.data[i,i]
    sys.set_energies(en)

    #
    # Set transition dipole moments
    #
    dge_wr = D.data[0:Ns[0],Ns[0]:Ns[0]+Ns[1],:]
    def_wr = D.data[Ns[0]:Ns[0]+Ns[1],(Ns[0]+Ns[1]):(Ns[0]+Ns[1]+Ns[2]),:]

    dge = numpy.zeros((3,Ns[0],Ns[1]), dtype=numpy.float64)
    deff = numpy.zeros((3,Ns[1],Ns[2]), dtype=numpy.float64)
    
    for i in range(3):
        dge[i,:,:] = dge_wr[:,:,i]
        deff[i,:,:] = def_wr[:,:,i]
    sys.set_dipoles(0,1,dge)
    sys.set_dipoles(1,2,deff)


#
# Relaxation rates
#
KK = agg.get_RedfieldRateMatrix()

# relaxation rate in single exciton band
Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]]*10.0
#print(1.0/Kr)

sys.init_dephasing_rates()
sys.set_relaxation_rates(1,Kr)


#
# Lineshape functions
#
sbi = agg.get_SystemBathInteraction()
cfm = sbi.CC
cfm.create_double_integral()


#
# Transformation matrices
#
SS = H.diagonalize()
SS1 = SS[1:Ns[1]+1,1:Ns[1]+1]
SS2 = SS[Ns[1]+1:,Ns[1]+1:]
H.undiagonalize()

sys.set_gofts(cfm._gofts)    # line shape functions
sys.set_sitep(cfm.cpointer)  # pointer to sites
sys.set_transcoef(1,SS1)      # matrix of transformation coefficients  
sys.set_transcoef(2,SS2)      # matrix of transformation coefficients  



#
# define lab settings
#
lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(X, X, X, X)

#
# Other parameters
#

dt = ta.step
t1s = ta.data 
t3s = ta.data 

# TimeAxis for t2 waiting time
t2s = TimeAxis(0.0, 5, 100.0)

#
# Finding population evolution matrix
#
prop = PopulationPropagator(ta, Kr)
Uee, Uc0 = prop.get_PropagationMatrix(t2s, corrections=True)

tc = 0
teetoos = t2s.data
for tt2 in teetoos:

    #
    # Initialize response storage
    #
    resp_r = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')
    resp_n = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')

    (it2, err) = ta.locate(tt2) 
    print("t2 = ", tt2, "fs (it2 = ", it2,")")
    
    #
    # calcute response
    #
    print("calculating response: ")
    rmin = 0.0001
    t1 = time.time()
    
    print(" - ground state bleach")
    # GSB
    nr3td.nr3_r3g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r) 
    nr3td.nr3_r4g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)

    print(" - stimulated emission")
    # SE
    nr3td.nr3_r1g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    nr3td.nr3_r2g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    
    print(" - excited state absorption")
    # ESA
    nr3td.nr3_r1fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    nr3td.nr3_r2fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    
    # Transfer
    sys.set_population_propagation_matrix(Uee[:,:,tc]-Uc0[:,:,tc])
    
    print(" - stimulated emission with transfer")    
    # SE
    nr3td.nr3_r1g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    nr3td.nr3_r2g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)

    print(" - excited state absorption with transfer") 
    # ESA
    nr3td.nr3_r1fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    nr3td.nr3_r2fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    
    
    t2 = time.time()
    print("... calculated in ",t2-t1, " sec")
    
    
    #
    # Calculate corresponding 2D spectrum
    #
    
    ftresp = numpy.fft.fft(resp_r,axis=1)
    ftresp = numpy.fft.ifft(ftresp,axis=0)
    reph2D = numpy.fft.fftshift(ftresp)
    
    ftresp = numpy.fft.ifft(resp_n,axis=1)
    ftresp = numpy.fft.ifft(ftresp,axis=0)*ftresp.shape[1]
    nonr2D = numpy.fft.fftshift(ftresp)
    
    atype = ta.atype
    ta.atype = 'complete'
    oa1 = ta.get_FrequencyAxis() 
    oa1.data += rwa
    oa1.start += rwa
    ta.atype = atype
    
    atype = ta.atype
    ta.atype = 'complete'
    oa3 = ta.get_FrequencyAxis() 
    oa3.data += rwa
    oa3.start += rwa
    ta.atype = atype
    
        

    #
    # Show 2D specrtrum
    #

    tot2D_real = numpy.real(reph2D) + numpy.real(nonr2D)
    
    if tc == 0:
        max_tot = numpy.max(tot2D_real)
            
    w1_min = 11000.0
    w1_max = 13000.0
    w3_min = 11000.0
    w3_max = 13000.0
    
    with energy_units("1/cm"):
        #
        # Plotting with given units on axes
        #
        #print("Frequency span w1: ", oa1.min, oa1.max)      
        (i1_min, dist) = oa1.locate(w1_min)
        (i1_max, dist) = oa1.locate(w1_max)
        #print("Frequency span w3: ", oa3.min, oa3.max)      
        (i3_min, dist) = oa3.locate(w3_min)
        (i3_max, dist) = oa3.locate(w3_max)    
    
        realout = numpy.real(tot2D_real[i1_min:i1_max,i3_min:i3_max])/max_tot
    
        Ncontour = 100
        plt.contourf(oa1.data[i1_min:i1_max],oa3.data[i3_min:i3_max],
                     realout, Ncontour)
    
    
    
    
    
    t_end = time.time()
    print("Finished at ", t_end-t_start, " secs")
    figname = "fig"+str(round(tt2))+".png"
    print(figname)
    plt.savefig(figname)
    
    #plt.show()

    tc += 1
