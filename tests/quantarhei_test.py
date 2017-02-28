# -*- coding: utf-8 -*-
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction
from quantarhei import eigenbasis_of

import numpy
import scipy
import time

from aceto.lab_settings import lab_settings
from aceto.band_system import band_system

import aceto.nr3td as nr3td

t_start = time.time()

###############################################################################
#
# Create band_system from quantarhei classes
#
###############################################################################

#
# define molecules
#
with energy_units("1/cm"):
    mol1 = Molecule(elenergies=[0.0, 12000.0])
    mol2 = Molecule(elenergies=[0.0, 12400.0])
#    mol3 = Molecule(elenergies=[0.0, 12200.0])
#    mol4 = Molecule(elenergies=[0.0, 12200.0])
mol1.position = [0.0, 0.0, 0.0]
mol2.position = [0.0, 10.0, 0.0]
#mol3.position = [10.0, 0.0, 0.0]
#mol4.position = [0.0, 0.0, 10.0]
mol1.set_dipole(0,1,[2.0, 0.0, 0.0])
mol2.set_dipole(0,1,[1.0, 0.0, 0.0])
#mol3.set_dipole(0,1,[0.0, 0.0, 1.0])
#mol4.set_dipole(0,1,[0.3, 0.3, 0.3])

# rwa frequency as an average transition frequency
rwa = (mol1.elenergies[1]+mol2.elenergies[1])/2.0

#
# System-bath interaction
#      
ta = TimeAxis(0.0, 1000, 1.0)
params = dict(ftype="OverdampedBrownian", T=300, reorg=20.0, cortime=100.0)
with energy_units('1/cm'):
    cf = CorrelationFunction(ta, params)
mol1.set_transition_environment((0,1),cf)
mol2.set_transition_environment((0,1),cf)
#mol3.set_transition_environment((0,1),cf)
#mol4.set_transition_environment((0,1),cf)

#
# Creating aggregate
#      
#agg = Aggregate("Dimer", molecules=[mol1, mol2, mol3, mol4])
agg = Aggregate("Dimer", molecules=[mol1, mol2])

with energy_units("1/cm"):
    agg.set_resonance_coupling(0,1, 300.0)
#    agg.set_resonance_coupling(1,2, 300.0)
#    agg.set_resonance_coupling(0,2, 600.0)
#    agg.set_resonance_coupling(1,3, 600.0)
agg.build(mult=2)

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

    #raise Exception()
#    for i in range(Ns[1]):
#        for j in range(Ns[2]):
#            print(i,j, sys.dd12[i,j])
#
# Relaxation rates
#
KK = agg.get_RedfieldRateMatrix()
Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]]
#print(Kr)
sys.init_dephasing_rates()
sys.set_relaxation_rates(1,Kr)


#
# Lineshape functions
#
sbi = agg.get_SystemBathInteraction()
cfm = sbi.CC
cfm.create_double_integral()
#print(cfm.cpointer)
#print(cfm._gofts.shape)

#
# Transformation matrices
#
SS = H.diagonalize()
SS1 = SS[1:Ns[1]+1,1:Ns[1]+1]
SS2 = SS[Ns[1]+1:,Ns[1]+1:]
H.undiagonalize()
#print(SS)

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
# Initialize response storage
#
Nr = 1000
resp_r = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')
resp_n = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')

#
# Other parameters
#
t2 = 100.0
dt = ta.step
t1s = ta.data 
t3s = ta.data 

(it2, err) = ta.locate(t2) 
print("it2 = ", it2)

#
# calcute response
#
print("calculating response: ")
rmin = 0.01
t1 = time.time()

nr3td.nr3_r2g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
nr3td.nr3_r3g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)

nr3td.nr3_r1g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
nr3td.nr3_r4g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)

nr3td.nr3_r1fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
nr3td.nr3_r2fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)


t2 = time.time()
print("... calculated in ",t2-t1, " sec")

#
# Show corresponding 2D spectrum
#
import matplotlib.pyplot as plt

ftresp = numpy.fft.fft(resp_r,axis=0)
ftresp = numpy.fft.ifft(ftresp,axis=1)
ftresp_r = numpy.fft.fftshift(ftresp)

ftresp = numpy.fft.ifft(resp_n,axis=0)
ftresp = numpy.fft.ifft(ftresp,axis=1)*ftresp.shape[1]
ftresp_n = numpy.fft.fftshift(ftresp)

om1 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t1s), d=dt)) + rwa
om3 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t3s), d=dt)) + rwa

print("max = ", numpy.max(numpy.real(ftresp_r)))
print("max = ", numpy.max(numpy.real(ftresp_n)))
                                     
tots = numpy.real(ftresp_r) + numpy.real(ftresp_n)
#tots = numpy.real(ftresp_n)


#plt.contourf(om1,om3,numpy.real(ftresp_r),100)
#plt.contourf(om1,om3,numpy.real(ftresp_n),100)
Np = resp_r.shape[0]
Nst = 44*Np//100 
print(Nst, Np)
Nfi = Np - Nst
plt.contourf(om1[Nst:Nfi],om3[Nst:Nfi],numpy.real(tots[Nst:Nfi,Nst:Nfi]),100)

t_end = time.time()
print("Finished in ", t_end-t_start, " secs")

#plt.savefig("fig"+str(it2)+".jpg")

plt.show()

