# -*- coding: utf-8 -*-
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction
from quantarhei import eigenbasis_of

from quantarhei import PopulationPropagator

import numpy
import scipy
import time

from aceto.lab_settings import lab_settings
from aceto.band_system import band_system

import aceto.nr3td as nr3td

import matplotlib.pyplot as plt

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

Nr = 500

#
# System-bath interaction
#      
ta = TimeAxis(0.0, Nr, 2.0)

params = dict(ftype="OverdampedBrownian", T=300, reorg=200.0, cortime=100.0)
with energy_units('1/cm'):
    cf = CorrelationFunction(ta, params)
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

#print(D.data.shape)
#
#for n in range(Ns[1]):
#    m = 0
#    for k in range(Ns[1]):
#        for l in range(k+1,Ns[1]):
#            print(n,m+1+Ns[1],"(",k,l,")", D.data[1+n, 1+Ns[1]+m,:])
#            m += 1

#raise Exception()

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

if False:
    AK = numpy.zeros((D.data.shape[0],D.data.shape[1]))
    A1 = numpy.zeros((1,Ns[1]))
    A2 = numpy.zeros((Ns[1],Ns[2]))
    
    for n in range(D.data.shape[0]):
        for m in range(D.data.shape[1]):
            v = D.data[n,m,:]
            v2 = numpy.dot(v,v)
            AK[n,m] = v2
    print(AK)
    for n in range(Ns[1]):
        v = dge[:,0,n]
        v2 = numpy.dot(v,v)  
        A1[0,n] = v2
    for n in range(Ns[1]):
        for m in range(Ns[2]):
            v = deff[:,n,m]
            v2 = numpy.dot(v,v)  
            A2[n,m] = v2        
    print(A1)
    print(A2)

    k = 0
    for n in range(Ns[1]):
        for m in range(n+1,Ns[1]):
            print("*: ", n,m)
            d1 = dge[:,0,m]
            d2 = deff[:,n,k]
            
            print("0", m, " (g->S) = ", d1)
            print(n,k, " (S->D) = ", d2)
            print(n," -> ",k, "(",n,",",m,") is ",m,": difference is ", d2-d1)
            print(AK[1+n,Ns[1]+k],AK[0,1+m])
            k += 1
            

    raise Exception()




#
# Relaxation rates
#
KK = agg.get_RedfieldRateMatrix()
Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]]*10.0
#Kr[0,0] = -0.0000001
#Kr[1,0] = 0.0000001
print(1.0/Kr)
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

#NN1 = SS1.shape[0]
#SS1 = numpy.eye(NN1)
#NN2 = SS2.shape[0]
#SS2 = numpy.eye(NN2)

sys.set_gofts(cfm._gofts)    # line shape functions
sys.set_sitep(cfm.cpointer)  # pointer to sites
sys.set_transcoef(1,SS1)      # matrix of transformation coefficients  
sys.set_transcoef(2,SS2)      # matrix of transformation coefficients  


#sys._check_twoex_dipoles()


#
# define lab settings
#
lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(X, X, X, X)

#
# Initialize response storage
#

#
# Other parameters
#

dt = ta.step
t1s = ta.data 
t3s = ta.data 

#t2 = 200.0

teetoos = [i*100.0 for i in range(5)]
t2s = TimeAxis(0.0, 5, 100.0)

prop = PopulationPropagator(ta, Kr)
Uee = prop.get_PropagationMatrix(t2s)

tc = 0
for tt2 in teetoos:

    resp_r = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')
    resp_n = numpy.zeros((Nr, Nr), dtype=numpy.complex128, order='F')

    (it2, err) = ta.locate(tt2) 
    print("it2 = ", it2)
    
    #
    # calcute response
    #
    print("calculating response: ")
    rmin = 0.0001
    t1 = time.time()
    
    # GSB
    nr3td.nr3_r3g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r) 
    nr3td.nr3_r4g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)

    # SE
    nr3td.nr3_r1g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    nr3td.nr3_r2g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    
    # ESA
    nr3td.nr3_r1fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    nr3td.nr3_r2fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    
    # Transfer
    sys.set_population_propagation_matrix(Uee[:,:,tc])
    
    # SE
    nr3td.nr3_r1g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    nr3td.nr3_r2g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)

    # ESA
    nr3td.nr3_r1fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
    nr3td.nr3_r2fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
    
    
    t2 = time.time()
    print("... calculated in ",t2-t1, " sec")
    
    
    #
    # Show corresponding 2D spectrum
    #
    
    ftresp = numpy.fft.fft(resp_r,axis=1)
    ftresp = numpy.fft.ifft(ftresp,axis=0)
    ftresp_r = numpy.fft.fftshift(ftresp)
    
    ftresp = numpy.fft.ifft(resp_n,axis=1)
    ftresp = numpy.fft.ifft(ftresp,axis=0)*ftresp.shape[1]
    ftresp_n = numpy.fft.fftshift(ftresp)
    
    om1 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t1s),
                                                            d=dt)) + rwa
    om3 = 2.0*scipy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(len(t3s),
                                                            d=dt)) + rwa
    
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
    plt.contourf(om1[Nst:Nfi],om3[Nst:Nfi],
                 numpy.real(tots[Nst:Nfi,Nst:Nfi]),100)
    
    t_end = time.time()
    print("Finished in ", t_end-t_start, " secs")
    figname = "fig"+str(round(tt2))+".png"
    print(figname)
    plt.savefig(figname)
    
    #plt.show()

    tc += 1
