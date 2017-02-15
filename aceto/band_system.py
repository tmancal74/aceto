# -*- coding: utf-8 -*-
import numpy

class band_system:
    """Python implementation of acetosys band_system class


    """    
    def __init__(self, Nb, Ns):
        self.Nb = Nb
        self.Ns = numpy.array(Ns, dtype=numpy.int)
        self.Ne = numpy.sum(Ns)
        
        self.en = None
        self.om01 = None
        self.om12 = None
        self.nn01 = None
        self.nn12 = None
        self.dd01 = None
        self.dd12 = None
        self.Kr11 = None
        self.Kr22 = None
        self.Kd01 = None
        self.Kd11 = None
        self.Kd12 = None
        self.gofts = None
        self.ptn = None
        self.SS1 = None
        self.SS2 = None
        
        
    def set_energies(self, en):
        self.en = en
        self.om01 = numpy.zeros((self.Ns[0],self.Ns[1]), 
                                dtype=numpy.float64, order='F')
        self.om12 = numpy.zeros((self.Ns[1],self.Ns[2]), 
                                dtype=numpy.float64, order='F')
        for i in range(self.Ns[0]):
            for j in range(self.Ns[1]):
                self.om01[i,j] = self.en[i] - self.en[self.Ns[0]+j]
        for i in range(self.Ns[1]):
            for j in range(self.Ns[2]):
                self.om12[i,j] = (self.en[self.Ns[0]+i] 
                                - self.en[self.Ns[0]+self.Ns[1]+j])
                
    def set_dipoles(self, Nbi, Nbf, dab):
        if not (dab.shape == (3,self.Ns[Nbi],self.Ns[Nbf])):
            raise Exception("Wrong shape of dipole matrix")
        if (Nbi == 0) and (Nbf == 1):
            self.nn01 = dab
            self.dd01 = numpy.zeros((self.Ns[0],self.Ns[1]),
                                    dtype=numpy.float64, order='F')
            N1 = self.Ns[0]
            N2 = self.Ns[1]
            for i in range(N1):
                for j in range(N2):
                    dd = numpy.sqrt(numpy.dot(self.nn01[:,i,j],
                                              self.nn01[:,i,j]))
                    self.dd01[i,j] = dd
                    if dd != 0.0:
                        self.nn01[:,i,j] = self.nn01[:,i,j]/dd
        elif (Nbi == 1) and (Nbf == 2):
            self.nn12 = dab
            self.dd12 = numpy.zeros((self.Ns[1],self.Ns[2]),
                                    dtype=numpy.float64, order='F')
            N1 = self.Ns[1]
            N2 = self.Ns[2]
            for i in range(N1):
                for j in range(N2):
                    dd = numpy.sqrt(numpy.dot(self.nn12[:,i,j],
                                              self.nn12[:,i,j]))
                    self.dd12[i,j] = dd
                    if dd != 0.0:
                        self.nn12[:,i,j] = self.nn12[:,i,j]/dd
        else:
            raise Exception("Attempt to assing unsupported dipole block")
            

    def set_gofts(self,gofts):
        self.gofts = gofts
        
    def set_sitep(self, ptn):
        self.ptn = ptn
        self.fptn = self.ptn + 1
        
    def set_transcoef(self, Nb, SS):
        if Nb == 1:
            self.SS1 = SS
        elif Nb == 2:
            self.SS2 = SS
        else:
            raise Exception("Attempt to assign unsupported block")

    
    def set_relaxation_rates(self, Nb, RR):
        if Nb == 1:
            self.Kr11 = RR
        elif Nb == 2:
            self.Kr22 = RR
        else:
            raise Exception("Attempt to set usupported rate block")
        self.update_dephasing_rates(Nb)
    
    def update_dephasing_rates(self, Nb):
        if Nb == 1:
            for i in range(self.Ns[0]):
                for j in range(self.Ns[1]):
                    self.Kd01[i,j] += self.Kr11[j,j]/2.0
            for i in range(self.Ns[1]):
                for j in range(self.Ns[1]):
                    self.Kd11[i,j] += (self.Kr11[i,i] + self.Kr11[j,j])/2.0
            for i in range(self.Ns[1]):
                for j in range(self.Ns[2]):
                    self.Kd12[i,j] += self.Kr11[i,i]/2.0
        elif Nb == 2:
            for i in range(self.Ns[1]):
                for j in range(self.Ns[2]):
                    self.Kd12[i,j] += self.Kr22[j,j]/2.0            
        else:
            raise Exception("Attempt to update unsupported "+
                            "dephasing rate block")
            
    
    def init_dephasing_rates(self):
        self.Kd01 = numpy.zeros((self.Ns[0], self.Ns[1]), 
                                    dtype=numpy.float64, order='F')
        self.Kd11 = numpy.zeros((self.Ns[1], self.Ns[1]), 
                                    dtype=numpy.float64, order='F')
        self.Kd12 = numpy.zeros((self.Ns[1], self.Ns[2]), 
                                    dtype=numpy.float64, order='F')
    
    
        
        
