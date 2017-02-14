# -*- coding: utf-8 -*-
"""
    Time-dependent third order non-linear response of a three band
    multi-level system


"""

import nr3td_fi

def nr3_r1g(lab, sys, t2, t1s, t3s, resp):
    """ Calculates R2g response function
    
    """
    pass    


def nr3_r2g(lab, sys, t2, t1s, t3s, rwa, resp):
    """ Calculates R2g response function
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    t2 : float
        The so-called waiting or population time of non-linear spectroscopic
        techniques
        
    t1s : float array
        Values of t1 time for which response should be calculated
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    resp : complex 2d array
        Non-linear response 
        
    """
#
#    For debugging, check if all arrays are fortran continuous
#
#    print(lab.orient_aver.flags['F_CONTIGUOUS'])
#    print(sys.Ns.flags['F_CONTIGUOUS'])
#    print(sys.om01.flags['F_CONTIGUOUS'])
#    print(sys.nn01.flags['F_CONTIGUOUS'])
#    print(sys.dd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd11.flags['F_CONTIGUOUS'])
#    print(t1s.flags['F_CONTIGUOUS'])
#    print(t3s.flags['F_CONTIGUOUS'])
#    print(resp.flags['F_CONTIGUOUS'])
    
    
    nr3td_fi.nr3_r2g_fi(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, t2, t1s, t3s, rwa, resp)


def nr3_r3g(lab, sys, t2, t1s, t3s, resp):
    """ Calculates R2g response function
    
    
    """
    pass

def nr3_r4g(lab, sys, t2, t1s, t3s, resp):
    """ Calculates R2g response function
    
    
    """
    pass

def nr3_r1f(lab, sys, t2, t1s, t3s, resp):
    """ Calculates R2g response function
    
    """
    pass

def nr3_r2f(lab, sys, t2, t1s, t3s, resp):
    """ Calculates R2g response function
    
    """
    pass

