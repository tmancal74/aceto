Aceto library code structure
============================

Fortran
-------

acetolib.f95
    classes and their code
        acetocnf.f95 aceto
        acetosys.f95 band_system
        acetolab.f95 lab_settings
                
    types and constants
        acetodef.f95
    
    auxiliary routines
        acetoaux.f95

    abbreviated routines
        e.g. nr3td.f95 

    code of the full interface routines
        e.g. nr3td_code.f95
        
    
acetoaux.f95
    auxiliary routines
    

aceto_fi.f95
    interfaces to "full interface" routines
        e.g. nr3td_fi.f95
        

Libraries
---------
 
    libaceto.so (to be used from fortran)
        acetocfn.f95
        acetosys.f95
        acetolab.f95
        acetodef.f95
        acetoaux.f95
        aceto_fi.f95
        nr3td.f95 etc.
        nr3td_fi.f95 etc.
        nr3td_code.f95 etc.
                
        
How to use it from Fortran
--------------------------

To use aceto from Fortran in standard way, you should "use" the module¨
"acetolib", and instantiate an object of the type "aceto" with required 
properties 


    use acetolib

    type(aceto) :: ac_cuda
    integer     :: err
    
    call ac_cuda%init(ACETO_GPU_CUDA, err)


The object is used to call the library routines with required settings


    call ac_cuda%nr3td_r2g(sys, lab, t2, t1s, t2s, rwa, resp)
    

The line above calls a CUDA accelerated version of the routine "nr3td_r2g".    

It is also possible to directly access the "full interface" routines by 
using the module "aceto_fi"


    use aceto_fi
    
    call nr3td_r2g_cuda(... a lot of arguments ...)
        
        
    
Python
------

aceto modules
    classes and their code
        aceto_conf.py 
        acetosys.py band_system
        acetolab.py lab_settings
    
    types and constants
        acetodef.py
        
    auxiliary routines
        acetodef.py
        
    abbreviated routines
        e.g. nr3td.py
        
    
Python uses the standard ctype package to load the "libaceto.so" file 
and calls the "full interface" routines directly.
 