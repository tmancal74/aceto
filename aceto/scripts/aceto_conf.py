# -*- coding: utf-8 -*-

from shutil import copyfile
import os, sys
from pathlib import Path
import pkg_resources
from sys import platform as _platform


def main():
    
    print("Post installation tasks: copying binary shared libraries")
    #shlib_dir = os.path.join(os.path.dirname(aceto.__file__), 'shared_lib')

    # library filename
    version = pkg_resources.require("aceto")[0].version
    libname = "aceto-"+version+"-"+_platform
    libfile = "lib"+libname+".so"    
    
    if sys.version_info[0] < 3:
        raise Exception("This package needs Python 3!")
        
    if sys.version_info[1] > 4:
        # home
        home = str(Path.home())
    else:
        from os.path import expanduser
        home = expanduser("~")
        
    
    # working directory
    realPath = os.path.realpath(__file__) 
    here = os.path.dirname(realPath)
    
    # shared lib current location
    shlib = os.path.join(here,"..", "..","lib",libfile)
    
    # target location
    dest = os.path.join(home,"lib",libfile)
    
    ddir = os.path.dirname(dest)      
    
    go_copy = False
    while True:
        yn = input("Copy shared librarty into "+dest+"? (y/n)")
        
        if yn == "y" or yn == "Y" :
            ddir = os.path.dirname(dest)
            if os.path.exists(ddir):
                if os.path.isdir(ddir):
                    go_copy = True
                else:
                    raise Exception(dest+" is not a directory")
                    
            else:
                
                yn = input("Directory "+ddir+" does not exist. Create? (y/n)")
                if yn == "n" or yn == "N":
                    
                    yn = input("Would you like to specify a different destinatiob directory? (y/n)")
                    if yn == "y" or yn == "Y":
                        ddir = input("Shared library destination: ")
                        dest = os.path.join(ddir,"libaceto.so")
                        go_copy = False
                    else:
                        print("Installation not finished ... aborting!")
                        break
                    
                elif yn == "y" or yn == "Y":
                    os.mkdir(ddir)
                    go_copy = True
               
            if go_copy:
                print("from ", shlib, " to ", dest)
                copyfile(shlib,dest)
                print("Post installation tasks finished")
                break
            
        if yn == "n" or yn == "N":
            yn = input("Would you like to specify a different destinatiob directory? (y/n)")
            if yn == "y" or yn == "Y":
                ddir = input("Shared library destination: ")
                dest = os.path.join(ddir,libfile)
                go_copy = False
            else:
                print("Installation not finished ... aborting!")
                break
    
    
    
if __name__ == "__main__":
    
    main()
    
    