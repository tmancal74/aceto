# -*- coding: utf-8 -*-

from shutil import copyfile
import os
from pathlib import Path
import aceto
from os import path

def main():
    
    print("Post installation tasks: copying binary shared libraries")
    #shlib_dir = os.path.join(os.path.dirname(aceto.__file__), 'shared_lib')

    home = str(Path.home())
    realPath = os.path.realpath(__file__)  # /home/user/test/my_script.py
    here = os.path.dirname(realPath)
    shlib = os.path.join(here,"lib","libaceto.so")
    dest = os.path.join(home,"lib","libaceto.so")
    
    while True:
        yn = input("Copy shared librarty into "+dest+"? (y/n)")
        if yn == "y" or yn == "Y" :
            print("from ", shlib, " to ", dest)
            copyfile(shlib,dest)
            break
        if yn == "n" or yn == "N":
            print("Installation not finished ... aborting!")
            break
    
    
    
if __name__ == "__main__":
    
    main()
    
    