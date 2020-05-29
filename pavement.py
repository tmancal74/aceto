"""
    Quantarhei Package Tests
    ========================
    
    To run the tests, the following packages have to be installed
    
    paver (pip install paver)
    nose (pip install nose)
    aloe (pip install aloe)
    behave (pip install behave)
    pylint (pip install pylint)
    
    About our tests
    ---------------
    
    Quantarhei is tested against unit tests run by the `nose` tool, against
    acceptance tests run by `aloe` tool and against doc tests (run by `nose`).
    We also calculate coverege using the `coverage` tool and we plan to 
    check the code by `pylint`.
    
    
    The ultimit goal of testing in 100 % coverage and perfect code style.
    
    
    The following tasks are available. Run them with `paver` as
    
    .. code:: bash
    
        $ paver TASK_NAME
        
    To test the whole suit, just run `paver` without any task.

    .. code:: bash

        $ paver
        
        
        
    Doctest Tasks
    -------------
    
    doc_tests
    doc_tests_v
    doc_tests_vs
    doc_tests_cov
    doc_tests_cov_v
    doc_tests_cov_vs
        
    Runs test doc tests. _cov means `--with-coverage` option, _v means
    with `-v` option for verbose output and _vs means with `-vs` option, i.e.
    verbose and printing standard output of the tests.
        
        
    Unit tests
    ----------
    
    unit_tests_vs
    unit_tests_cov_vs


    Acceptance tests
    ----------------
    
    aloe_tests_vs
    aloe_tests_cov_vs
    
    
    Tests to be run during development
    ----------------------------------
    
    doc_dev
    unit_dev
    
    Edit the list of files or directories that you want to test during
    development.
    

"""
import contextlib
import os
import subprocess
import platform


from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh


version = "0.0.8"

sys_name = platform.system()



pip = 'pip'
python = 'python'

# 
# Commands for deleting files and directories silently and without error codes
#
if sys_name == "Darwin" or sys_name == "Linux":
    deldir = 'rm -r -f '
    delfile = 'rm -r -f '
elif sys_name == "Windows":
    deldir = 'rmdir /s /q '
    delfile = 'del /s /q '
else:
    raise Exception("Unknown system")

#
# The Mother of all repositories
#
repository = 'https://github.com/tmancal74/aceto'

    

#
# Context manager for getting into subdirectories
#
@contextlib.contextmanager
def cd(path):
   old_path = os.getcwd()
   os.chdir(path)
   try:
       yield
   finally:
       os.chdir(old_path)
       
def rm_rf(path):
    """Removal of files and directories, recursively and silently
    
    """
    pass
    
#
###############################################################################
#
#      Standard developer tasks
#
###############################################################################
#      
###############################################################################
#  Creates distribution
###############################################################################   
@task
def wheel():
    sh(python+" setup.py bdist_wheel")

###############################################################################
# Installs Quantarhei from the source
###############################################################################
@needs('wheel')
@task
def inst():
    
    sh(pip+' install dist/a*.whl')
    
#
#  The same as above
#
@needs('inst')
@task
def install():
    pass

###############################################################################
#  Uninstall local installation of Quantarhei
###############################################################################
@needs('uninst')
@task
def uninstall():
    pass

#
# The same as above  
#
@task
def uninst():
	sh(pip+' uninstall -y aceto')
    

###############################################################################
#  Upload to pypi 
###############################################################################
@needs('wheel')
@task
def upload():
	pass
    #sh('twine upload dist/quantarhei-'+version+'.tar.gz')


###############################################################################
#  Clean-up 
###############################################################################
@task
def clean():
    try:
        sh(deldir+'dist')
    except:
        print("Directory not present - no problem")
    try:    
        sh(delfile+'aceto.egg-info')
    except:
        print("File not present - no problem")
    try:
        sh(deldir+'result_images')
    except:
        print("Directory not present - no problem")
    try:
        sh("make clean")
    except:
        print("Make clean failed")


###############################################################################
# Reinstallation with clean-up 
###############################################################################
@needs('clean','uninst', 'inst')
@task
def reinst():
    pass

###############################################################################
# Local tests: this will reinstall Quantarhei and run tests 
###############################################################################
@needs('reinst')
@task
def local_tests():
	sh('paver')


###############################################################################
#
###############################################################################
@task
def tasks():
	print("")
	print("Quantarhei Paver Tasks")
	print("======================")
	print("")
	print("Essential tasks: ")
	print("----------------")
	print("inst        ... install quantarhei from this source code")
	print("reinst      ... uninstall and install from this source code")
	print("wheel       ... create distribution")
	print("clean       ... clean the repository")
	print("")


