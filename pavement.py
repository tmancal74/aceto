"""
    Aceto Package Managemant and Tests
    ==================================

    To compile aceto conveniently, the following packages have to be installed

    paver (pip install paver)

    The following tasks are available. Run them with `paver` as

    .. code:: bash

        $ paver TASK_NAME

    To test the whole suit, just run `paver` without any task.

    .. code:: bash

        $ paver


"""
import contextlib
import os
import subprocess
import platform


from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh


version = "0.0.9"

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
@task
def bootstrap():
    sh('./bootstrap')

@task
def configure():
    sh('./configure')

@needs('configure')
@task
def compile():
    sh('make')

###############################################################################
#  Creates distribution
###############################################################################
@needs('compile')
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
