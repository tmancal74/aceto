#! /bin/sh

# path to shared libraries for different platforms
#lib1=LD_LIBRARY_PATH=${HOME}/lib/:$LD_LIBRARY_PATH
#lib2=DYLD_LIBRARY_PATH=${HOME}/../lib/:$DYLD_LIBRARY_PATH

# exporting variables
#export ${lib1}
#export ${lib2}

# reporting
#echo ${lib1}
#echo ${lib2}

#exec python aceto_2D_test.py
 
echo "Testing integration of ACETO with Quantarhei"
echo "Running without failure means passing the test"
echo "Short calculation will be started below by 'qrhei' script"
echo " "
exec qrhei aceto_2D_test.py



