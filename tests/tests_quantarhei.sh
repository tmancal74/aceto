#! /bin/sh

# path to shared libraries for different platforms
lib1=LD_LIBRARY_PATH=`pwd`/../lib/:$LD_LIBRARY_PATH
lib2=DYLD_LIBRARY_PATH=`pwd`/../lib/:$DYLD_LIBRARY_PATH

# exporting variables
export ${lib1}
#export ${lib2}

# reporting
echo ${lib1}
echo ${lib2}

exec python quantarhei_test.py
 


