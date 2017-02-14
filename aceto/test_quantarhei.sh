#! /bin/sh

export LD_LIBRARY_PATH=`pwd`/../lib/:`pwd`/../src/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=`pwd`/../lib/:`pwd`/../src/:$DYLD_LIBRARY_PATH
echo LD_LIBRARY_PATH=`pwd`/../lib/:`pwd`/../src/:$LD_LIBRARY_PATH
exec python quantarhei_test.py



