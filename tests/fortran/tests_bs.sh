#! /bin/sh

export LD_LIBRARY_PATH=`pwd`/../../lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=`pwd`/../../lib/:$DYLD_LIBRARY_PATH
echo LD_LIBRARY_PATH=`pwd`/../../lib/:$LD_LIBRARY_PATH
exec ./band_system_test.x


