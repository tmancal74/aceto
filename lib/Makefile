include ../conf/conf.in
#include ../conf/${COMPILER_SETTINGS}
include ../conf/pgi.in


COMP = ${FC} ${FFLAGS}
LINK = ${FC} ${LFLAGS}

#----------------------------------------------------------
# Library routines
#----------------------------------------------------------
trp2.o: trp2.f03
	${COMP} trp2.o  trp2.f03
	
	
	
#
# Predefined tasks
#

clean: 
	rm -rf *.o *.mod

delete: clean
	rm -rf libaceto*
	
