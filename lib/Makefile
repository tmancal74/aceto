include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}
#include ../conf/pgi.in


COMP = ${FC} ${PICFLAGS} ${FFLAGS} 
LINK = ${CC} ${DYNFLAGS}

# List of library routines
LIBO = trp2.o

libaceto.so: ${LIBO}
	${LINK} libaceto.so ${LIBO} ${DYNLIBS}

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
	
