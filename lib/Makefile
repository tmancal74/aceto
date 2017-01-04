include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}


COMP = ${FC} ${PICFLAGS} ${FFLAGS} 
LINK = ${CC} ${DYNFLAGS}

# List of library routines
LIBO = trp2.o

libaceto.so: echo_building_library ${LIBO}
	${LINK} libaceto.so ${LIBO} ${DYNLIBS}
	@echo "...library built"

echo_building_library:
	@echo "Building library ..."
	
#----------------------------------------------------------
# Library routines
#----------------------------------------------------------
trp2.o: trp2.f95
	${COMP} trp2.o  trp2.f95
	
	
	
#
# Predefined tasks
#

.PHONY: clean delete

clean: 
	rm -rf *.o *.mod

delete: clean
	rm -rf libaceto*
	
