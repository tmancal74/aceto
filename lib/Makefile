include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}


COMP = ${FC} -I../src  ${PICFLAGS} ${FFLAGS}  
LINK = ${CC} ${DYNFLAGS}

# List of library routines
LIBO = trp2.o nresp3o.o nr3td.o 

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

nresp3o.o: nresp3o.f95
	${COMP} nresp3o.o nresp3o.f95


nr3td_fi.o: nr3td_fi.f95
	${COMP} nr3td_fi.o nr3td_fi.f95

nr3td.o: nr3td.f95 nr3td_fi.o
	${COMP} nr3td.o nr3td.f95


#
# Predefined tasks
#

.PHONY: clean delete

clean: 
	rm -rf *.o *.mod

delete: clean
	rm -rf libaceto*

