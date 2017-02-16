include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}


COMP = ${FC} -I../src  ${PICFLAGS} ${FFLAGS}  
LINK = ${CC} -L../src  ${DYNFLAGS}

# List of library routines
LIBO = acetosys.o acetodef.o acetolab.o trp2.o nr3td.o nr3td_fi.o

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


nr3td_fi.o: nr3td_fi.f95
	${COMP} nr3td_fi.o nr3td_fi.f95

nr3td.o: nr3td.f95 nr3td_fi.o
	${COMP} nr3td.o nr3td.f95

acetosys.o: acetosys.f95 acetodef.o
	@echo "Building systems module"
	${COMP} acetosys.o acetosys.f95
	@echo "... systems module built"
	

acetodef.o: acetodef.f95
	${COMP} acetodef.o acetodef.f95

acetolab.o: acetolab.f95
	${COMP} acetolab.o acetolab.f95


#
# Predefined tasks
#

.PHONY: clean delete

clean: 
	rm -rf *.mod *.o *.so 

delete: clean
	rm -rf libaceto* 

