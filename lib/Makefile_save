include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

LNAME = aceto # this is default, LNAME is defined on command line
LFILE = lib${LNAME}.a
COMP = ${FC} -I../src  ${PICFLAGS} ${FFLAGS}  
LINK = ${CC} -L../src  ${DYNFLAGS}

AR = ar -r

# List of library routines
LIBO = acetosys.o acetodef.o acetolab.o acetoaux.o trp2.o nr3td.o nr3td_fic.o

libaceto.a: echo_building_library ${LIBO}
	${AR} ${LFILE} ${LIBO} ${DYNLIBS}
	@echo "...library built"

echo_building_library:
	@echo "Building library ..."

#----------------------------------------------------------
# Library routines
#----------------------------------------------------------
trp2.o: trp2.f95
	${COMP} trp2.o  trp2.f95


nr3td_fic.o: nr3td_fic.f90 acetoaux.o
	${COMP} nr3td_fic.o nr3td_fic.f90

nr3td_int.o: nr3td_int.f95
	${COMP} nr3td_int.o nr3td_int.f95
     
nr3td.o: nr3td.f95 nr3td_int.o nr3td_fic.o
	${COMP} nr3td.o nr3td.f95

acetosys.o: acetosys.f95 acetodef.o
	@echo "Building systems module"
	${COMP} acetosys.o acetosys.f95
	@echo "... systems module built"
	

acetodef.o: acetodef.f95
	${COMP} acetodef.o acetodef.f95

acetolab.o: acetolab.f95
	${COMP} acetolab.o acetolab.f95

acetoaux.o: acetoaux.f95
	${COMP} acetoaux.o acetoaux.f95


#
# Predefined tasks
#

.PHONY: clean delete

clean: 
	rm -rf *.mod *.o *.so *.a

delete: clean
	rm -rf libaceto* 

