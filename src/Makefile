# select suitable predefined variables for compilation
include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

COMP = ${FC} -fPIC ${FFLAGS}
LINK = ${CC} ${DYNFLAGS}


LIBO =  aceto.o acetosys.o acetodef.o acetolab.o

all: ${LIBO} libsrcaceto.so

libsrcaceto.so: ${LIBO}
	${LINK} libsrcaceto.so ${LIBO} ${DYNLIBS}



#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f95 acetodef.o
	@echo "Building control module"
	${COMP} aceto.o  aceto.f95
	@echo "... control module built"

#-----------------------------------------------------------
# Library classes 
#-----------------------------------------------------------
acetosys.o: acetosys.f95 acetodef.o
	@echo "Building systems module"
	${COMP} acetosys.o acetosys.f95
	@echo "... systems module built"
	

acetodef.o: acetodef.f95
	${COMP} acetodef.o acetodef.f95

acetolab.o: acetolab.f95
	${COMP} acetolab.o acetolab.f95

.PHONY: clean delete

clean:
	rm -rf *.o *.mod 

delete: clean
