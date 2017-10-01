# select suitable predefined variables for compilation
include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

COMP = ${FC} -fPIC  -I../lib ${FFLAGS}
LINK = ${CC} -L../lib ${DYNFLAGS} -laceto

all: aceto.o 


#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f95 
	@echo "Building control module"
	${COMP} aceto.o  aceto.f95
	@echo "... control module built"


.PHONY: clean delete

clean:
	rm -rf *.o *.mod 

delete: clean
