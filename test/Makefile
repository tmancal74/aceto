# select suitable predefined variables for compilation
include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

COMP = ${FC} -I../src/ ${FFLAGS} 
LINK = ${FC} -L../lib/ -L../src/ -I../src/ ${LFLAGS}

# List of library routines
LIBO = ../src/aceto.o -laceto

all: aceto_test.x

#-----------------------------------------------------------
# Test driver
#-----------------------------------------------------------

aceto_test.x: aceto_test.o
	@echo "Building test driver ..."
	${LINK} aceto_test.x  aceto_test.o ${LIBO}
	@echo "...test driver built"
	
aceto_test.o: aceto_test.f95
	${COMP} aceto_test.o  aceto_test.f95

#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
test: aceto_test.x tests.sh
	./tests.sh

.PHONY: clean delete

clean:
	rm -rf *.o *.mod 

delete: clean
	rm -rf aceto_test.x