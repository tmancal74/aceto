# select suitable predefined variables for compilation
include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

COMP = ${FC} -I../src/ -I../lib ${FFLAGS} 
LINK = ${FC} -L../lib/ -L../src/ -I../src/ ${LFLAGS}

# List of library routines
LIBO = ../src/aceto.o ../src/acetosys.o ../src/acetolab.o -laceto

all: aceto_test.x band_system_test.x

#-----------------------------------------------------------
# Test drivers
#-----------------------------------------------------------

aceto_test.x: aceto_test.o
	@echo "Building test driver ..."
	${LINK} aceto_test.x  aceto_test.o ${LIBO}
	@echo "...test driver built"
	
aceto_test.o: aceto_test.f95
	${COMP} aceto_test.o  aceto_test.f95

band_system_test.x: band_system_test.o
	${LINK} band_system_test.x band_system_test.o ${LIBO}

band_system_test.o: band_system_test.f95
	${COMP} band_system_test.o band_system_test.f95


#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
test: test_aceto test_bs


test_aceto: aceto_test.x tests_aceto.sh
	./tests_aceto.sh

test_bs: band_system_test.x tests_bs.sh
	./tests_bs.sh


.PHONY: clean delete

clean:
	rm -rf *.o *.mod 

delete: clean
	rm -rf *.x
