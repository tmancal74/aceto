FC = gfortran
FFLAGS = -c -W -Wall -o
LFLAGS = -o

COMP = ${FC} -I../../src/ -I../../src/lib/src ${FFLAGS} 
LINK = ${FC} -L../../src/ -L../../src/lib/src ${LFLAGS}

# List of library routines
LLINK = -lfaceto -laceto



all: aceto_test.x  band_system_test.x 

#-----------------------------------------------------------
# Test drivers
#-----------------------------------------------------------

aceto_test.x: aceto_test.o
	@echo "Building test driver ..."
	${LINK} aceto_test.x  aceto_test.o ${LLINK}
	@echo "...test driver built"

aceto_test.o: aceto_test.f95
	${COMP} aceto_test.o  aceto_test.f95

band_system_test.x: band_system_test.o
	${LINK} band_system_test.x band_system_test.o ${LLINK}

band_system_test.o: band_system_test.f90
	${COMP} band_system_test.o band_system_test.f90


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
	rm -rf *.o *.mod *.so *.a

delete: clean
	rm -rf *.x
