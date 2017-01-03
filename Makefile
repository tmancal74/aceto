#
#  Aceto library makefile
#
#
#
#

# select suitable predefined variables for compilation
include gfortran.in


#
# 
#

COMP = ${FC} ${FFLAGS}
LINK = ${FC} ${LFLAGS}

main: run_tests

# List of library routines
LIBO = trp2.o

#-----------------------------------------------------------
# Test driver
#-----------------------------------------------------------
run_tests: aceto_test.o
	${LINK} run_tests  aceto_test.o aceto.o ${LIBO}

aceto_test.o: aceto_test.f03 aceto.o
	${COMP} aceto_test.o  aceto_test.f03

#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f03 ${LIBO}
	${COMP} aceto.o  aceto.f03
 
#----------------------------------------------------------
# Library routines
#----------------------------------------------------------
trp2.o: trp2.f03
	${COMP} trp2.o   trp2.f03



#
# Predefined tasks
#
run: run_tests
	./run_tests

clean: 
	rm -rf *.o *.mod

delete: clean
	rm -rf run_tests


