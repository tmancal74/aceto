#
#  Aceto library makefile
#
#
#
#

# select suitable predefined variables for compilation
include conf/conf.in
include conf/${COMPILER_SETTINGS}
include conf/pgi_lib.in



#
# 
#

COMP = ${FC} ${FFLAGS}
LINK = ${FC} ${LFLAGS}


main: run_tests

# List of library routines
LIBO = lib/trp2.o

#-----------------------------------------------------------
# Test driver
#-----------------------------------------------------------
run_tests: aceto_test.o
	${LINK} run_tests ${FLIBDIR} aceto_test.o aceto.o ${LIBO} ${FLIB}

aceto_test.o: aceto_test.f03 aceto.o
	${COMP} aceto_test.o  aceto_test.f03

#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f03 ${LIBO}
	${COMP} aceto.o  aceto.f03
 




#
# Predefined tasks
#
run: run_tests
	./run_tests

clean: 
	rm -rf *.o *.mod

delete: clean
	rm -rf run_tests


