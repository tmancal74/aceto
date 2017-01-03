#
#  Aceto library makefile
#
#
#
#

# select suitable predefined variables for compilation
include conf/conf.in
include conf/${COMPILER_SETTINGS}


#
# 
#

COMP = ${FC} ${FFLAGS}
LINK = ${FC} -L./lib/ ${LFLAGS}

main: aceto_test.x


# List of library routines
LIBO = -laceto

#-----------------------------------------------------------
# Test driver
#-----------------------------------------------------------
aceto_test.x: aceto_test.o
	cd lib/; make
	${LINK} aceto_test.x  aceto_test.o aceto.o ${LIBO}

aceto_test.o: aceto_test.f03 aceto.o
	${COMP} aceto_test.o  aceto_test.f03

#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f03 
	${COMP} aceto.o  aceto.f03
 




#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
test: aceto_test.x tests.sh
	./tests.sh

clean:
	cd lib/; make clean
	rm -rf *.o *.mod 

delete: clean
	cd lib/; make delete
	rm -rf aceto_test.x



