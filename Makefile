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

all: library src aceto_test

.PHONY: library src aceto_test

library: 
	cd lib; make 

src:
	cd src; make
	
aceto_test:
	cd test; make





#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
.PHONY: test

test: 
	cd test; make aceto_test.x; ./tests.sh


.PHONY: clean delete

clean:
	cd lib/; make clean
	cd src/; make clean
	cd test/; make clean
	cd aceto; make clean

delete: clean
	cd lib/; make delete
	cd src/; make delete
	cd test/; make delete
	cd aceto/; make clean





