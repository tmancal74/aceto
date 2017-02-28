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

all: library src  aceto_test

.PHONY: library src aceto_test

library: 
	cd lib; make 

src:
	cd src; make
	
aceto_test:
	cd tests; make

install: 
	python setup.py install




#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
.PHONY: test quantarhei_test

test:  
	cd tests; ./tests_quantarhei.sh


quantarhei_test:
	cd tests; python quantarhei_test.py



.PHONY: clean delete

clean:
	rm -rf aceto.egg*
	cd lib/; make clean
	cd src/; make clean
	cd tests/; make clean
	cd aceto; make clean
	rm -rf build

delete: clean
	rm -rf dist
	cd lib/; make delete
	cd src/; make delete
	cd tests/; make delete
	cd aceto/; make clean
     





