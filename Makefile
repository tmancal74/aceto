#
#  Aceto library makefile
#
#
#
#

#
# select suitable predefined variables for compilation
#
include conf/conf.in
include conf/${COMPILER_SETTINGS}


#
# compiler flags are read from conf.in file 
#
COMP = ${FC} ${FFLAGS}
LINK = ${FC} -L./lib/ ${LFLAGS}
LIBNAME = "aceto-0.0.4-darwin"

all: library src aceto_test

.PHONY: library src aceto_test

library: 
	cd lib; make LNAME=${LIBNAME}

src:
	cd src; make
	
aceto_test:
	cd tests; make LNAME=${LIBNAME}

egg: library
	python setup.py bdist_egg

install: library src aceto_test 
	python setup.py install
	aceto_conf


post: 
	aceto_conf
	



#-----------------------------------------------------------
# Predefined tasks
#-----------------------------------------------------------
.PHONY: test quantarhei_test

test:  
	cd tests; qrhei test.py


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

     
uninst: 
	pip uninstall -y aceto
	
reinst: uninst install
	@echo "Reinstallation complete"




