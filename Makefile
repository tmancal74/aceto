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

all: library src aceto_test

.PHONY: library src aceto_test

library: 
	cd lib; make 

src:
	cd src; make
	
aceto_test:
	cd tests; make

egg:
	python setup.py bdist_egg

install: library src aceto_test 
	python setup.py install
	python postinstall.py


post: 
	python postinstall.py



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




