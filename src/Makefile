# select suitable predefined variables for compilation
include ../conf/conf.in
include ../conf/${COMPILER_SETTINGS}

COMP = ${FC} ${FFLAGS}

#-----------------------------------------------------------
# Library control module
#-----------------------------------------------------------
aceto.o: aceto.f03 
	@echo "Building control module"
	${COMP} aceto.o  aceto.f03
	@echo "...control module built"
	
.PHONY: clean delete

clean:
	rm -rf *.o *.mod 

delete: clean
