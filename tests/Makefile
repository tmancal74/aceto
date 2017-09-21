LNAME = aceto

all: ftrn


ftrn: 
	cd fortran; make LNAME=${LNAME}

quantarhei: ${LNAME}
	python quantarhei_test.py
    
    


${LNAME}:
	cd ../lib; make;
	cp ../lib/lib${LNAME}.so .

clean:
	rm -rf *.so *.png 
	cd fortran; make clean

delete: clean
	cd fortran; make delete 
