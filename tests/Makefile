
all: quantarhei


quantarhei: libaceto.so
	python quantarhei_test.py
    
    


libaceto.so:
	cd ../lib; make;
	cp ../lib/libaceto.so .

clean:
	rm -rf *.so 
	cd fortran; make clean
    
