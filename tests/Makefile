
all: ftrn


ftrn: 
	cd fortran; make

quantarhei: libaceto.so
	python quantarhei_test.py
    
    


libaceto.so:
	cd ../lib; make;
	cp ../lib/libaceto.so .

clean:
	rm -rf *.so *.png 
	cd fortran; make clean

delete: clean
	cd fortran; make delete 
