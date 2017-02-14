
all: f2py



#f2py:
#	f2py --fcompiler=pg -m nresp3o -c ../lib/nresp3o.f95


f2py:
	f2py -I../lib/ -I../src/ -L../lib -L../src/ -m nr3td_fi -c ../lib/nr3td_fi.f95 -laceto  -lsrcaceto

    
clean:
	rm -rf *.so


re:
	cd ../src; make clean; make
	cd ../lib; make clean; make
	make clean; make

cleanall:
	cd ../src; make clean
	cd ../lib; make clean
	make clean


