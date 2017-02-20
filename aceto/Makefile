
all: f2py

#f2py:
#	f2py --fcompiler=pg -m nresp3o -c ../lib/nresp3o.f95


f2py:
	f2py -I../lib/ -L../lib -m nr3td_fic -c ../lib/nr3td_fic.f95 -laceto

    
clean:
	rm -rf *.so


re:
	cd ../lib; make clean; make
	make clean; make

cleanall:
	cd ../lib; make clean
	make clean


