# Top level FLUX makefile -- for now this just branches down into lib
# and makes libflux.a as necessary.
# 
# This would be much better implemented using automake but I'm too lazy
# to figure it out -- CED 18-Aug-2004

all: libinstall
	cd pdl
	perl Makefile.PL
	make all
	cd ..

install: libinstall pdlinstall

libinstall:
	/bin/sh -c 'cd lib; make; make install';

pdlinstall: 
	cd pdl
	perl Makefile.PL
	make all
	cd ..
	cd pdl
	make install
	cd ..

clean:
	rm -f *~ \#* 
	cd lib
	make clean
	cd ..
	cd pdl
	make clean
	cd ..
