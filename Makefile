# Top level FLUX makefile -- for now this just branches down into lib
# and makes libflux.a as necessary.
# 
# This would be much better implemented using automake but I'm too lazy
# to figure it out -- CED 18-Aug-2004

all: libinstall
	pushd pdl; perl Makefile.PL; make all; popd

install: libinstall pdlinstall

libinstall:
	pushd lib; make install; popd;

pdlinstall:
	pushd pdl; perl Makefile.PL; make all; popd
	pushd pdl/PDL; make install; popd

clean:
	rm -f *~ \#* 
	cd lib; make clean
	cd pdl; make clean
