# Top level FLUX makefile -- for now this just branches down into lib
# and makes libflux.a as necessary.
# 
# This would be much better implemented using automake but I'm too lazy
# to figure it out -- CED 18-Aug-2004

all: 
	cd lib; make all

install: libinstall

libinstall:
	cd lib; make install

clean:
	rm -f *~ \#* 
	cd lib; make clean
	cd pdl; make clean
