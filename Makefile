# Top level FLUX makefile -- for now this just branches down into lib
# and makes libflux.a as necessary.
#
# This would be much better implemented using automake but I'm too lazy
# to figure it out -- CED 18-Aug-2004

FL_PREFIX ?= /usr/local
export FL_PREFIX

everything: libbuild install

install: libinstall pdlbuild pdltest pdlinstall

libbuild:
	/bin/sh -c 'cd lib; FL_PREFIX=$(FL_PREFIX) make';

libinstall: libbuild
	/bin/sh -c 'cd lib; FL_PREFIX=$(FL_PREFIX) make install';

pdlbuild:
	@echo "Using $(PL_PREFIX) in INSTALL_BASE"; \
	cd pdl ; \
	perl Makefile.PL INSTALL_BASE=$(PL_PREFIX); \
	PERL_INSTALL_QUIET=1 make all ; \
	cd .. ;

pdltest:
	{ cd pdl; } && { make test; } && { cd ..; }

pdlinstall:
	cd pdl ; \
	make install ; \
	cd .. ;

clean:
	rm -f *~ \#* ; \
	cd lib; \
	make clean; \
	cd ..; \
	cd pdl; \
	make clean; \
	cd .. ;

realclean: clean
	rm -rf sandbox

uninstall:
	rm -r $(FL_PREFIX)/lib/libflux.a $(FL_PREFIX)/include/flux ;
