# Top level FLUX makefile -- for now this just branches down into lib
# and makes libflux.a as necessary.
#
# This would be much better implemented using automake but I'm too lazy
# to figure it out -- CED 18-Aug-2004

# FL_PREFIX ?= $(FL_PREFIX)
# export FL_PREFIX

everything: libbuild install done

done:
	@echo "\n\nFlux is now installed!\n\n";

install: libinstall pdlbuild pdltest pdlinstall

paths:
	@/bin/sh -c 'fluxpype/PREFIX_PATHS.sh';

libbuild:
	/bin/sh -c 'cd lib; FL_PREFIX=$(FL_PREFIX) make';

libinstall: libbuild
	/bin/sh -c 'cd lib; FL_PREFIX=$(FL_PREFIX) make install';

pdlMakefile:
	/bin/sh -c 'if [ -n "$(PL_PREFIX)" ]; then echo "\n\nWill install Flux perl modules into $(PL_PREFIX). Make sure this is in your @INC. See README.\n\n"; else echo "INFO: PL_PREFIX is not defined (this is probably OK)."; fi'; \
	/bin/sh -c 'cd pdl; if [ ! -f Makefile ]; then perl Makefile.PL INSTALL_BASE=$(PL_PREFIX); fi';

pdlbuild: pdlMakefile
	/bin/sh -c 'cd pdl; if !(PERL_INSTALL_QUIET=1 make all); then (echo "\n\nRunning make for you. You are welcome!\n"; PERL_INSTALL_QUIET=1 make all); fi';

pdltest:
	/bin/sh -c 'cd pdl; make test ;';

pdlinstall:
	/bin/sh -c 'cd pdl; make install;';

realclean: clean
	rm -rf sandbox

clean:
	rm -f *~ \#* ; \
	rm -f pdl/*.png; \
	cd lib; \
	make clean; \
	cd ..; \
	cd pdl; \
	make clean; \
	cd .. ; \


uninstall:
	@echo "\nUninstalling FLUX...";
	@-rm -rf $(FL_PREFIX)/lib/libflux.a || true;
	@-rm -rf $(FL_PREFIX)/include/flux || true;
	@echo "\tFlux uninstall complete.\n";
