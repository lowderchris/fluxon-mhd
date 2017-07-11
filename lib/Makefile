######################################################################
# This file is part of FLUX, the Field Line Universal relaXer.
# Copyright (c) Craig DeForest 2004-2007
#
# You may modify and/or distribute this software under the terms of
# the Gnu Public License, version 2.  You should have received a copy
# of the license with this software, in the file COPYING to be found
# in the top level directory of the distribution.  You may obtain
# additional copies of the license via http://www.gnu.org or by 
# writing to the Free Software Foundation, 59 Temple Place - Suite 330
# Boston, MA 02111-1307 USA.
#

INCLUDE_DIR = ../include
ID = $(INCLUDE_DIR)
FL_PREFIX ?= /usr/local
SANDBOX = ../sandbox

# DFLAGS = -DUSE_DEBUGGING_MALLOC
# DFLAGS = -I/usr/lib/perl5/5.8.5/i686-linux/CORE -DUSE_PERL_MALLOC
# DFLAGS = -DUSE_PADDED_MALLOC

CC = cc

#CFLAGS = -I$(INCLUDE_DIR) -g -fPIC -O0 --coverage
CFLAGS = -I$(INCLUDE_DIR) -g -fPIC
#CFLAGS = -I$(INCLUDE_DIR) -O3

#LIBS not actually used anywhere
#LIBS = -lm -lgcov --coverage


HEADERS = $(ID)/data.h $(ID)/io.h $(ID)/geometry.h $(ID)/model.h $(ID)/physics.h $(ID)/fluxperl.h

SOURCES = data.c io.c geometry.c model.c physics.c
MODULES = data.o io.o geometry.o model.o physics.o

libflux.a: $(MODULES) $(HEADERS) Makefile
	rm -f libflux.a
	ar ru libflux.a $(MODULES)
	ranlib libflux.a
	install -d $(SANDBOX)/lib $(SANDBOX)/include/flux
	install libflux.a $(SANDBOX)/lib
	install $(ID)/*.h $(SANDBOX)/include/flux

install: $(HEADERS) libflux.a
	install -d $(FL_PREFIX)/lib
	install libflux.a $(FL_PREFIX)/lib
#	ranlib $(FL_PREFIX)/lib/libflux.a
	install -d $(FL_PREFIX)/include/flux/
	install $(ID)/*.h $(FL_PREFIX)/include/flux/

lib: libflux.a
	echo "libflux.a is now up-to-date"

all: 	lib
	echo made all

clean:
	rm -f *~ *.a *.o \#* core a.out $(OBJECTS) *.OOGL sample_output.txt *.cf.c

%.o:	%.c $(HEADERS) Makefile
	rm -f $<f.c
	./subst_filename $< $<f.c
	$(CC) $(CFLAGS) -o $@ -c $(DFLAGS) $<f.c
