######################################################################
# This file is part of FLUX, the Field Line Universal relaXer.
# Copyright (c) Craig DeForest 2004
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

CFLAGS = -I$(INCLUDE_DIR) -g
LIBS = -lm


HEADERS = $(ID)/data.h $(ID)/io.h $(ID)/geometry.h $(ID)/model.h $(ID)/physics.h

SOURCES = data.c io.c geometry.c model.c physics.c
MODULES = data.o io.o geometry.o model.o physics.o

libflux.a: $(MODULES) $(HEADERS) Makefile
	rm -f libflux.a
	ar ru libflux.a $(MODULES)

install: $(HEADERS) libflux.a
	install libflux.a /usr/local/lib
	install -d /usr/local/include/flux/
	install ../include/*.h /usr/local/include/flux/

lib: libflux.a 
	echo "libflux.a is now up-to-date"

all: 	lib 
	echo made all

clean: 
	rm -f *~ *.a *.o \#* core a.out $(OBJECTS) *.OOGL sample_output.txt

.c.o:	$(HEADERS)
	cc $(CFLAGS) -c $?



