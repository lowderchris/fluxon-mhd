INCLUDE_DIR = ../include
ID = $(INCLUDE_DIR)

CFLAGS = -I$(INCLUDE_DIR) -g
LIBS = -lm


HEADERS = $(ID)/data.h $(ID)/io.h $(ID)/geometry.h $(ID)/model.h $(ID)/geomview.h $(ID)/hvectext.h $(ID)/physics.h

SOURCES = data.c io.c geometry.c model.c geomview.c hvectext.c physics.c
MODULES = data.o io.o geometry.o model.o geomview.o hvectext.o physics.o

libflux.a: $(MODULES)	
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



