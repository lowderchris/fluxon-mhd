use ExtUtils::MakeMaker;

WriteMakefile( DIR=>[ 'Flux' ],
#	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
