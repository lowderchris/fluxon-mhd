use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::Fluxon',
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
