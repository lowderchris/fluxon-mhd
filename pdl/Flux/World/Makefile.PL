use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::World',
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
