use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'FluxSystem',
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
