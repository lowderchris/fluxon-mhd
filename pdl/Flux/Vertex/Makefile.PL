use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::Vertex',
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
