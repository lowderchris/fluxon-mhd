use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux',
	       DIR => [qw/ World Fluxon Vertex Concentration /],
#	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
