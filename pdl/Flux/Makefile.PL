use ExtUtils::MakeMaker;


WriteMakefile( NAME => 'Flux',
	       DIR => [qw/ Core World Fluxon Vertex Concentration /],
	       INC => join(" ",map { "-I$_"} @inc),
	       LIBS=>$libs,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
