use ExtUtils::MakeMaker;


WriteMakefile( NAME => 'Flux::Vertex',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );

