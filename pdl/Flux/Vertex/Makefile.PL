use ExtUtils::MakeMaker;


WriteMakefile( NAME => 'Flux::Vertex',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       LIBS=>$libs,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );

