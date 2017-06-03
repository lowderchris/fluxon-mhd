use ExtUtils::MakeMaker;


WriteMakefile( NAME => 'Flux::Concentration',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       LIBS=>$libs,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
