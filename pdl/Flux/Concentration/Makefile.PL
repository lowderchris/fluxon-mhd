use ExtUtils::MakeMaker;


WriteMakefile( NAME => 'Flux::Concentration',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
