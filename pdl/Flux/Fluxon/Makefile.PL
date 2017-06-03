use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::Fluxon',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       LIBS=>$libs,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
