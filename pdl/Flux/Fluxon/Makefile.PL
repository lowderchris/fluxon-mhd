use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::Fluxon',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
