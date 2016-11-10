use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::World',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
