use ExtUtils::MakeMaker;

WriteMakefile( NAME => 'Flux::World',
	       INC=>$inc,
	       CCFLAGS=>$ccflags,
	       LIBS=>$libs,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT)'
	       );
